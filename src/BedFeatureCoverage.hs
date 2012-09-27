{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main
  where

import Control.Applicative
import Control.Arrow
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as LBS
import Data.IORef
import Data.List
import Data.Maybe
import Numeric
import System.Console.GetOpt
import System.IO

import qualified Data.Iteratee.IO as IterIO
import qualified Data.Iteratee.ListLike as IterLL
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Statistics.Function as SF
import Statistics.Quantile
import Statistics.Sample
import System.Console.CmdTheLine

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.BamIndex as BamIdx
import qualified Bio.SamTools.Iteratee as BamIter
import qualified Bio.SeqLoc.Bed as Bed
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

main :: IO ()
main = run ( termCoverage, info)
  where info = defTI { termName = "bed-feature-coverage"
                     , version = "0.0"
                     , termDoc = "Compute per-feature coverage statistics for a BED file"
                     }
        termCoverage = coverage <$> bamfile <*> bedfile <*> output
        coverage bam bed out = BamIdx.withIndex bam $ \hidx -> 
          withFile out WriteMode $ \hout ->
          processBed hidx hout bed

processBed :: BamIdx.IdxHandle -> Handle -> FilePath -> IO ()
processBed hidx hout = IterIO.fileDriver trxIter 
  where trxIter = Bed.bedTranscriptEnum $ IterLL.mapM_ trxOne
        trxOne trx = trxLine hidx trx >>= BS.hPutStr hout

trxLine :: BamIdx.IdxHandle -> Transcript -> IO (BS.ByteString)
trxLine hidx trx = liftM (maybe "" profLine) $ countTrxReads hidx trx
  where profLine prof = let fields = statFields prof
                        in BS.unlines . (: []) $
                           BS.intercalate "\t" $
                           (unSeqLabel . geneId $ trx) : fields
                           
statFields :: U.Vector Int -> [BS.ByteString]
statFields zprof = map BS.pack [ show $ U.length zprof
                               , show total
                               , showFFloat (Just 2) (mean xprof) ""
                               , showFFloat (Just 2) med ""
                               , show npos
                               , showFFloat (Just 2) (avgpos) ""
                               , show n50
                               , showFFloat (Just 4) f50 ""
                               , showFFloat (Just 4) lorenz ""
                               ]
  where xprof = U.map fromIntegral zprof
        total = U.sum zprof -- Total mass
        med = continuousBy medianUnbiased 1 2 xprof
        npos = U.length $ U.filter (> 0) zprof
        avgpos = (fromIntegral total) / (fromIntegral npos)
        xsrt = U.reverse $ SF.sort xprof
        cprof = U.map (* (recip . fromIntegral $ total)) . U.postscanl' (+) 0.0 $ xsrt
        n50 = fromMaybe 0 $ U.findIndex (>= 0.50) cprof
        (f50 :: Double) = (fromIntegral n50) / (fromIntegral $ U.length xprof)
        lorenz = let unif = U.generate (U.length cprof) unifAt
                     unifAt i = (fromIntegral $ i + 1) / (fromIntegral $ U.length cprof)
                 in (U.sum $ U.zipWith (-) cprof unif) / (fromIntegral $ U.length cprof)
        
countTrxReads :: BamIdx.IdxHandle -> Transcript -> IO (Maybe (U.Vector Int))
countTrxReads hidx trx = maybe (return Nothing) (liftM Just . withTarget) $ 
                         Bam.lookupTarget header $ unSeqLabel trxref
  where header = BamIdx.idxHeader hidx
        (OnSeq trxref trxsploc) = location trx
        trxbnds = (fromIntegral *** fromIntegral) $ Loc.bounds trxsploc
        withTarget tid = do prof <- UM.replicate (fromIntegral $ Loc.length $ trxsploc) 0
                            let countReads = IterLL.mapM_ $ countRead prof trxsploc
                                targetIter = IterLL.joinIM $ BamIter.enumIndexRegion hidx tid trxbnds countReads
                            IterLL.run targetIter
                            U.freeze prof
                     
countRead :: (MonadIO m) => UM.IOVector Int -> SpLoc.SpliceLoc -> Bam.Bam1 -> m ()
countRead prof trxsploc b = case Bam.refSpLoc b >>= liftM Loc.startPos . flip SpLoc.locInto trxsploc of
  Just (Pos.Pos off Plus) -> countOne prof $ fromIntegral off
  _ -> return ()

countOne :: (MonadIO m) => UM.IOVector Int -> Int -> m ()
countOne prof off | off >= 0 && off < UM.length prof = liftIO $ do n0 <- UM.read prof off
                                                                   UM.write prof off $! succ n0
                  | otherwise = return ()

output :: Term String
output = required $ opt Nothing $ (optInfo [ "o", "output" ])
  { optName = "OUT", optDoc = "Output filename" }
         
bedfile :: Term String
bedfile = required $ opt Nothing $ (optInfo [ "b", "bed" ])
  { optName = "BED", optDoc = "Bed-format annotation filename" }

bamfile :: Term String
bamfile = required $ pos 0 Nothing $ posInfo { posName = "BAM", posDoc = "BAM-format alignments" }