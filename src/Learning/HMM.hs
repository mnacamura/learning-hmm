module Learning.HMM (
    HMM (..)
  , LogLikelihood
  , new
  , viterbi
  , baumWelch
  ) where

import Control.Arrow ((***), first)
import Data.Random.Distribution (pdf)
import Data.Random.Distribution.Categorical (Categorical)
import qualified Data.Random.Distribution.Categorical as C (
    fromList, fromWeightedList, normalizeCategoricalPs
  )
import Data.Random.Distribution.Categorical.Util ()
import Data.List (genericLength)
import Data.Number.LogFloat (fromLogFloat, logFloat, logFromLogFloat)
import Data.Vector ((!))
import qualified Data.Vector as V (elemIndex, fromList, map, toList, zip)
import qualified Data.Vector.Util.LinearAlgebra as V (transpose)
import Learning.HMM.Internal

type LogLikelihood = Double

-- | Parameter set of the hidden Markov model. Direct use of the
--   constructor is not recommended. Instead, call 'new'.
data HMM s o = HMM { states  :: [s] -- ^ Hidden states
                   , outputs :: [o] -- ^ Observed outputs
                   , initialStateDist :: Categorical Double s
                     -- ^ Categorical distribusion of initial states
                   , transitionDist :: s -> Categorical Double s
                     -- ^ Categorical distribution of next states
                     --   conditioned by the previous states
                   , emissionDist :: s -> Categorical Double o
                     -- ^ Categorical distribution of outputs conditioned
                     --   by the hidden states
                   }

instance (Show s, Show o) => Show (HMM s o) where
  show = showHMM

showHMM :: (Show s, Show o) => HMM s o -> String
showHMM hmm = "HMM {states = "           ++ show ss
              ++ ", outputs = "          ++ show os
              ++ ", initialStateDist = " ++ show pi0
              ++ ", transitionDist = "   ++ show [(w s, s) | s <- ss]
              ++ ", emissionDist = "     ++ show [(phi s, s) | s <- ss]
              ++ "}"
  where
    ss  = states hmm
    os  = outputs hmm
    pi0 = initialStateDist hmm
    w   = transitionDist hmm
    phi = emissionDist hmm

-- | Construct a 'HMM' from the given states and outputs. The
--   'initialStateDist' and 'emissionDist' are set to be uniform
--   distributions. The 'transitionDist' is specified as follows: with
--   probability 1/2, move to the same state, otherwise, move to a random
--   state (which might be the same state).
new :: (Ord s, Ord o) => [s] -> [o] -> HMM s o
new ss os = HMM { states           = ss
                , outputs          = os
                , initialStateDist = pi0
                , transitionDist   = w
                , emissionDist     = phi
                }
  where
    pi0 = C.fromWeightedList [(1, s) | s <- ss]
    w s | s `notElem` ss = C.fromList []
        | otherwise      = C.fromList [(p s', s') | s' <- ss]
      where
        k = genericLength ss
        p s' | s' == s   = 1/2 * (1 + 1/k)
             | otherwise = 1/2 / k
    phi s | s `notElem` ss = C.fromList []
          | otherwise      = C.fromWeightedList [(1, o) | o <- os]

-- | Perform the Viterbi algorithm and return the most likely state path
--   and its log likelihood.
viterbi :: (Eq s, Eq o) => HMM s o -> [o] -> ([s], LogLikelihood)
viterbi model xs =
  checkModelIn "viterbi" model `seq`
  checkDataIn "viterbi" model xs `seq`
  (V.toList *** logFromLogFloat) $ viterbi' model' xs'
  where
    model' = toHMM' model
    xs'    = V.fromList xs

-- | Perform the Baum-Welch EM algorithm steps iteratively and return
--   a list of updated 'HMM' parameters and their corresponding log
--   likelihoods.
baumWelch :: (Eq s, Eq o) => HMM s o -> [o] -> [(HMM s o, LogLikelihood)]
baumWelch model xs =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs `seq`
  map (fromHMM' *** logFromLogFloat) $ baumWelch' model' xs'
  where
    model' = toHMM' model
    xs'    = V.fromList xs

-- | Check if the 'HMM' is valid in the sense of whether the 'states' and
--   'outputs' are not empty.
checkModelIn :: String -> HMM s o -> ()
checkModelIn fun hmm
  | null ss   = err "empty states"
  | null os   = err "empty outputs"
  | otherwise = ()
  where
    ss = states hmm
    os = outputs hmm
    err = errorIn fun

-- | Check if all the elements of the data are contained in the 'outputs'
--   of the 'HMM'.
checkDataIn :: Eq o => String -> HMM s o -> [o] -> ()
checkDataIn fun hmm xs
  | all (`elem` os) xs = ()
  | otherwise          = err "illegal data"
  where
    os = outputs hmm
    err = errorIn fun

-- | Convert 'HMM'' to 'HMM'.
fromHMM' :: (Eq s, Eq o) => HMM' s o -> HMM s o
fromHMM' hmm' = HMM { states           = V.toList ss
                    , outputs          = V.toList os
                    , initialStateDist = C.fromList pi0'
                    , transitionDist   = \s -> case V.elemIndex s ss of
                                                 Nothing -> C.fromList []
                                                 Just i  -> C.fromList $ w' i
                    , emissionDist     = \s -> case V.elemIndex s ss of
                                                 Nothing -> C.fromList []
                                                 Just i  -> C.fromList $ phi' i
                    }
  where
    ss  = states' hmm'
    os  = outputs' hmm'
    pi0 = initialStateDist' hmm'
    w   = transitionDist' hmm'
    phi = V.transpose $ emissionDistT' hmm'
    pi0'   = V.toList $ V.map (first fromLogFloat) $ V.zip pi0 ss
    w' i   = V.toList $ V.map (first fromLogFloat) $ V.zip (w ! i) ss
    phi' i = V.toList $ V.map (first fromLogFloat) $ V.zip (phi ! i) os

-- | Convert 'HMM' to 'HMM''. The 'initialStateDist'', 'transisionDist'',
--   and 'emissionDistT'' are normalized.
toHMM' :: (Eq s, Eq o) => HMM s o -> HMM' s o
toHMM' hmm = HMM' { states'           = V.fromList ss
                  , outputs'          = V.fromList os
                  , initialStateDist' = V.fromList pi0'
                  , transitionDist'   = V.fromList w'
                  , emissionDistT'    = V.fromList phi'
                  }
  where
    ss  = states hmm
    os  = outputs hmm
    pi0 = C.normalizeCategoricalPs $ initialStateDist hmm
    w   = C.normalizeCategoricalPs . transitionDist hmm
    phi = C.normalizeCategoricalPs . emissionDist hmm
    pi0' = [logFloat $ pdf pi0 s | s <- ss]
    w'   = [V.fromList [logFloat $ pdf (w s) s' | s' <- ss] | s <- ss]
    phi' = [V.fromList [logFloat $ pdf (phi s) o | s <- ss] | o <- os]

errorIn :: String -> String -> a
errorIn fun msg = error $ "Learning.HMM." ++ fun ++ ": " ++ msg
