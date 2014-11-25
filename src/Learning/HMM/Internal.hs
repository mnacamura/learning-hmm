module Learning.HMM.Internal (
    HMM' (..)
  , Likelihood
  , Probability
  , init'
  , withEmission'
  , viterbi'
  , baumWelch'
  -- , baumWelch1'
  -- , forward'
  -- , backward'
  ) where

import Control.Applicative ((<$>))
import Control.Monad (forM_, replicateM)
import Control.Monad.ST (runST)
import qualified Data.Map.Strict as M (findWithDefault)
import Data.Number.LogFloat (LogFloat, logFloat)
import Data.Random.RVar (RVar)
import Data.Random.Distribution.Simplex (stdSimplex)
import Data.Random.Distribution.Uniform.Util ()
import Data.Vector (Vector, (!))
import qualified Data.Vector as V (
    filter, foldl1', freeze, fromList, generate, last, length, map, maximum
  , maxIndex, replicate, sum, tail, zip, zipWith, zipWith3, zipWith4
  )
import qualified Data.Vector.Mutable as MV (new, read, write)
import qualified Data.Vector.Util as V (frequencies)
import Data.Vector.Util.LinearAlgebra (
    (>+>), (>.>), (>/>), (#+#), (.>), (>/), (#/), (<.>), (#.>), (<.#)
  )
import qualified Data.Vector.Util.LinearAlgebra as V (transpose)

type Likelihood  = LogFloat
type Probability = LogFloat

-- | More efficient data structure of the 'HMM' model. The 'states' and
--   'outputs' in 'HMM' are represented by their indices. The
--   'initialStateDist', 'transitionDist', and 'emissionDist' are
--   represented by matrices. The 'emissionDistT'' is a transposed matrix
--   in order to simplify the calculation.
data HMM' = HMM' { nStates'          :: Int -- ^ Number of states
                 , nOutputs'         :: Int -- ^ Number of outputs
                 , initialStateDist' :: Vector Probability
                 , transitionDist'   :: Vector (Vector Probability)
                 , emissionDistT'    :: Vector (Vector Probability)
                 }

init' :: Int -> Int -> RVar HMM'
init' n m = do
  pi0 <- V.fromList <$> stdSimplex (n-1)
  w   <- V.fromList <$> replicateM n (V.fromList <$> stdSimplex (n-1))
  phi <- V.fromList <$> replicateM n (V.fromList <$> stdSimplex (m-1))
  return HMM' { nStates'          = n
              , nOutputs'         = m
              , initialStateDist' = pi0
              , transitionDist'   = w
              , emissionDistT'    = V.transpose phi
              }

withEmission' :: HMM' -> Vector Int -> HMM'
withEmission' model xs = model { emissionDistT' = phi' }
  where
    ss   = V.generate (nStates' model) id
    os   = V.generate (nOutputs' model) id
    phi' = let (path, _) = viterbi' model xs
               freqs     = V.frequencies $ V.zip path xs
               hists     = V.map (\s -> V.map (\o ->
                                 M.findWithDefault 0 (s, o) freqs) os) ss
           in V.transpose $ V.map (\f -> f >/ V.sum f) hists

viterbi' :: HMM' -> Vector Int -> (Vector Int, Likelihood)
viterbi' model xs = (path, likelihood)
  where
    -- The following procedure is based on
    -- http://ibisforest.org/index.php?cmd=read&page=Viterbi%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0&word=Viterbi
    path = runST $ do
      ix <- MV.new n
      ix `MV.write` (n-1) $ V.maxIndex $ deltas ! (n-1)
      forM_ (reverse [1..(n-1)]) $ \t -> do
        i <- ix `MV.read` t
        ix `MV.write` (t-1) $ psis ! t ! i
      V.freeze ix

    likelihood = V.maximum $ deltas ! (n-1)

    deltas :: Vector (Vector Probability)
    psis   :: Vector (Vector Int)
    (deltas, psis) = runST $ do
      ds <- MV.new n
      ps <- MV.new n
      ds `MV.write` 0 $ (phi' ! (xs ! 0)) >.> pi0
      ps `MV.write` 0 $ V.replicate k (0 :: Int)
      forM_ [1..(n-1)] $ \t -> do
        d <- ds `MV.read` (t-1)
        let dws = V.map (d >.>) w'
        ds `MV.write` t $ (phi' ! (xs ! t)) >.> V.map V.maximum dws
        ps `MV.write` t $ V.map V.maxIndex dws
      ds' <- V.freeze ds
      ps' <- V.freeze ps
      return (ds', ps')
      where
        k    = nStates' model
        pi0  = initialStateDist' model
        w'   = V.transpose $ transitionDist' model
        phi' = emissionDistT' model

    -- Here we assumed that
    n = V.length xs

baumWelch' :: HMM' -> Vector Int -> [(HMM', Likelihood)]
baumWelch' model xs = zip ms $ tail ells
  where
    n = V.length xs
    (ms, ells) = unzip $ iterate step (model, undefined)
    step (m, _) = m `seq` baumWelch1' m n xs

-- | Perform one step of the Baum-Welch algorithm and return the updated
--   model and the likelihood of the old model.
baumWelch1' :: HMM' -> Int -> Vector Int -> (HMM', Likelihood)
baumWelch1' model n xs = (model', likelihood)
  where
    model' = model { initialStateDist' = pi0
                   , transitionDist'   = w
                   , emissionDistT'    = phi'
                   }
    likelihood = V.last ells

    -- First, we calculate the alpha and beta values using the
    -- forward-backward algorithm.
    alphas = forward' model n xs
    betas  = backward' model n xs

    -- Then, we obtain the likelihoods for each time. This should be
    -- constant over time.
    ells = V.zipWith (<.>) alphas betas

    -- Based on the alpha, beta, and likelihood values, we calculate the
    -- gamma and xi values.
    gammas = V.zipWith3 (\a b l -> a >.> b >/ l) alphas betas ells
    xis    = V.zipWith4 (\a b l x -> let w1 = V.zipWith (.>) a w0
                                         w2 = V.map (phi0 ! x >.> b >.>) w1
                                     in w2 #/ l)
               alphas (V.tail betas) (V.tail ells) (V.tail xs)
      where
        w0   = transitionDist' model
        phi0 = emissionDistT' model

    -- Using the gamma and xi values, we finally obtain the optimal initial
    -- state probability vector, transition probability matrix, and
    -- emission probability matrix.
    pi0  = let gs = gammas ! 0
           in gs >/ V.sum gs
    w    = let ws = V.foldl1' (#+#) xis
               zs = V.map V.sum ws
           in V.zipWith (>/) ws zs
    phi' = let gs' o = V.map snd $ V.filter ((== o) . fst) $ V.zip xs gammas
               phis  = V.foldl1' (>+>) . gs'
               zs    = V.foldl1' (>+>) gammas
           in V.map (\o -> phis o >/> zs) os

    -- Here we assumed that
    os = V.generate (nOutputs' model) id

forward' :: HMM' -> Int -> Vector Int -> Vector (Vector Probability)
{-# INLINE forward' #-}
forward' model n xs = runST $ do
  v <- MV.new n
  v `MV.write` 0 $ (phi' ! (xs ! 0)) >.> pi0
  forM_ [1..(n-1)] $ \t -> do
    a <- v `MV.read` (t-1)
    v `MV.write` t $ (phi' ! (xs ! t)) >.> (a <.# w)
  V.freeze v
  where
    pi0  = initialStateDist' model
    w    = transitionDist' model
    phi' = emissionDistT' model

backward' :: HMM' -> Int -> Vector Int -> Vector (Vector Probability)
{-# INLINE backward' #-}
backward' model n xs = runST $ do
  v <- MV.new n
  v `MV.write` (n-1) $ V.replicate k $ logFloat (1 :: Double)
  forM_ (reverse [1..(n-1)]) $ \t -> do
    b <- v `MV.read` t
    v `MV.write` (t-1) $ w #.> ((phi' ! (xs ! t)) >.> b)
  V.freeze v
  where
    k    = nStates' model
    w    = transitionDist' model
    phi' = emissionDistT' model
