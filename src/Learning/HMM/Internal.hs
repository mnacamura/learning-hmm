module Learning.HMM.Internal (
    HMM' (..)
  , Likelihood
  , Probability
  , viterbi'
  , baumWelch'
  -- , baumWelch1'
  -- , forward'
  -- , backward'
  ) where

import Control.Monad (forM_)
import Control.Monad.ST (runST)
import Data.Number.LogFloat (LogFloat, logFloat)
import Data.Vector (Vector, (!))
import qualified Data.Vector as V (
    filter, foldl1', freeze, last, length, map, maximum, maxIndex
  , replicate, sum , tail, zip , zipWith, zipWith3, zipWith4
  )
import qualified Data.Vector.Mutable as M (new, read, write)
import qualified Data.Vector.Util as V (unsafeElemIndex)
import Data.Vector.Util.LinearAlgebra (
    (>+>), (>.>), (>/>), (#+#), (.>), (>/), (#/), (<.>), (#.>), (<.#)
  )
import qualified Data.Vector.Util.LinearAlgebra as V (transpose)

type Likelihood  = LogFloat
type Probability = LogFloat

-- | More efficient data structure of the HMM parameters. This should be
--   only used internally. The 'emissionDistT'' is a transposed matrix in
--   order to simplify the calculation.
data HMM' s o = HMM' { states'           :: Vector s
                     , outputs'          :: Vector o
                     , initialStateDist' :: Vector Probability
                     , transitionDist'   :: Vector (Vector Probability)
                     , emissionDistT'    :: Vector (Vector Probability)
                     }

-- | Perform the Viterbi algorithm and return the most likely state path
--   and its likelihood.
viterbi' :: Eq o => HMM' s o -> Vector o -> (Vector s, Likelihood)
viterbi' model xs = (path, likelihood)
  where
    n = V.length xs

    path = V.map (ss !) $ runST $ do
      ix <- M.new n
      ix `M.write` (n-1) $ V.maxIndex $ deltas ! (n-1)
      forM_ (reverse [0..(n-2)]) $ \i -> do
        j <- ix `M.read` (i+1)
        ix `M.write` i $ psis ! (i+1) ! j
      V.freeze ix
      where
        ss = states' model
    likelihood = V.maximum $ deltas ! (n-1)

    deltas :: Vector (Vector Probability)
    psis   :: Vector (Vector Int)
    (deltas, psis) = runST $ do
      ds <- M.new n
      ps <- M.new n
      ds `M.write` 0 $ (phi' ! x 0) >.> pi0
      ps `M.write` 0 $ V.replicate k (0 :: Int)
      forM_ [1..(n-1)] $ \i -> do
        d <- ds `M.read` (i-1)
        let dws = V.map (d >.>) w'
        ds `M.write` i $ phi' ! x i >.> V.map V.maximum dws
        ps `M.write` i $ V.map V.maxIndex dws
      ds' <- V.freeze ds
      ps' <- V.freeze ps
      return (ds', ps')
      where
        k   = V.length $ states' model
        x i = let os  = outputs' model
                  xs' = V.map (`V.unsafeElemIndex` os) xs
              in xs' ! i
        pi0  = initialStateDist' model
        w'   = V.transpose $ transitionDist' model
        phi' = emissionDistT' model

-- | Perform the Baum-Welch EM algorithm steps iteratively and return
--   a list of updated 'HMM'' parameters and their corresponding likelihoods.
baumWelch' :: (Eq s, Eq o) => HMM' s o -> Vector o -> [(HMM' s o, Likelihood)]
baumWelch' model xs = zip ms $ tail ells
  where
    (ms, ells) = unzip $ iterate ((`baumWelch1'` xs) . fst) (model, undefined)

-- | Perform one step of the Baum-Welch EM algorithm and return the updated
--   'HMM'' parameters and the likelihood of the old parameters.
baumWelch1' :: (Eq s, Eq o) => HMM' s o -> Vector o -> (HMM' s o, Likelihood)
baumWelch1' model xs = (model', likelihood)
  where
    model' = HMM' { states'           = ss
                  , outputs'          = os
                  , initialStateDist' = pi0
                  , transitionDist'   = w
                  , emissionDistT'    = phi'
                  }
    likelihood = V.last ells

    -- First, we calculate the alpha and beta values using the
    -- forward-backward algorithm.
    alphas = forward' model xs
    betas  = backward' model xs

    -- Then, we obtain the likelihoods for each time. This should be
    -- constant over time.
    ells = V.zipWith (<.>) alphas betas

    -- Based on the alpha, beta, and likelihood values, we calculate the
    -- gamma and xi values.
    gammas = V.zipWith3 (\a b l -> a >.> b >/ l) alphas betas ells
    xis    = V.zipWith4 (\a b l x -> let w1 = V.zipWith (.>) a w0
                                         w2 = V.map (phi0 ! x >.> b >.>) w1
                                     in w2 #/ l)
               alphas (V.tail betas) (V.tail ells) (V.tail xs')
      where
        xs'  = V.map (`V.unsafeElemIndex` os) xs
        w0   = transitionDist' model
        phi0 = emissionDistT' model

    -- Using the gamma and xi values, we finally obtain the optimal initial
    -- state probability vector, transition probability matrix, and
    -- emission probability matrix.
    ss   = states' model
    os   = outputs' model
    pi0  = let gs = gammas ! 0
           in gs >/ V.sum gs
    w    = let ws = V.foldl1' (#+#) xis
               zs = V.map V.sum ws
           in V.zipWith (>/) ws zs
    phi' = let gs' o = V.map snd $ V.filter ((== o) . fst) $ V.zip xs gammas
               phis  = V.foldl1' (>+>) . gs'
               zs    = V.foldl1' (>+>) gammas
           in V.map (\o -> phis o >/> zs) os

-- | Baum-Welch forward algorithm that generates α values
forward' :: Eq o => HMM' s o -> Vector o -> Vector (Vector Probability)
forward' model xs = runST $ do
  v <- M.new n
  v `M.write` 0 $ (phi' ! x 0) >.> pi0
  forM_ [1..(n-1)] $ \i -> do
    a <- v `M.read` (i-1)
    v `M.write` i $ (phi' ! x i) >.> (a <.# w)
  V.freeze v
  where
    n   = V.length xs
    x i = let os  = outputs' model
              xs' = V.map (`V.unsafeElemIndex` os) xs
          in xs' ! i
    pi0  = initialStateDist' model
    w    = transitionDist' model
    phi' = emissionDistT' model

-- | Baum-Welch backward algorithm that generates β values
backward' :: Eq o => HMM' s o -> Vector o -> Vector (Vector Probability)
backward' model xs = runST $ do
  v <- M.new n
  v `M.write` (n-1) $ V.replicate k $ logFloat (1 :: Double)
  forM_ (reverse [0..(n-2)]) $ \i -> do
    b <- v `M.read` (i+1)
    v `M.write` i $ w #.> ((phi' ! x (i+1)) >.> b)
  V.freeze v
  where
    n   = V.length xs
    k   = V.length $ states' model
    x i = let os  = outputs' model
              xs' = V.map (`V.unsafeElemIndex` os) xs
          in xs' ! i
    w    = transitionDist' model
    phi' = emissionDistT' model
