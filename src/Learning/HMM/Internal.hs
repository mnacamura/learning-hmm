module Learning.HMM.Internal (
    HMM (..)
  , LogLikelihood
  , init
  , withEmission
  , viterbi
  , baumWelch
  , baumWelch'
  -- , baumWelch1
  -- , forward
  -- , backward
  -- , posterior
  ) where

import Prelude hiding (init)
import Control.Applicative ((<$>))
import Control.DeepSeq (NFData, force, rnf)
import Control.Monad (forM_, replicateM)
import Control.Monad.ST (runST)
import qualified Data.Map.Strict as M (findWithDefault)
import Data.Random.RVar (RVar)
import Data.Random.Distribution.Simplex (stdSimplex)
import qualified Data.Vector as V (
    Vector, filter, foldl1', map, unsafeFreeze, unsafeIndex, unsafeTail
  , zip, zipWith3
  )
import qualified Data.Vector.Generic as G (convert)
import qualified Data.Vector.Generic.Extra as G (frequencies)
import qualified Data.Vector.Mutable as MV (
    unsafeNew, unsafeRead, unsafeWrite
  )
import qualified Data.Vector.Unboxed as U (
    Vector, fromList, length, map, sum, unsafeFreeze, unsafeIndex
  , unsafeTail, zip
  )
import qualified Data.Vector.Unboxed.Mutable as MU (
    unsafeNew, unsafeRead, unsafeWrite
  )
import qualified Numeric.LinearAlgebra.Data as H (
    (!), Matrix, Vector, diag, fromColumns, fromList, fromLists
  , fromRows, konst, maxElement, maxIndex, toColumns, tr
  )
import qualified Numeric.LinearAlgebra.HMatrix as H (
    (<>), (#>), sumElements
  )

type LogLikelihood = Double

-- | More efficient data structure of the 'HMM' model. The 'states' and
--   'outputs' in 'HMM' are represented by their indices. The
--   'initialStateDist', 'transitionDist', and 'emissionDist' are
--   represented by matrices. The 'emissionDistT' is a transposed matrix
--   in order to simplify the calculation.
data HMM = HMM { nStates          :: Int -- ^ Number of states
               , nOutputs         :: Int -- ^ Number of outputs
               , initialStateDist :: H.Vector Double
               , transitionDist   :: H.Matrix Double
               , emissionDistT    :: H.Matrix Double
               }

instance NFData HMM where
  rnf hmm = rnf k `seq` rnf l `seq` rnf pi0 `seq` rnf w `seq` rnf phi'
    where
      k    = nStates hmm
      l    = nOutputs hmm
      pi0  = initialStateDist hmm
      w    = transitionDist hmm
      phi' = emissionDistT hmm

init :: Int -> Int -> RVar HMM
init k l = do
  pi0 <- H.fromList <$> stdSimplex (k-1)
  w   <- H.fromLists <$> replicateM k (stdSimplex (k-1))
  phi <- H.fromLists <$> replicateM k (stdSimplex (l-1))
  return HMM { nStates          = k
             , nOutputs         = l
             , initialStateDist = pi0
             , transitionDist   = w
             , emissionDistT    = H.tr phi
             }

withEmission :: HMM -> U.Vector Int -> HMM
withEmission model xs = model'
  where
    n = U.length xs
    k = nStates model
    l = nOutputs model
    ss = [0..(k-1)]
    os = [0..(l-1)]

    step m = fst $ baumWelch1 (m { emissionDistT = H.tr phi }) n xs
      where
        phi :: H.Matrix Double
        phi = let zs  = fst $ viterbi m xs
                  fs  = G.frequencies $ U.zip zs xs
                  hs  = H.fromLists $ map (\s -> map (\o ->
                          M.findWithDefault 0 (s, o) fs) os) ss
                  -- hs' is needed to not yield NaN vectors
                  hs' = hs + H.konst 1e-9 (k, l)
                  ns  = hs' H.#> H.konst 1 k
              in hs' / H.fromColumns (replicate l ns)

    ms  = iterate step model
    ms' = tail ms
    ds  = zipWith euclideanDistance ms ms'

    model' = fst $ head $ dropWhile ((> 1e-9) . snd) $ zip ms' ds

-- | Return the Euclidean distance between two models.
euclideanDistance :: HMM -> HMM -> Double
euclideanDistance model model' =
  sqrt $ H.sumElements ((w - w') ** 2) + H.sumElements ((phi - phi') ** 2)
  where
    w    = transitionDist model
    w'   = transitionDist model'
    phi  = emissionDistT model
    phi' = emissionDistT model'

viterbi :: HMM -> U.Vector Int -> (U.Vector Int, LogLikelihood)
viterbi model xs = (path, logL)
  where
    n = U.length xs

    -- First, we calculate the value function and the state maximizing it
    -- for each time.
    deltas :: V.Vector (H.Vector Double)
    psis   :: V.Vector (U.Vector Int)
    (deltas, psis) = runST $ do
      ds <- MV.unsafeNew n
      ps <- MV.unsafeNew n
      let x0 = U.unsafeIndex xs 0
      MV.unsafeWrite ds 0 $ log (phi' H.! x0) + log pi0
      forM_ [1..(n-1)] $ \t -> do
        d <- MV.unsafeRead ds (t-1)
        let x   = U.unsafeIndex xs t
            dws = map (\wj -> d + log wj) w'
        MV.unsafeWrite ds t $ log (phi' H.! x) + H.fromList (map H.maxElement dws)
        MV.unsafeWrite ps t $ U.fromList (map H.maxIndex dws)
      ds' <- V.unsafeFreeze ds
      ps' <- V.unsafeFreeze ps
      return (ds', ps')
      where
        pi0  = initialStateDist model
        w'   = H.toColumns $ transitionDist model
        phi' = emissionDistT model

    deltaE = V.unsafeIndex deltas (n-1)

    -- The most likely path and corresponding log likelihood are as follows.
    path = runST $ do
      ix <- MU.unsafeNew n
      MU.unsafeWrite ix (n-1) $ H.maxIndex deltaE
      forM_ [n-l | l <- [1..(n-1)]] $ \t -> do
        i <- MU.unsafeRead ix t
        let psi = V.unsafeIndex psis t
        MU.unsafeWrite ix (t-1) $ U.unsafeIndex psi i
      U.unsafeFreeze ix
    logL = H.maxElement deltaE

baumWelch :: HMM -> U.Vector Int -> [(HMM, LogLikelihood)]
baumWelch model xs = zip models (tail logLs)
  where
    n = U.length xs
    step (m, _)     = baumWelch1 m n xs
    (models, logLs) = unzip $ iterate step (model, undefined)

baumWelch' :: HMM -> U.Vector Int -> (HMM, LogLikelihood)
baumWelch' model xs = go (undefined, -1/0) (baumWelch1 model n xs)
  where
    n = U.length xs
    go (m, l) (m', l')
      | l' - l > 1.0e-9 = go (m', l') (baumWelch1 m' n xs)
      | otherwise       = (m, l')

-- | Perform one step of the Baum-Welch algorithm and return the updated
--   model and the likelihood of the old model.
baumWelch1 :: HMM -> Int -> U.Vector Int -> (HMM, LogLikelihood)
baumWelch1 model n xs = force (model', logL)
  where
    k = nStates model
    l = nOutputs model

    -- First, we calculate the alpha, beta, and scaling values using the
    -- forward-backward algorithm.
    (alphas, cs) = forward model n xs
    betas        = backward model n xs cs

    -- Based on the alpha, beta, and scaling values, we calculate the
    -- posterior distribution, i.e., gamma and xi values.
    (gammas, xis) = posterior model n xs alphas betas cs

    -- Using the gamma and xi values, we obtain the optimal initial state
    -- probability vector, transition probability matrix, and emission
    -- probability matrix.
    pi0  = V.unsafeIndex gammas 0
    w    = let ds = V.foldl1' (+) xis   -- denominators
               ns = ds H.#> H.konst 1 k -- numerators
           in H.diag (H.konst 1 k / ns) H.<> ds
    phi' = let gs' o = V.map snd $ V.filter ((== o) . fst) $ V.zip (G.convert xs) gammas
               ds    = V.foldl1' (+) . gs'  -- denominators
               ns    = V.foldl1' (+) gammas -- numerators
           in H.fromRows $ map (\o -> ds o / ns) [0..(l-1)]

    -- We finally obtain the new model and the likelihood for the old model.
    model' = model { initialStateDist = pi0
                   , transitionDist   = w
                   , emissionDistT    = phi'
                   }
    logL = - (U.sum $ U.map log cs)

-- | Return alphas and scaling variables.
forward :: HMM -> Int -> U.Vector Int -> (V.Vector (H.Vector Double), U.Vector Double)
{-# INLINE forward #-}
forward model n xs = runST $ do
  as <- MV.unsafeNew n
  cs <- MU.unsafeNew n
  let x0 = U.unsafeIndex xs 0
      a0 = (phi' H.! x0) * pi0
      c0 = 1 / H.sumElements a0
  MV.unsafeWrite as 0 (H.konst c0 k * a0)
  MU.unsafeWrite cs 0 c0
  forM_ [1..(n-1)] $ \t -> do
    a <- MV.unsafeRead as (t-1)
    let x  = U.unsafeIndex xs t
        a' = (phi' H.! x) * (w' H.#> a)
        c' = 1 / H.sumElements a'
    MV.unsafeWrite as t (H.konst c' k * a')
    MU.unsafeWrite cs t c'
  as' <- V.unsafeFreeze as
  cs' <- U.unsafeFreeze cs
  return (as', cs')
  where
    k    = nStates model
    pi0  = initialStateDist model
    w'   = H.tr $ transitionDist model
    phi' = emissionDistT model

-- | Return betas using scaling variables.
backward :: HMM -> Int -> U.Vector Int -> U.Vector Double -> V.Vector (H.Vector Double)
{-# INLINE backward #-}
backward model n xs cs = runST $ do
  bs <- MV.unsafeNew n
  let bE = H.konst 1 k
      cE = U.unsafeIndex cs (n-1)
  MV.unsafeWrite bs (n-1) (H.konst cE k * bE)
  forM_ [n-l | l <- [1..(n-1)]] $ \t -> do
    b <- MV.unsafeRead bs t
    let x  = U.unsafeIndex xs t
        b' = w H.#> ((phi' H.! x) * b)
        c' = U.unsafeIndex cs (t-1)
    MV.unsafeWrite bs (t-1) (H.konst c' k * b')
  V.unsafeFreeze bs
  where
    k    = nStates model
    w    = transitionDist model
    phi' = emissionDistT model

-- | Return the posterior distribution.
posterior :: HMM -> Int -> U.Vector Int -> V.Vector (H.Vector Double) -> V.Vector (H.Vector Double) -> U.Vector Double -> (V.Vector (H.Vector Double), V.Vector (H.Matrix Double))
{-# INLINE posterior #-}
posterior model _ xs alphas betas cs = (gammas, xis)
  where
    gammas = V.zipWith3 (\a b c -> a * b / H.konst c k)
               alphas betas (G.convert cs)
    xis    = V.zipWith3 (\a b x -> H.diag a H.<> w H.<> H.diag (b * (phi' H.! x)))
               alphas (V.unsafeTail betas) (G.convert $ U.unsafeTail xs)
    k    = nStates model
    w    = transitionDist model
    phi' = emissionDistT model
