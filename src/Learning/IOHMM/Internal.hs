module Learning.IOHMM.Internal (
    IOHMM (..)
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
    Vector, filter, foldl1', generate, map, replicateM, unsafeFreeze
  , unsafeIndex , unsafeTail , zip, zipWith3
  )
import qualified Data.Vector.Generic as G (convert)
import qualified Data.Vector.Generic.Util as G (frequencies)
import qualified Data.Vector.Mutable as MV (
    unsafeNew, unsafeRead, unsafeWrite
  )
import qualified Data.Vector.Unboxed as U (
    Vector, fromList, length, map, sum, unsafeFreeze, unsafeIndex
  , unsafeTail, unzip, zip
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

-- | More efficient data structure of the 'IOHMM' model. The 'inputs',
--   'states', and 'outputs' in 'IOHMM' are represented by their indices.
--   The 'initialStateDist', 'transitionDist', and 'emissionDist' are
--   represented by matrices. The 'emissionDistT' is a transposed matrix
--   in order to simplify the calculation.
data IOHMM = IOHMM { nInputs          :: Int -- ^ Number of inputs
                   , nStates          :: Int -- ^ Number of states
                   , nOutputs         :: Int -- ^ Number of outputs
                   , initialStateDist :: H.Vector Double
                   , transitionDist   :: V.Vector (H.Matrix Double)
                   , emissionDistT    :: H.Matrix Double
                   }

instance NFData IOHMM where
  rnf hmm = rnf m `seq` rnf k `seq` rnf l `seq` rnf pi0 `seq` rnf w `seq` rnf phi'
    where
      m    = nInputs hmm
      k    = nStates hmm
      l    = nOutputs hmm
      pi0  = initialStateDist hmm
      w    = transitionDist hmm
      phi' = emissionDistT hmm

init :: Int -> Int -> Int -> RVar IOHMM
init m k l = do
  pi0 <- H.fromList <$> stdSimplex (k-1)
  w   <- V.replicateM m (H.fromLists <$> replicateM k (stdSimplex (k-1)))
  phi <- H.fromLists <$> replicateM k (stdSimplex (l-1))
  return IOHMM { nInputs          = m
               , nStates          = k
               , nOutputs         = l
               , initialStateDist = pi0
               , transitionDist   = w
               , emissionDistT    = H.tr phi
               }

withEmission :: IOHMM -> U.Vector (Int, Int) -> IOHMM
withEmission model xys = model'
  where
    n = U.length xys
    k = nStates model
    l = nOutputs model
    ss = [0..(k-1)]
    os = [0..(l-1)]
    ys = U.map snd xys

    step m = fst $ baumWelch1 (m { emissionDistT = H.tr phi }) n xys
      where
        phi :: H.Matrix Double
        phi = let zs  = fst $ viterbi m xys
                  fs  = G.frequencies $ U.zip zs ys
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
euclideanDistance :: IOHMM -> IOHMM -> Double
euclideanDistance model model' =
  sqrt $ sum $ H.sumElements ((phi - phi') ** 2) :
               map (\i -> H.sumElements ((w i - w' i) ** 2)) is
  where
    is   = [0..(nInputs model - 1)]
    w    = V.unsafeIndex (transitionDist model)
    w'   = V.unsafeIndex (transitionDist model')
    phi  = emissionDistT model
    phi' = emissionDistT model'

viterbi :: IOHMM -> U.Vector (Int, Int) -> (U.Vector Int, LogLikelihood)
viterbi model xys = (path, logL)
  where
    n = U.length xys

    -- First, we calculate the value function and the state maximizing it
    -- for each time.
    deltas :: V.Vector (H.Vector Double)
    psis   :: V.Vector (U.Vector Int)
    (deltas, psis) = runST $ do
      ds <- MV.unsafeNew n
      ps <- MV.unsafeNew n
      let (_, y0) = U.unsafeIndex xys 0
      MV.unsafeWrite ds 0 $ log (phi' H.! y0) + log pi0
      forM_ [1..(n-1)] $ \t -> do
        d <- MV.unsafeRead ds (t-1)
        let (x, y) = U.unsafeIndex xys t
            dws    = map (\wj -> d + log wj) (w' x)
        MV.unsafeWrite ds t $ log (phi' H.! y) + H.fromList (map H.maxElement dws)
        MV.unsafeWrite ps t $ U.fromList (map H.maxIndex dws)
      ds' <- V.unsafeFreeze ds
      ps' <- V.unsafeFreeze ps
      return (ds', ps')
      where
        pi0  = initialStateDist model
        w'   = H.toColumns . V.unsafeIndex (transitionDist model)
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

baumWelch :: IOHMM -> U.Vector (Int, Int) -> [(IOHMM, LogLikelihood)]
baumWelch model xys = zip models (tail logLs)
  where
    n = U.length xys
    step (m, _)     = baumWelch1 m n xys
    (models, logLs) = unzip $ iterate step (model, undefined)

baumWelch' :: IOHMM -> U.Vector (Int, Int) -> (IOHMM, LogLikelihood)
baumWelch' model xys = go (undefined, -1/0) (baumWelch1 model n xys)
  where
    n = U.length xys
    go (m, l) (m', l')
      | l' - l > 1.0e-9 = go (m', l') (baumWelch1 m' n xys)
      | otherwise       = (m, l')

-- | Perform one step of the Baum-Welch algorithm and return the updated
--   model and the likelihood of the old model.
baumWelch1 :: IOHMM -> Int -> U.Vector (Int, Int) -> (IOHMM, LogLikelihood)
baumWelch1 model n xys = force (model', logL)
  where
    m = nInputs model
    k = nStates model
    l = nOutputs model
    (xs, ys) = U.unzip xys

    -- First, we calculate the alpha, beta, and scaling values using the
    -- forward-backward algorithm.
    (alphas, cs) = forward model n xys
    betas        = backward model n xys cs

    -- Based on the alpha, beta, and scaling values, we calculate the
    -- posterior distribution, i.e., gamma and xi values.
    (gammas, xis) = posterior model n xys alphas betas cs

    -- Using the gamma and xi values, we obtain the optimal initial state
    -- probability vector, transition probability matrix, and emission
    -- probability matrix.
    pi0  = V.unsafeIndex gammas 0
    w    = let xis' i = V.map snd $ V.filter ((== i) . fst) $ V.zip (G.convert $ U.unsafeTail xs) xis
               ds     = V.foldl1' (+) . xis'  -- denominators
               ns i   = ds i H.#> H.konst 1 k -- numerators
           in V.map (\i -> H.diag (H.konst 1 k / ns i) H.<> ds i) (V.generate m id)
    phi' = let gs' o = V.map snd $ V.filter ((== o) . fst) $ V.zip (G.convert ys) gammas
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
forward :: IOHMM -> Int -> U.Vector (Int, Int) -> (V.Vector (H.Vector Double), U.Vector Double)
{-# INLINE forward #-}
forward model n xys = runST $ do
  as <- MV.unsafeNew n
  cs <- MU.unsafeNew n
  let (_, y0) = U.unsafeIndex xys 0
      a0      = (phi' H.! y0) * pi0
      c0      = 1 / H.sumElements a0
  MV.unsafeWrite as 0 (H.konst c0 k * a0)
  MU.unsafeWrite cs 0 c0
  forM_ [1..(n-1)] $ \t -> do
    a <- MV.unsafeRead as (t-1)
    let (x, y) = U.unsafeIndex xys t
        a'     = (phi' H.! y) * (w' x H.#> a)
        c'     = 1 / H.sumElements a'
    MV.unsafeWrite as t (H.konst c' k * a')
    MU.unsafeWrite cs t c'
  as' <- V.unsafeFreeze as
  cs' <- U.unsafeFreeze cs
  return (as', cs')
  where
    k    = nStates model
    pi0  = initialStateDist model
    w'   = H.tr . V.unsafeIndex (transitionDist model)
    phi' = emissionDistT model

-- | Return betas using scaling variables.
backward :: IOHMM -> Int -> U.Vector (Int, Int) -> U.Vector Double -> V.Vector (H.Vector Double)
{-# INLINE backward #-}
backward model n xys cs = runST $ do
  bs <- MV.unsafeNew n
  let bE = H.konst 1 k
      cE = U.unsafeIndex cs (n-1)
  MV.unsafeWrite bs (n-1) (H.konst cE k * bE)
  forM_ [n-l | l <- [1..(n-1)]] $ \t -> do
    b <- MV.unsafeRead bs t
    let (x, y) = U.unsafeIndex xys t
        b'     = w x H.#> ((phi' H.! y) * b)
        c'     = U.unsafeIndex cs (t-1)
    MV.unsafeWrite bs (t-1) (H.konst c' k * b')
  V.unsafeFreeze bs
  where
    k    = nStates model
    w    = V.unsafeIndex (transitionDist model)
    phi' = emissionDistT model

-- | Return the posterior distribution.
posterior :: IOHMM -> Int -> U.Vector (Int, Int) -> V.Vector (H.Vector Double) -> V.Vector (H.Vector Double) -> U.Vector Double -> (V.Vector (H.Vector Double), V.Vector (H.Matrix Double))
{-# INLINE posterior #-}
posterior model _ xys alphas betas cs = (gammas, xis)
  where
    gammas = V.zipWith3 (\a b c -> a * b / H.konst c k)
               alphas betas (G.convert cs)
    xis    = V.zipWith3 (\a b (x, y) -> H.diag a H.<> w x H.<> H.diag (b * (phi' H.! y)))
               alphas (V.unsafeTail betas) (G.convert $ U.unsafeTail xys)
    k    = nStates model
    w    = V.unsafeIndex (transitionDist model)
    phi' = emissionDistT model
