-- | Operators commonly used in the basic linear algebra. Note that all the
--   functions defined here do not check the dimension/length of
--   vectors/matrices.
module Data.Vector.Generic.Util.LinearAlgebra (
  -- * Pairwise operators
    (>+>)
  -- , (>->)
  , (>.>)
  , (>/>)
  , (#+#)
  -- , (#-#)
  -- , (#.#)
  -- , (#/#)

  -- * Scalar-vector/vector-scalar operators
  -- , (+>)
  -- , (->)
  , (.>)
  -- , (/>)
  -- , (>+)
  -- , (>-)
  -- , (>.)
  , (>/)

  -- * Scalar-matrix/matrix-scalar operators
  -- , (+#)
  -- , (-#)
  -- , (.#)
  -- , (/#)
  -- , (#+)
  -- , (#-)
  -- , (#.)
  , (#/)

  -- * Dot and matrix-vector/vector-matrix products
  , (<.>)
  , (#.>)
  , (<.#)

  -- * Unary operators
  , transpose
  ) where

import Prelude hiding (any, map, null, sum, zipWith)
import Data.Vector.Generic (
    Vector, any, cons, convert, empty, map, null, sum, unsafeHead, unsafeTail
  , zipWith
  )

-- $setup
-- >>> :module + Data.Vector

-- | Pairwise addition between two vectors
--
-- >>> fromList [1, 2] >+> fromList [3, 4 :: Int]
-- fromList [4,6]
(>+>) :: (Num a, Vector v a) => v a -> v a -> v a
{-# INLINE (>+>) #-}
u >+> v = zipWith (+) u v

-- | Pairwise product between two vectors
--
-- >>> fromList [1, 2] >.> fromList [3, 4 :: Double]
-- fromList [3.0,8.0]
(>.>) :: (Num a, Vector v a) => v a -> v a -> v a
{-# INLINE (>.>) #-}
u >.> v = zipWith (*) u v

-- | Pairwise division between two vectors
--
-- >>> fromList [1, 2] >/> fromList [3, 4 :: Double]
-- fromList [0.3333333333333333,0.5]
(>/>) :: (Fractional a, Vector v a) => v a -> v a -> v a
{-# INLINE (>/>) #-}
u >/> v = zipWith (/) u v

-- | Pairwise addition between two matrices
--
-- >>> fromList [fromList [1, 2], fromList [3, 4]] #+# fromList [fromList [5, 6], fromList [7, 8 :: Int]]
-- fromList [fromList [6,8],fromList [10,12]]
(#+#) :: (Num a, Vector v a, Vector w (v a)) => w (v a) -> w (v a) -> w (v a)
{-# INLINE (#+#) #-}
m #+# n = zipWith (>+>) m n

-- | Scalar-vector product
--
-- >>> 2 .> fromList [1, 2 :: Integer]
-- fromList [2,4]
(.>) :: (Num a, Vector v a) => a -> v a -> v a
{-# INLINE (.>) #-}
s .> v = map (s *) v

-- | Vector-scalar division
--
-- >>> fromList [1, 2 :: Double] >/ 2
-- fromList [0.5,1.0]
(>/) :: (Fractional a, Vector v a) => v a -> a -> v a
{-# INLINE (>/) #-}
v >/ s = map (/ s) v

-- | Matrix-scalar division
--
-- >>> fromList [fromList [1, 2], fromList [3, 4 :: Double]] #/ 2
-- fromList [fromList [0.5,1.0],fromList [1.5,2.0]]
(#/) :: (Fractional a, Vector v a, Vector w (v a)) => w (v a) -> a -> w (v a)
{-# INLINE (#/) #-}
m #/ s = map (>/ s) m

-- | Dot product
--
-- >>> fromList [1, 2] <.> fromList [3, 4 :: Int]
-- 11
(<.>) :: (Num a, Vector v a) => v a -> v a -> a
{-# INLINE (<.>) #-}
u <.> v = sum $ u >.> v

-- | Matrix-vector product
--
-- >>> fromList [fromList [1, 2], fromList [3, 4]] #.> fromList [1, 2 :: Double]
-- fromList [5.0,11.0]
(#.>) :: (Num a, Vector v a, Vector w (v a), Vector w a) => w (v a) -> v a -> v a
{-# INLINE (#.>) #-}
m #.> v = convert $ map (<.> v) m

-- | Vector-matrix product
--
-- >>> fromList [1, 2 :: Double] <.# fromList [fromList [1, 2], fromList [3, 4]]
-- fromList [7.0,10.0]
(<.#) :: (Num a, Vector v a, Vector w (v a), Vector w a) => v a -> w (v a) -> v a
{-# INLINE (<.#) #-}
v <.# m | any null m = empty
        | otherwise  = hd `cons` tl
  where
    hd = v <.> convert (map unsafeHead m)
    tl = v <.# map unsafeTail m

-- | Matrix transpose
--
-- >>> transpose $ fromList [fromList "ab", fromList "cd"]
-- fromList [fromList "ac",fromList "bd"]
transpose :: (Vector v a, Vector w (v a), Vector w a) => w (v a) -> w (v a)
{-# INLINE transpose #-}
transpose m
  | any null m = empty
  | otherwise  = hd `cons` tl
  where
    hd = convert $ map unsafeHead m
    tl = transpose $ map unsafeTail m
