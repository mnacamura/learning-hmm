-- | Miscellaneous utility functions for "Data.Vector"
module Data.Vector.Util (
    unsafeElemIndex
  ) where

import Data.Maybe (fromJust)
import Data.Vector (Vector, elemIndex)

-- | Return the index of the first occurrence of the given element or throw
--   an error if no such occurrence.
{-# INLINE unsafeElemIndex #-}
unsafeElemIndex :: Eq a => a -> Vector a -> Int
unsafeElemIndex e = fromJust . elemIndex e
