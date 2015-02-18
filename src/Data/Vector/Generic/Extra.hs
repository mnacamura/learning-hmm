-- | Extra functions for "Data.Vector"
module Data.Vector.Generic.Extra (
    frequencies
  ) where

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as M (empty, insertWith)
import Data.Vector.Generic (Vector, foldl')

-- $setup
-- >>> :module + Data.Vector

-- | @frequencies xs@ returns a 'Map' from distinct items in @xs@ to
-- the number of times they appear.
--
-- >>> frequencies $ fromList "bra bra bar"
-- fromList [(' ',2),('a',3),('b',3),('r',3)]
frequencies :: (Ord a, Vector v a, Num n) => v a -> Map a n
frequencies = foldl' (\m k -> M.insertWith (+) k 1 m) M.empty
