module Main (main) where

import Test.DocTest

main :: IO ()
main = doctest [ "-isrc"
               , "src/Data/Vector/Generic/Util/LinearAlgebra.hs"
               , "src/Learning/HMM.hs"
               ]
