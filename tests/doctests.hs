module Main (main) where

import Test.DocTest

main :: IO ()
main = doctest [ "-isrc"
               , "src/Data/Vector/Util/LinearAlgebra.hs"
               , "src/Learning/HMM.hs"
               ]
