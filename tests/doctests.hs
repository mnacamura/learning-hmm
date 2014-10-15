import Test.DocTest

main :: IO ()
main = doctest [ "-isrc"
               , "src/Data/Vector/Util/LinearAlgebra.hs"
               ]
