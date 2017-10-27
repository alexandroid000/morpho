module Morpho where

import qualified    Data.Graph.Inductive as Graph
import              Data.Graph.Inductive (Gr, match)

import Data.HashMap
import Graphics.Triangulation.Delaunay (triangulate)

import System.Random


mkPoints :: Int -> IO [(Double, Double)]
mkPoints n = getStdGen >>= return . (rndPoints n)


rndPoints :: Int -> StdGen -> [(Double, Double)]
rndPoints n sg = take n $ zip x y
    where   (sg1, sg2) = split sg
            x = randomRs (-1.0, 1.0) sg1
            y = randomRs (-1.0, 1.0) sg2

someFunc :: IO ()
someFunc = putStrLn "someFunc"
