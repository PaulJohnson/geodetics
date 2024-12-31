{- |
The following is based on equations in Section 1.4.7.1 in 
OGP Surveying and Positioning Guidance Note number 7, part 2 – August 2006
<http://ftp.stu.edu.tw/BSD/NetBSD/pkgsrc/distfiles/epsg-6.11/G7-2.pdf>
-}

{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
module Geodetics.Stereographic (
   GridStereo (gridTangent, gridOrigin, gridScale),
   mkGridStereo
) where

import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid

import qualified Data.Stream as Stream

-- | A stereographic projection with its origin at an arbitrary point on Earth, other than the poles.
data GridStereo e = GridStereo {
      gridTangent :: Geodetic e, -- ^ Point where the plane of projection touches the ellipsoid. Often known as the Natural Origin.
      gridOrigin :: GridOffset,  -- ^ Grid position of the tangent point. Often known as the False Origin.
      gridScale :: Double, -- ^ Scaling factor that balances the distortion between the center and the edges.
                                         -- Should be slightly less than unity.
      
      -- Memoised parameters derived from the tangent point.
      gridR :: Double,
      gridN, gridC, gridSin, gridCos :: Double,
      gridLatC :: Double,
      gridG, gridH :: Double
   } deriving (Show)
   
-- | Create a stereographic projection. The tangency point must not be one of the poles.  
mkGridStereo :: (Ellipsoid e) => Geodetic e -> GridOffset -> Double -> GridStereo e
mkGridStereo tangent origin scale = GridStereo {
      gridTangent = tangent,
      gridOrigin = origin,
      gridScale = scale,
      gridR = r,
      gridN = n,
      gridC = c,
      gridSin = sinLatC1,
      gridCos = sqrt $ 1 - sinLatC1 * sinLatC1,
      gridLatC = asin sinLatC1,
      gridG = g,
      gridH = h
   }
   where 
      -- The reference seems to use χO to refer to two slightly different values. 
      -- Here these will be called LatC0 and LatC1.
      ellipse = ellipsoid tangent
      op :: Num a => a -> a    -- Values of longitude, tangent longitude, E and N
      op = if latitude tangent < 0 then negate else id  -- must be negated in the southern hemisphere.
      lat0 = op $ latitude tangent
      sinLat0 = sin lat0
      e2 = eccentricity2 ellipse
      e = sqrt e2
      r = sqrt $ meridianRadius ellipse lat0 * primeVerticalRadius ellipse lat0
      n = sqrt $ 1 + ((e2 * cos lat0 ^ _4)/(1 - e2))
      s1 = (1 + sinLat0) / (1 - sinLat0)
      s2 = (1 - e * sinLat0) / (1 + e * sinLat0)
      w1 = (s1 * s2 ** e) ** n
      sinLatC0 = (w1 - 1)/(w1 + 1)
      c = ((n + sin lat0) * (1 - sinLatC0)) / ((n - sin lat0) * (1 + sinLatC0))
      w2 = c * w1
      sinLatC1 = (w2 - 1)/(w2 + 1)
      g = 2 * r * scale * tan (pi/4 - latC1/2)
      h = 4 * r * scale * tan latC1 + g
      latC1 = asin sinLatC1
      

instance (Ellipsoid e) => GridClass (GridStereo e) e where
   toGrid grid geo = applyOffset (gridOrigin grid) $ GridPoint east north (geoAlt geo) grid
      where
         op :: Num a => a -> a    -- Values of longitude, tangent longitude, E and N
         op = if latitude (gridTangent grid) < 0 then negate else id  -- must be negated in the southern hemisphere.
         sinLatC = (w - 1)/(w + 1)
         cosLatC = sqrt $ 1 - sinLatC * sinLatC
         longC = gridN grid * (op (longitude geo) - long0) + long0
         w = gridC grid * (sA * sB ** e) ** gridN grid
         sA = (1+sinLat) / (1 - sinLat)
         sB = (1 - e*sinLat) / (1 + e*sinLat)
         sinLat = sin $ op $ latitude geo
         e = sqrt $ eccentricity2 $ ellipsoid geo
         long0 = op $ longitude $ gridTangent grid
         b = 1 + sinLatC * gridSin grid + cosLatC * gridCos grid * cos (longC - long0)
         east = 2 * gridR grid * gridScale grid * cosLatC * sin (longC - long0) / b
         north = 2 * gridR grid * gridScale grid * (sinLatC * gridCos grid - cosLatC * gridSin grid * cos (longC - long0)) / b
   
   fromGrid gp = 
      {- trace (    -- Remove comment brackets for debugging.
         "fromGrid values:\n   i = " ++ show i ++ "\n   j = " ++ show j ++
         "\n   longC = " ++ show longC ++ "\n   long = " ++ show long ++
         "\n   latC = " ++ show latC ++
         "\n   lat1 = " ++ show lat1 ++ "\n   latN = " ++ show latN ) $ -}
         Geodetic (op latN) (op long) height $ gridEllipsoid grid
      where
         op :: Num a => a -> a                   -- Values of longitude, tangent longitude, E and N
         op = if latitude (gridTangent grid) < 0 then negate else id  -- must be negated in the southern hemisphere.
         GridPoint east north height _ = applyOffset (offsetNegate $ gridOrigin grid) gp
         east' = east
         north' = north
         grid = gridBasis gp
         long0 = op $ longitude $ gridTangent grid
         i = atan2 east' (gridH grid + north')
         j = atan2 east' (gridG grid - north') - i
         latC = gridLatC grid + 2 * atan2 (north' - east' * tan (j/2)) (2 * gridR grid * gridScale grid)
         longC = j + 2 * i + long0
         sinLatC = sin latC
         long = (longC - long0) / gridN grid + long0
         isoLat = log ((1 + sinLatC) / (gridC grid * (1 - sinLatC))) / (2 * gridN grid)
         lat1 = 2 * atan (exp isoLat) - pi/2
         next lat = lat - (isoN - isoLat) * cos lat * (1 - e2 * sin lat ^ _2) / (1 - e2)
            where isoN = isometricLatitude (gridEllipsoid grid) lat
                  e2 = eccentricity2 $ gridEllipsoid grid
         lats = Stream.iterate next lat1
         latN = snd $ Stream.head $ Stream.dropWhile (\(v1, v2) -> abs (v1-v2) > 0.01 * arcsecond) $ Stream.zip lats $ Stream.drop 1 lats
   gridEllipsoid = ellipsoid . gridTangent
