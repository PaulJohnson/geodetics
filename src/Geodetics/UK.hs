{-# LANGUAGE MultiParamTypeClasses #-}

-- | Distinguished coordinate systems for the United Kingdom.
module Geodetics.UK (
   OSGB36 (..),
   UkNationalGrid (..),
   fromUkGridReference,
   ukTrueOrigin,
   ukFalseOrigin,
   fromUkGridLetters,
   fromUkGridDigits,
   toUkGridReference
) where

import Control.Applicative
import Data.Array
import Data.Char
import Data.Monoid
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P

import Geodetics.Coordinates
import Geodetics.Ellipsoids


-- | Ellipsoid definition for Great Britain. Airy 1830 offset from the centre of the Earth 
-- and rotated slightly.

-- The Helmert parameters are from the Ordnance Survey document 
-- \"A Guide to Coordinate Systems in Great Britain\", which notes that it
-- can be in error by as much as 5 meters and should not be used in applications
-- requiring greater accuracy.  A more precise conversion requires a large table 
-- of corrections for historical inaccuracies in the triangulation of the UK.
data OSGB36 = OSGB36 deriving (Eq, Show)

instance Ellipsoid OSGB36 where
   majorRadius _ = 6377563.396 *~ meter
   flatR _ = 299.3249646 *~ one
   helmert _ = Helmert {
      cX = 446.448 *~ meter, cY = (-125.157) *~ meter, cZ = 542.06 *~ meter,
      helmertScale = (-20.4894) *~ one,
      rX = 0.1502 *~ arcsecond, rY = 0.247 *~ arcsecond, rZ = 0.8421 *~ arcsecond }


-- | The UK National is a Transverse Mercator projection with a true origin at
-- 42 degrees North, 2 degrees West on OSGB36, and a false origin 400km West and 100 km North of
-- the true origin. The scale factor is 0.9996012717.
data UkNationalGrid = UkNationalGrid deriving (Eq, Show)

instance GridClass UkNationalGrid OSGB36 where
   toGrid _ = fromPrivateGrid . toGrid ukGrid
   fromGrid = fromGrid . toPrivateGrid
   gridEllipsoid _ = OSGB36


-- ------------------------------------------------------------------------------
-- Private definitions used in the instance of UkNationalGrid.

toPrivateGrid :: GridPoint UkNationalGrid -> GridPoint (GridTM OSGB36)
toPrivateGrid p = GridPoint (eastings p) (northings p) (altitude p) ukGrid

fromPrivateGrid :: GridPoint (GridTM OSGB36) -> GridPoint UkNationalGrid
fromPrivateGrid p = GridPoint (eastings p) (northings p) (altitude p) UkNationalGrid



ukTrueOrigin :: Geodetic OSGB36
ukTrueOrigin = Geodetic {
   latitude = 42 *~ degree,
   longitude = (-2) *~ degree,
   geoAlt = 0 *~ meter,
   ellipsoid = OSGB36
}

ukFalseOrigin :: GridOffset 
ukFalseOrigin = GridOffset (400 *~ kilo meter) ((-100) *~ kilo meter) (0 *~ meter)

ukGrid :: GridTM OSGB36
ukGrid = mkGridTM ukTrueOrigin ukFalseOrigin (0.9996012717 *~ one)


-- | Convert a grid reference to a position, if the reference is valid. 
-- This actually returns the position of the south-west corner of the nominated 
-- grid square and an offset to its centre. Altitude is set to zero.
fromUkGridReference :: String -> Maybe (GridPoint UkNationalGrid, GridOffset)
fromUkGridReference str = if length str < 2 then Nothing else do
      let c1:c2:ds = str
      (sw1, dc) <- fromUkGridDigits ds
      sw2 <- applyOffset sw1 <$> fromUkGridLetters c1 c2
      return (sw2, dc)

      


-- | Convert an even number of digits into a pair of grid offsets. The first result
-- is the south-west corner of the square, and the second is a further offset from
-- the corner to the centre.
fromUkGridDigits :: String -> Maybe (GridOffset, GridOffset)
fromUkGridDigits str = do
      let n = length str
          halfSquare = 50 *~ kilo meter * (0.1 *~ one) ** (fromIntegral (n `div` 2) *~ one)
      (s1,s2) <- if even n && all isDigit str then Just $ splitAt (n `div` 2) str else Nothing
      return $ (GridOffset (grid s1) (grid s2) (0 *~ meter),
                GridOffset halfSquare halfSquare (0 *~ meter))
   where
      grid :: String -> Length Double
      grid s = sum $ zipWith (*) 
         (map ((*~ one) . fromIntegral . digitToInt) s) 
         (iterate (/ (10 *~ one)) (10 *~ kilo meter)) 

-- | The south west corner of the nominated grid square, if it is a legal square.
-- This function works for all pairs of letters except 'I' (as that is not used).
-- In practice only those pairs covering the UK are actually used.
fromUkGridLetters :: Char -> Char -> Maybe (GridPoint UkNationalGrid)
fromUkGridLetters c1 c2 = applyOffset <$> (mappend <$> g1 <*> g2) <*> letterOrigin
   where
      letterOrigin = Just $ GridPoint ((-1000) *~ kilo meter) ((-500) *~ kilo meter) m0 UkNationalGrid
      gridIndex c = 
         if inRange ('A', 'H') c then Just $ ord c P.- ord 'A'  -- 'I' is not used.
         else if inRange ('J', 'Z') c then Just $ ord c P.- ord 'B'
         else Nothing
      gridSquare c = do -- Maybe monad
         g <- gridIndex c
         let (y,x) = g `divMod` 5 
         return (fromIntegral x *~ one, _4 - fromIntegral y *~ one)
      g1 = do
         (x,y) <- gridSquare c1
         return $ GridOffset (x * (500 *~ kilo meter)) (y * (500 *~ kilo meter)) m0
      g2 = do
         (x,y) <- gridSquare c2
         return $ GridOffset (x * (100 *~ kilo meter)) (y * (100 *~ kilo meter)) m0
      m0 = 0 *~ meter


-- | Find the nearest UK grid reference point to a specified position. The Int argument is the number of
-- digits precision, so 2 for a 4-figure reference and 3 for a 6-figure reference, although any value 
-- between 0 and 5 can be used (giving a 1 meter precision).
-- Altitude is ignored. If the result is outside the area defined by the two letter grid codes then
-- @Nothing@ is returned.
toUkGridReference :: Int -> GridPoint UkNationalGrid -> Maybe String
toUkGridReference n p
   | n < 0         = error "toUkGridReference: precision argument must not be negative."
   | n > 5         = error "toUkGridReference: precision argument must be <= 5."
   | not griddable = Nothing 
   | otherwise     = Just $ c1 : c2 : paddedShow dx ++ paddedShow dy
   where
      griddable = and [
         vx >= 0,
         vx < 2500000,
         vy >= 0,
         vy < 2500000]
      precI :: Integer
      precI = 100000 `div` (10 P.^ n)
      coord :: Length Double -> (Int, Int, Integer, Integer)
      coord v = (fromIntegral g1, fromIntegral g2, v1, v2)
         where
            -- v: Distance north or east from VV000000 (most southwesterly griddable point).
            -- v1: v rounded to grid precision, in integer meters.
            -- v2: Northing or easting of v1 within grid square.
            -- g2: Base of grid square below d, for second letter.
            -- g1: Base of larger (5x5) grid square, for first letter.
            v1 = round ((v /~ meter) P./ fromIntegral precI) P.* precI
            (g,  v2) = v1 `divMod` 100000
            (g1, g2) = g `divMod` 5
      (gx1,gx2,vx,dx) = coord $ eastings p + 1000 *~ kilo meter
      (gy1,gy2,vy,dy) = coord $ northings p + 500 *~ kilo meter
      c1 = gridSquare gx1 gy1
      c2 = gridSquare gx2 gy2
      gridSquare x y = letters ! (4 P.- y, x)
      letters :: Array (Int, Int) Char
      letters = listArray ((0,0),(4,4)) $ ['A'..'H'] ++ ['J'..'Z']
      paddedShow d = take n $ replicate (5 P.- length s) '0' ++ s
         where s = show d

      
   