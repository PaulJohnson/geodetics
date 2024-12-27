{-# LANGUAGE MultiParamTypeClasses #-}

-- | Distinguished coordinate systems for the United Kingdom.
module Geodetics.UK (
   OSGB36 (..),
   UkNationalGrid (..),
   ukGrid,
   fromUkGridReference,
   toUkGridReference
) where

import Control.Monad
import Data.Array
import Data.Char
import Geodetics.Geodetic
import Geodetics.Grid
import Geodetics.Ellipsoids
import Geodetics.TransverseMercator



-- | Ellipsoid definition for Great Britain. Airy 1830 offset from the centre of the Earth
-- and rotated slightly.
--
-- The Helmert parameters are from the Ordnance Survey document
-- \"A Guide to Coordinate Systems in Great Britain\", which notes that it
-- can be in error by as much as 5 meters and should not be used in applications
-- requiring greater accuracy.  A more precise conversion requires a large table
-- of corrections for historical inaccuracies in the triangulation of the UK.
data OSGB36 = OSGB36 deriving (Eq, Show)

instance Ellipsoid OSGB36 where
   majorRadius _ = 6377563.396
   flatR _ = 299.3249646
   helmert _ = Helmert {
      cX = 446.448, cY = (-125.157), cZ = 542.06,
      helmertScale = (-20.4894),
      rX = 0.1502 * arcsecond, rY = 0.247 * arcsecond, rZ = 0.8421 * arcsecond }

-- | The UK National Grid is a Transverse Mercator projection with a true origin at
-- 49 degrees North, 2 degrees West on OSGB36, and a false origin 400km West and 100 km North of
-- the true origin. The scale factor is defined as @10**(0.9998268 - 1)@.
data UkNationalGrid = UkNationalGrid deriving (Eq, Show)

instance GridClass UkNationalGrid OSGB36 where
   toGrid _ = unsafeGridCoerce UkNationalGrid . toGrid ukGrid
   fromGrid = fromGrid . unsafeGridCoerce ukGrid
   gridEllipsoid _ = OSGB36



ukTrueOrigin :: Geodetic OSGB36
ukTrueOrigin = Geodetic {
   latitude = 49 * degree,
   longitude = (-2) * degree,
   geoAlt = 0,
   ellipsoid = OSGB36
}

ukFalseOrigin :: GridOffset
ukFalseOrigin = GridOffset ((-400) * kilometer) (100 * kilometer) (0 * kilometer)


-- | Numerical definition of the UK national grid.
ukGrid :: GridTM OSGB36
ukGrid = mkGridTM ukTrueOrigin ukFalseOrigin
   (10 ** (0.9998268 - 1))


-- | Size of a UK letter-pair grid square.
ukGridSquare :: Double
ukGridSquare = 100 * kilometer


-- | Convert a grid reference to a position, if the reference is valid.
-- This returns the position of the south-west corner of the nominated
-- grid square and an offset to its centre. Altitude is set to zero.
fromUkGridReference :: String -> Maybe (GridPoint UkNationalGrid, GridOffset)
fromUkGridReference str =
      case str of
         c1:c2:ds -> do
            let n = length ds
            guard $ even n
            let (dsE, dsN) = splitAt (n `div` 2) ds
            (east, sq) <- fromGridDigits ukGridSquare dsE
            (north, _) <- fromGridDigits ukGridSquare dsN
            base <- fromUkGridLetters c1 c2
            let half = sq / 2
            return (applyOffset (GridOffset east north 0) base,
                  GridOffset half half 0)
         _ -> Nothing




-- | The south west corner of the nominated grid square, if it is a legal square.
-- This function works for all pairs of letters except 'I' (as that is not used).
-- In practice only those pairs covering the UK are actually considered meaningful.
fromUkGridLetters :: Char -> Char -> Maybe (GridPoint UkNationalGrid)
fromUkGridLetters c1 c2 = applyOffset <$> (mappend <$> g1 <*> g2) <*> letterOrigin
   where
      letterOrigin = Just $ GridPoint ((-1000) * kilometer) ((-500) * kilometer) m0 UkNationalGrid
      gridIndex c =
         if inRange ('A', 'H') c then Just $ ord c - ord 'A'  -- 'I' is not used.
         else if inRange ('J', 'Z') c then Just $ ord c - ord 'B'
         else Nothing
      gridSquare c = do -- Maybe monad
         g <- gridIndex c
         let (y,x) = g `divMod` 5
         return (fromIntegral x, 4 - fromIntegral y)
      g1 = do
         (x,y) <- gridSquare c1
         return $ GridOffset (x * (500 * kilometer)) (y * (500 * kilometer)) m0
      g2 = do
         (x,y) <- gridSquare c2
         return $ GridOffset (x * (100 * kilometer)) (y * (100 * kilometer)) m0
      m0 = 0


-- | Find the nearest UK grid reference point to a specified position. The Int argument is the number of
-- digits precision, so 2 for a 4-figure reference and 3 for a 6-figure reference, although any value
-- between 0 and 5 can be used (giving a 1 meter precision).
-- Altitude is ignored. If the result is outside the area defined by the two letter grid codes then
-- @Nothing@ is returned.
toUkGridReference :: Int -> GridPoint UkNationalGrid -> Maybe String
toUkGridReference n p
   | n < 0         = error "toUkGridReference: precision argument must not be negative."
   | otherwise     = do
      (gx, strEast) <- toGridDigits ukGridSquare n $ eastings p + 1000 * kilometer
      (gy, strNorth) <- toGridDigits ukGridSquare n $ northings p + 500 * kilometer
      let (gx1, gx2) = fromIntegral gx `divMod` 5
          (gy1, gy2) = fromIntegral gy `divMod` 5
      guard (gx1 < 5 && gy1 < 5)
      let c1 = gridSquare gx1 gy1
          c2 = gridSquare gx2 gy2
      return $ c1 : c2 : strEast ++ strNorth
   where
      gridSquare x y = letters ! (4 - y, x)
      letters :: Array (Int, Int) Char
      letters = listArray ((0,0),(4,4)) $ ['A'..'H'] ++ ['J'..'Z']
