{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Geodetics.PolarStereographic (
   Pole (..),
   PolarStereographic (trueOrigin, falseOrigin, polarEllipsoid, gridScale),
   mkGridPolarStereographic
) where

import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid

-- | Polar stereographic grids are defined for true origins at the north and south poles.
data Pole = NorthPole | SouthPole deriving (Show, Ord, Eq, Enum, Bounded)


{- | Polar Stereographic Grids

Formulae are taken from
/The Universal Grids: Univerersal Transverse Mercator (UTM) and Universal Polar Stereographic (UPS)/
DMA Technical Manual 8358.2, Defense Mapping Agency, Fairfax, VA. https://apps.dtic.mil/sti/tr/pdf/ADA266497.pdf
-}
data PolarStereographic e = PolarStereographic {
   trueOrigin :: Pole,
   falseOrigin :: GridOffset,
      -- ^ The negation of the grid position of the true origin. Used to avoid negative coordinates over the area
      -- of interest. The altitude gives a vertical offset from the ellipsoid.
   polarEllipsoid :: e,
      -- ^ The ellipsoid for the projection. Arguments passed to `toGrid` *must* use this ellipsoid.
      -- The type system cannot verify this for `LocalEllipsoid`.
   gridScale :: Double,
      -- ^ The scaling factor applied at the pole. This balances the distortion between the center
      -- and edges of the projection.

   -- Remaining elements are memoised parameters computed from the ellipsoid.
   gridA, gridB, gridC, gridD, gridC0 :: !Double
} deriving (Show)

instance (Ellipsoid e) => GridClass (PolarStereographic e) e where
   fromGrid p = Geodetic lat long (altGP p) (polarEllipsoid gb)
      where
         gridZero = GridPoint 0 0 0 gb
         gb = gridBasis p
         p' = gridZero `gridOffset` (falseOrigin gb `applyOffset` p)
         radius = offsetDistance p'
         isoColat = 2 * atan (radius / (gridScale gb * gridC0 gb))
         isoLat = pi/2 - isoColat
         lat1 = isoLat +
            gridA gb * sin (2*isoLat) +
            gridB gb * sin (4*isoLat) +
            gridC gb * sin (6*isoLat) +
            gridD gb * sin (8*isoLat)
         lat = case trueOrigin gb of
            NorthPole -> lat1
            SouthPole -> negate lat1
         long = case trueOrigin gb of
            NorthPole -> offsetBearing p' { deltaNorth = negate $ deltaNorth p'}
            SouthPole -> offsetBearing p'

   toGrid r geo = offsetNegate (falseOrigin r) `applyOffset` GridPoint east north 0 r
      where
         absLat = abs $ latitude geo
         e = sqrt (eccentricity2 $ polarEllipsoid r)
         eSinLat = e * sin absLat
         tz2 = ((1 + eSinLat)/(1-eSinLat))**(e/2) * tan (pi/4 - absLat / 2)
         radius = gridScale r * gridC0 r * tz2
         north = case trueOrigin r of
            NorthPole -> negate $ radius * cos (longitude geo)
            SouthPole -> radius * cos (longitude geo)
         east = radius * sin (longitude geo)
   gridEllipsoid = polarEllipsoid


mkGridPolarStereographic :: (Ellipsoid e) =>
   Pole          -- ^ True origin at north or south pole.
   -> e          -- ^ The ellipsoid used for the projection.
   -> GridOffset -- ^ Vector from true origin to the false origin.
   -> Double     -- ^ Scale factor.
   -> PolarStereographic e
mkGridPolarStereographic pole ellip offset scale =
   PolarStereographic {
      trueOrigin = pole,
      falseOrigin = offset,
      polarEllipsoid = ellip,
      gridScale = scale,
      gridA = e2/2 + (5/24)*e4 + e6/12 + (13/360)*e8,
      gridB = (7/48)*e4 + (29/240)*e6 + (811/11520)*e8,
      gridC = (7/120)*e6 + (81/1120)*e8,
      gridD = (4279/161280)*e8,
      gridC0 = (2 * majorRadius ellip / sqrt (1 - e2)) * ((1-e1)/(1+e1))**(e1/2)
   }
   where
      e1 = sqrt $ eccentricity2 ellip
      e2 = eccentricity2 ellip
      e4 = e2^_2
      e6 = e2^_3
      e8 = e2^_4
