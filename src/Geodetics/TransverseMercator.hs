{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}

module Geodetics.TransverseMercator(
   GridTM (trueOrigin, falseOrigin, gridScale),
   mkGridTM
) where

import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid

import qualified Data.Stream as Stream

-- | A Transverse Mercator projection gives an approximate mapping of the ellipsoid on to a 2-D grid. It models
-- a sheet curved around the ellipsoid so that it touches it at one north-south line (hence making it part of
-- a slightly elliptical cylinder).
--
-- The calculations here are based on \"Transverse Mercator Projection: Constants, Formulae and Methods\"
-- by the Ordnance Survey, March 1983.
-- Retrieved from http://www.threelittlemaids.co.uk/magdec/transverse_mercator_projection.pdf
data GridTM e = GridTM {
   trueOrigin :: Geodetic e,
      -- ^ A point on the line where the projection touches the ellipsoid (altitude is ignored).
   falseOrigin :: GridOffset,
      -- ^ The grid position of the true origin. Used to avoid negative coordinates over
      -- the area of interest. The altitude gives a vertical offset from the ellipsoid.
   gridScale :: Double,
      -- ^ A scaling factor that balances the distortion between the east & west edges and the middle
      -- of the projection.

   -- Remaining elements are memoised parameters computed from the ellipsoid underlying the true origin.
   gridN1, gridN2, gridN3, gridN4 :: Double
} deriving (Show)


-- | Create a Transverse Mercator grid.
mkGridTM :: (Ellipsoid e) =>
   Geodetic e               -- ^ True origin.
   -> GridOffset            -- ^ Vector from true origin to false origin.
   -> Double                -- ^ Scale factor.
   -> GridTM e
mkGridTM origin offset sf =
   GridTM {trueOrigin = origin,
           falseOrigin = offset,
           gridScale = sf,
           gridN1 = 1 + n + (5/4) * n**2 + (5/4) * n**3,
           gridN2 = 3 * n + 3 * n**2 + (21/8) * n**3,
           gridN3 = (15/8) * (n**2 + n**3),
           gridN4 = (35/24) * n**3
        }
    where
       f = flattening $ ellipsoid origin
       n = f / (2-f)  -- Equivalent to (a-b)/(a+b) where b = (1-f)*a


-- | Equation C3 from reference [1].
m :: (Ellipsoid e) => GridTM e -> Double -> Double
m grid lat = bF0 * (gridN1 grid * dLat
                    - gridN2 grid * sin dLat * cos sLat
                    + gridN3 grid * sin (2 * dLat) * cos (2 * sLat)
                    - gridN4 grid * sin (3 * dLat) * cos (3 * sLat))
   where
      dLat = lat - latitude (trueOrigin grid)
      sLat = lat + latitude (trueOrigin grid)
      bF0 = minorRadius (gridEllipsoid grid) * gridScale grid


instance (Ellipsoid e) => GridClass (GridTM e) e where
   fromGrid p = Geodetic
      (lat' - east' ** 2 * tanLat / (2 * rho * v)  -- Term VII
            + east' ** 4 * (tanLat / (24 * rho * v ** 3))
                           * (5 + 3 * tanLat ** 2 + eta2 - 9 * tanLat ** 2 * eta2)  -- Term VIII
            - east' ** 6 * (tanLat / (720 * rho * v ** 5))
                           * (61 + 90 * tanLat ** 2 + 45 * tanLat ** 4)) -- Term IX
      (longitude (trueOrigin grid)
            + east' / (cosLat * v)  -- Term X
            - (east' ** 3 / (6 * cosLat * v ** 3)) * (v / rho + 2 * tanLat ** 2)  -- Term XI
            + (east' ** 5 / (120 * cosLat * v ** 5))
                 * (5 + 28 * tanLat ** 2  + 24 * tanLat ** 4)  -- Term XII
            - (east' ** 5 * east' ** 2 / (5040 * cosLat * v ** 7))
                 * (61 + 662 * tanLat ** 2 + 1320 * tanLat ** 4 + 720 * tanLat ** 6)) -- Term XIIa
     0
     (gridEllipsoid grid)

      where
         GridPoint east' north' _ _ = falseOrigin grid `applyOffset` p
         lat' = fst $ Stream.head $ Stream.dropWhile ((> 1e-5) . snd)
               $ Stream.tail $ Stream.iterate next (latitude $ trueOrigin grid, 1)
            where
               next (phi, _) = let delta = north' - m grid phi in (phi + delta / aF0, delta)

         sinLat = sin lat'
         cosLat = cos lat'
         tanLat = tan lat'
         sinLat2 = sinLat * sinLat
         v = aF0 / sqrt (1 - e2 * sinLat2)
         rho = aF0 * (1 - e2) * (1 - e2 * sinLat2) ** (-1.5)
         eta2 = v / rho - 1

         aF0 = majorRadius (gridEllipsoid grid) * gridScale grid
         e2 = eccentricity2 $ gridEllipsoid grid
         grid = gridBasis p

   toGrid grid geo = applyOffset (off  `mappend` offsetNegate (falseOrigin grid)) $
                     GridPoint 0 0 0 grid
      where
         v = aF0 / sqrt (1 - e2 * sinLat2)
         rho = aF0 * (1 - e2) * (1 - e2 * sinLat2) ** (-1.5)
         eta2 = v / rho - 1
         off = GridOffset
                  (dLong * term_IV
                   + dLong ** 3 * term_V
                   + dLong ** 5 * term_VI)
                  (m grid lat + dLong ** 2 * term_II
                     + dLong ** 4 * term_III
                     + dLong ** 6 * term_IIIa)
                  0
         -- Terms defined in [1].
         term_II   = (v/2) * sinLat * cosLat
         term_III  = (v/24) * sinLat * cosLat ** 3
                     * (5 - tanLat ** 2 + 9 * eta2)
         term_IIIa = (v/720) * sinLat * cosLat ** 5
                     * (61 - 58 * tanLat ** 2 + tanLat ** 4)
         term_IV   = v * cosLat
         term_V    = (v/6) * cosLat ** 3 * (v/rho - tanLat ** 2)
         term_VI   = (v/120) * cosLat ** 5
                     * (5 - 18 * tanLat ** 2
                              + tanLat ** 4 + 14 * eta2
                              - 58 * tanLat ** 2 * eta2)
         {-
         -- Trace message for debugging. Uncomment this code for easy access to intermediate values.
         traceMsg = concat [
            "v    = ", show v, "\n",
            "rho  = ", show rho, "\n",
            "eta2 = ", show eta2, "\n",
            "M    = ", show $ m grid lat, "\n",
            "I    = ", show $ m grid lat + deltaNorth (falseOrigin grid), "\n",
            "II   = ", show term_II, "\n",
            "III  = ", show term_III, "\n",
            "IIIa = ", show term_IIIa, "\n",
            "IV   = ", show term_IV, "\n",
            "V    = ", show term_V, "\n",
            "VI   = ", show term_VI, "\n"]
         -}
         -- Common subexpressions
         lat = latitude geo
         long = longitude geo
         dLong = long - longitude (trueOrigin grid)
         sinLat = sin lat
         cosLat = cos lat
         tanLat = tan lat
         sinLat2 = sinLat * sinLat
         aF0 = majorRadius (gridEllipsoid grid) * gridScale grid
         e2 = eccentricity2 $ gridEllipsoid grid

   gridEllipsoid = ellipsoid . trueOrigin
