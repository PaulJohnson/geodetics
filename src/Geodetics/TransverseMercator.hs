{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}

module Geodetics.TransverseMercator(
   GridTM (trueOrigin, falseOrigin, gridScale),
   mkGridTM
) where

import Data.Function
import Data.Monoid
import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid
import Numeric.Units.Dimensional.Prelude
import Prelude ()

-- | A Transverse Mercator projection gives an approximate mapping of the ellipsoid on to a 2-D grid. It models
-- a sheet curved around the ellipsoid so that it touches it at one north-south line (hence making it part of
-- a slightly elliptical cylinder).
data GridTM e = GridTM {
   trueOrigin :: Geodetic e,
      -- ^ A point on the line where the projection touches the ellipsoid (altitude is ignored).
   falseOrigin :: GridOffset,
      -- ^ The grid position of the true origin. Used to avoid negative coordinates over 
      -- the area of interest. The altitude gives a vertical offset from the ellipsoid.
   gridScale :: Dimensionless Double,
      -- ^ A scaling factor that balances the distortion between the east & west edges and the middle 
      -- of the projection.
      
   -- Remaining elements are memoised parameters computed from the ellipsoid underlying the true origin.
   gridN1, gridN2, gridN3, gridN4 :: Dimensionless Double
} deriving (Show)


-- | Create a Transverse Mercator grid.
mkGridTM :: (Ellipsoid e) => 
   Geodetic e               -- ^ True origin.
   -> GridOffset            -- ^ Vector from true origin to false origin.
   -> Dimensionless Double  -- ^ Scale factor.
   -> GridTM e
mkGridTM origin offset sf =
   GridTM {trueOrigin = origin,
           falseOrigin = offset,
           gridScale = sf,
           gridN1 = _1 + n + (_5/_4) * n^pos2 + (_5/_4) * n^pos3,
           gridN2 = _3 * n + _3 * n^pos2 + ((21*~one)/_8) * n^pos3,
           gridN3 = ((15*~one)/_8) * (n^pos2 + n^pos3),
           gridN4 = ((35*~one)/(24*~one)) * n^pos3
        }
    where 
       f = flattening $ ellipsoid origin
       n = f / (_2-f)  -- Equivalent to (a-b)/(a+b) where b = (1-f)*a




-- | Equation C3 from reference [1].
m :: (Ellipsoid e) => GridTM e -> Dimensionless Double -> Length Double
m grid lat = bF0 * (gridN1 grid * dLat 
                    - gridN2 grid * sin dLat * cos sLat
                    + gridN3 grid * sin (_2 * dLat) * cos (_2 * sLat) 
                    - gridN4 grid * sin (_3 * dLat) * cos (_3 * sLat))
   where
      dLat = lat - latitude (trueOrigin grid)
      sLat = lat + latitude (trueOrigin grid)
      bF0 = minorRadius (gridEllipsoid grid) * gridScale grid


instance (Ellipsoid e) => GridClass (GridTM e) e where
   fromGrid p = Geodetic
      (lat' - east' ^ pos2 * tanLat / (_2 * rho * v)  -- Term VII
            + east' ^ pos4 * (tanLat / ((24 *~ one) * rho * v ^ pos3)) 
                           * (_5 + _3 * tanLat ^ pos2 + eta2 - _9 * tanLat ^ pos2 * eta2)  -- Term VIII
            - east' * east' ^ pos5 * (tanLat / ((720 *~ one) * rho * v ^ pos5))
                           * (61 *~ one + (90 *~ one) * tanLat ^ pos2 + (45 *~ one) * tanLat ^ pos4)) -- Term IX
      (longitude (trueOrigin grid) 
            + east' / (cosLat * v)  -- Term X
            - (east' ^ pos3 / (_6 * cosLat * v ^ pos3)) * (v / rho + _2 * tanLat ^ pos2)  -- Term XI
            + (east' ^ pos5 / ((120 *~ one) * cosLat * v ^ pos5)) 
                 * (_5 + (28 *~ one) * tanLat ^ pos2  + (24 *~ one) * tanLat ^ pos4)  -- Term XII
            - (east' ^ pos5 * east' ^ pos2 / ((5040 *~ one) * cosLat * v * v * v ^ pos5))
                 * ((61 *~ one) + (662 *~ one) * tanLat ^ pos2 + (1320 *~ one) * tanLat ^ pos4 + (720 *~ one) * tanLat * tanLat ^ pos5)) -- Term XIIa
     (0 *~ meter) (gridEllipsoid grid)
            
            
      where
         GridPoint east' north' _ _ = falseOrigin grid `applyOffset` p
         lat' = fst $ head $ dropWhile ((> 0.01 *~ milli meter) . snd) 
               $ tail $ iterate next (latitude $ trueOrigin grid, 1 *~ meter) 
            where
               next (phi, _) = let delta = north' - m grid phi in (phi + delta / aF0, delta) 
               -- head and tail are safe because iterate returns an infinite list.
          
         sinLat = sin lat'
         cosLat = cos lat'
         tanLat = tan lat'
         sinLat2 = sinLat ^ pos2
         v = aF0 / sqrt (_1 - e2 * sinLat2)
         rho = aF0 * (_1 - e2) * (_1 - e2 * sinLat2) ** ((-1.5) *~ one)
         eta2 = v / rho - _1
               
               
         aF0 = majorRadius (gridEllipsoid grid) * gridScale grid
         e2 = eccentricity2 $ gridEllipsoid grid
         grid = gridBasis p
         
   toGrid grid geo = applyOffset (off  `mappend` (offsetNegate $ falseOrigin grid)) $ 
                     GridPoint _0 _0 _0 grid
      where
         v = aF0 / sqrt (_1 - e2 * sinLat2)
         rho = aF0 * (_1 - e2) * (_1 - e2 * sinLat2) ** ((-1.5) *~ one)
         eta2 = v / rho - _1
         off = GridOffset
                  (dLong * term_IV
                   + dLong ^ pos3 * term_V
                   + dLong ^ pos5 * term_VI)
                  (m grid lat + dLong ^ pos2 * term_II
                     + dLong ^ pos4 * term_III 
                     + dLong * dLong ^ pos5 * term_IIIa)
                  (0 *~ meter)
         -- Terms defined in [1].
         term_II   = (v/_2) * sinLat * cosLat
         term_III  = (v/(24*~one)) * sinLat * cosLat ^ pos3 
                     * (_5 - tanLat ^ pos2 + _9 * eta2)
         term_IIIa = (v/(720*~one)) * sinLat * cosLat ^ pos5 
                     * ((61 *~ one) - (58 *~ one) * tanLat ^ pos2 + tanLat ^ pos4)
         term_IV   = v * cosLat
         term_V    = (v/_6) * cosLat ^ pos3 * (v/rho - tanLat ^ pos2)
         term_VI   = (v/(120*~one)) * cosLat ^ pos5 
                     * (_5 - (18*~one) * tanLat ^ pos2 
                              + tanLat ^ pos4 + (14*~one) * eta2
                              - (58*~one) * tanLat ^ pos2 * eta2)
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
         sinLat2 = sinLat ^ pos2
         aF0 = (majorRadius $ gridEllipsoid grid) * gridScale grid
         e2 = eccentricity2 $ gridEllipsoid grid
         
   gridEllipsoid = ellipsoid . trueOrigin
