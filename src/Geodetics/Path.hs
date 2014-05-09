-- | The implementation assumes IEEE 754 arithmetic.

module Geodetics.Path where

import Control.Monad
import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P

import Debug.Trace  -- Uncomment if using tracing to debug @intersect@.

-- | Lower and upper exclusive bounds within which a path is valid. 
type PathValidity = (Length Double, Length Double)

-- | A path is a parametric function of distance along the path. The result is the
-- position, and the bearing and elevation along the path at that point.
--
-- A well-behaved path should be continuous for all distances within its validity range,
-- and make the physical distance along the path from the origin to the result be 
-- approximately equal to the argument. Outside its validity the path function may
-- return anything or bottom.
data Path e = Path {
      pathFunc :: Length Double -> (Geodetic e, Angle Double, Angle Double),
      pathValidity :: PathValidity
   }
   
-- | Convenience value for paths that are valid for all distances.
alwaysValid :: PathValidity
alwaysValid = (negate inf, inf) where
   inf = (1.0 *~ meter) / (0 *~ one)


pathValidAt :: Path e -> Length Double -> Bool
pathValidAt path d = d > x1 && d < x2
   where (x1,x2) = pathValidity path


-- | Find where a path meets a given condition using bisection. This is not the
-- most efficient algorithm, but it will always produce a result if there is one
-- within the initial bounds. If there is more than one result then an arbitrary
-- one will be returned.
-- 
-- The initial bounds must return one GT or EQ value and one LT or EQ value. If they
-- do not then @Nothing@ is returned.
bisect :: 
   Path e 
   -> (Geodetic e -> Ordering)        -- ^ Evaluation function.
   -> Length Double                   -- ^ Required accuracy in terms of distance along the path.
   -> Length Double -> Length Double  -- ^ Initial bounds.
   -> Maybe (Length Double)
bisect path f t b1 b2 = do
      guard $ pathValidAt path b1
      guard $ pathValidAt path b2
      let r = pairResults b1 b2
      guard $ hasRoot r
      bisect1 $ sortPair r
   where
      f' d = let (p, _, _) = pathFunc path d in f p 
      pairResults d1 d2 = ((d1, f' d1), (d2, f' d2))
      hasRoot (v1, v2) = snd v1 <= EQ && EQ <= snd v2
      sortPair (v1, v2) = if snd v1 <= snd v2 then (v1, v2) else (v2, v1)
      bisect1 ((d1, r1), (d2, r2)) =
         let d3 = (d1 + d2) / _2
             r3 = f' d3
             c1 = ((d1, r1), (d3, r3))
             c2 = ((d3, r3), (d2, r2))
         in if abs (d1 - d2) <= t 
            then return d3
            else bisect1 $ if hasRoot c1 then c1 else c2


-- | Try to find the intersection point of two ground paths (i.e. ignoring 
-- altitude). Returns the distance of the intersection point along each path
-- using a modified Newton-Raphson method. If the two paths are 
-- fairly straight and not close to parallel then this will converge rapidly.
--
-- The algorithm projects great-circle paths forwards using the bearing at the 
-- estimate to find the estimated intersection, and then uses the distances to this 
-- intersection as the next estimate.
--
-- If the estimates depart from the path validity then @Nothing@ is returned.
intersect :: (Ellipsoid e) =>
   Length Double -> Length Double     -- ^ Starting estimates.
   -> Length Double                   -- ^ Required accuracy.
   -> Int                             -- ^ Iteration limit.  
   -> Path e -> Path e                -- ^ Paths to intersect.
   -> Maybe (Length Double, Length Double)
intersect d1 d2 accuracy n path1 path2
   | not $ pathValidAt path1 d1  = Nothing
   | not $ pathValidAt path2 d2  = Nothing
   | n <= 0                      = Nothing
   | abs (atan2 sinBase cosBase) * r <= accuracy     = Just (d1, d2)
   | otherwise = 
      trace ("From " ++ show (pt1, pt2)
            ++ ", base = " ++ show ((r * atan2 sinBase cosBase) /~ kilo meter)
            ++ ", new deltas " ++ show (d1' /~ kilo meter, d2' /~ kilo meter) 
            ++ ", bearings = " ++ show (b1 /~ degree, b2 /~ degree)
            ++ ", angles = " ++ show (theta1 /~ degree, theta2 /~ degree)) $
         intersect (d1 + d1') (d2 + d2') accuracy (pred n) path1 path2
   where 
      (pt1, h1, _) = pathFunc path1 d1
      (pt2, h2, _) = pathFunc path2 d2      
      theta0 = longitude pt2 - longitude pt1   
      sinTheta0 = sin theta0
      cosTheta0 = cos theta0
      cosLat1 = cos $ pi/_2 - latitude pt1  -- Distance from north pole on unit sphere
      cosLat2 = cos $ pi/_2 - latitude pt2  -- Ignoring trig identities for clarity.
      sinLat1 = sin $ pi/_2 - latitude pt1
      sinLat2 = sin $ pi/_2 - latitude pt2
      cosBase = cosLat1 * cosLat2 + sinLat1 * sinLat2 * cosTheta0  -- Cosine Rule
      sinBase = sqrt(_1 - cosBase * cosBase)
      tb = sinTheta0 / sinBase
      sinB1 = sinLat2 * tb           -- Bearing of pt2 from pt1, by Sine Rule
      sinB2 = negate $ sinLat1 * tb  -- and vice versa.  Negated because its the other side.   
      cosB1 = (cosLat2 - cosLat1 * cosBase) / (sinLat1 * sinBase)  -- And by Cosine Rule
      cosB2 = (cosLat1 - cosLat2 * cosBase) / (sinLat2 * sinBase)
      b1 = atan2 sinB1 cosB1         -- Gives correct bearing for all quadrants.
      b2 = atan2 sinB2 cosB2
      theta1 = b1 - h1
      theta2 = h2 - b2 -- Other side of triangle, so angles are the other way round.
      d1' = r * atan (sinBase / (cosBase * cos theta1 + sin theta1 / tan theta2))
      d2' = r * atan (sinBase / (cosBase * cos theta2 + sin theta2 / tan theta1))  -- Cotan rule.
      r = majorRadius $ ellipsoid pt1
      
{- Note on derivation

A spherical approximation is used to get the approximate distances to pt3,
the nearest intersection of the great circles that lie tangent to the 
paths at pt1 and pt2. All the arithmetic is performed on the unit sphere
except for d1 and d2, the distances along the paths.

pt1 and pt2 form a triangle with the north pole. The lines of longitude
from the pole to pt1 and pt2 are the known sides of the triangle, and the angle
theta0 is the difference in the longitudes of pt1 and pt2.

By the cosine rule we get the baseline from pt1 to pt2, and the
interior angles of this triangle give us the bearings (b1, b2) of the end points
of the baseline.

Hence we now have a triangle (pt1,pt2,pt3) where we know the distance pt1->pt2
and the interior angles theta1 & theta2 at pt1 & pt2 respectively. The cotangent
rule now gives us the distances from pt1->pt3 and pt2->pt3 (d1', d2').

The hairiest aspect of this is not the spherical trig, its making sure that
it does the Right Thing for every quadrant of every angle.

-}


-- | A ray from a point heading in a straight line in 3 dimensions. 
rayPath :: (Ellipsoid e) => 
   Geodetic e          -- ^ Start point.
   -> Angle Double     -- ^ Bearing.
   -> Angle Double     -- ^ Elevation.
   -> Path e
rayPath pt1 bearing elevation = Path ray alwaysValid
   where
      ray distance = (Geodetic lat long alt (ellipsoid pt1), bearing2, elevation2)
         where
            pt2' = pt1' `add3` (delta `scale3` distance)      -- ECEF of result point.
            (lat,long,alt) = earthToGeo (ellipsoid pt1) pt2'  -- Geodetic of result point.
            (dE,dN,dU) = transform3 (trans3 $ ecefMatrix lat long) delta  -- Direction of ray at result point.
            elevation2 = asin dU
            bearing2 = if dE == _0 && dN == _0 then bearing else atan2 dE dN  -- Allow for vertical elevation.
            
      ecefMatrix lat long =   -- Transform matrix for vectors from (East, North, Up) to (X,Y,Z).
         ((negate sinLong, negate cosLong*sinLat, cosLong*cosLat),
              --    East X      North X               Up X
          (       cosLong, negate sinLong*sinLat, sinLong*cosLat),
              --    East Y      North Y               Up Y
          (  _0           ,      cosLat         , sinLat))
              --    East Z      North Z               Up Z
         where
            sinLong = sin long
            cosLong = cos long
            sinLat = sin lat
            cosLat = cos lat
      
      direction = (sinB*cosE, cosB*cosE, sinE)  -- Direction of ray in ENU
      delta = transform3 (ecefMatrix (latitude pt1) (longitude pt1)) direction  -- Convert to ECEF
      pt1' = geoToEarth pt1    -- ECEF of origin point.
      sinB = sin bearing
      cosB = cos bearing
      sinE = sin elevation
      cosE = cos elevation
      
-- | Rhumb line: path following a constant course. Also known as a loxodrome.
--
-- The valid range stops a few arc-minutes short of the poles. 
--
-- Based on *Practical Sailing Formulas for Rhumb-Line Tracks on an Oblate Earth* 
-- by G.H. Kaplan, U.S. Naval Observatory. Except for points close to the poles 
-- the approximation is accurate to within a few meters over 1000km.
rhumbPath :: (Ellipsoid e) =>
   Geodetic e            -- ^ Start point.
   -> Angle Double       -- ^ Course.
   -> Path e
rhumbPath pt course = Path rhumb validity
   where
      rhumb distance = (Geodetic lat (properAngle lon) (0 *~ meter) (ellipsoid pt), course, _0)
         where
            lat' = lat0 + distance * cosC / m0   -- Kaplan Eq 13.
            lat = lat0 + (m0 / (a*(_1-e2))) * ((_1-_3*e2/_4)*(lat'-lat0)
                                              + (_3*e2/_8)*(sin (_2*lat') - sin (_2*lat0)))
            lon | abs cosC > 1e-7 *~ one 
                     = lon0 + tanC * (q lat - q0)     -- Kaplan Eq 16.
                | otherwise
                     = lon0 + distance * sinC / latitudeRadius (ellipsoid pt) ((lat0 + lat')/_2)
      validity
         | cosC > _0  = ((negate pi/_2 - latitude pt) * b / cosC, (pi/_2 - latitude pt) * b / cosC)
         | otherwise  = ((pi/_2 - latitude pt) * b / cosC, (negate pi/_2 - latitude pt) * b / cosC)
      q0 = q lat0
      q phi = log (tan (pi/_4+phi/_2)) + e * log ((_1-eSinPhi)/(_1+eSinPhi)) / _2
         where                                -- Factor out expression from Eq 16 of Kaplan
            eSinPhi = e * sin phi
      sinC = sin course
      cosC = cos course
      tanC = tan course
      lat0 = latitude pt
      lon0 = longitude pt
      e2 = eccentricity2 $ ellipsoid pt
      e = sqrt e2
      m0 = meridianRadius (ellipsoid pt) lat0
      a = majorRadius $ ellipsoid pt
      b = minorRadius $ ellipsoid pt
   

-- | A path following the line of latitude around the Earth eastwards.
--
-- This is equivalent to @rhumbPath pt (pi/2)@
latitudePath :: (Ellipsoid e) =>
   Geodetic e           -- ^ Start point.
   -> Path e
latitudePath pt = Path line alwaysValid
   where
      line distance = (pt2, pi/_2, _0) 
         where
            pt2 = Geodetic 
               (latitude pt) 
               (longitude pt + distance / r)
               (0 *~ meter)
               (ellipsoid pt)
      r = latitudeRadius (ellipsoid pt) (latitude pt)


-- | A path from the specified point to the North Pole. Use negative distances
-- for the southward path.
--
-- This is equivalent to @rhumbPath pt _0@
longitudePath :: (Ellipsoid e) =>
   Geodetic e    -- ^ Start point.
   -> Path e
longitudePath pt = rhumbPath pt _0
      


