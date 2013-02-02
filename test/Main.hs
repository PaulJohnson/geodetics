
module Main where

import Data.Maybe
import Data.Monoid
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P
import Test.Framework (Test, defaultMainWithOpts, testGroup)
import Test.Framework.Options (TestOptions, TestOptions'(..))
import Test.Framework.Runners.Options (RunnerOptions, RunnerOptions'(..))
import Test.Framework.Providers.HUnit
import Test.Framework.Providers.QuickCheck2 (testProperty)
-- import Test.QuickCheck
import qualified Test.HUnit as HU 

import ArbitraryInstances
import Geodetics.Coordinates
import Geodetics.Ellipsoids
import Geodetics.UK


main :: IO ()
main = do
   let empty_test_opts = mempty :: TestOptions
   let my_test_opts = empty_test_opts {
     topt_maximum_generated_tests = Just 1000
   }

   let empty_runner_opts = mempty :: RunnerOptions
   let my_runner_opts = empty_runner_opts {
     ropt_test_options = Just my_test_opts
   }

   defaultMainWithOpts tests my_runner_opts

tests :: [Test]
tests = [
   testGroup "Geodetic" [
      testProperty "WGS84 and back" prop_WGS84_and_back,
      testGroup "UK Points" $ map pointTest ukPoints],
   testGroup "Grid" [
      testProperty "Grid Offset 1" prop_offset1,
      testProperty "Grid Offset 2" prop_offset2,
      testProperty "Grid Offset 3" prop_offset3,
      testProperty "Grid 1" prop_grid1 ],
   testGroup "UK" [
      testProperty "UK Grid 1" prop_ukGrid1,
      testGroup "UK Grid 2" $ map ukGridTest2 ukSampleGrid,
      testGroup "UK Grid 3" $ map ukGridTest3 ukSampleGrid,
      testGroup "UK Grid 4" $ map ukGridTest4 ukSampleGrid,
      testGroup "UK Grid 5" $ map ukGridTest5 ukSampleGrid
      ]
   ]


-- | The positions are within 30 cm.
samePlace :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Bool
samePlace p1 p2 = geometricalDistance p1 p2 < 0.3 *~ meter


-- | The positions are within 10 m.
closeEnough :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Bool
closeEnough p1 p2 = geometricalDistance p1 p2 < 10 *~ meter


-- | The grid positions are within 1mm
sameGrid :: (GridClass r e) => GridPoint r -> GridPoint r -> Bool
sameGrid p1 p2 = check eastings && check northings && check altitude
   where check f = f p1 - f p2 < 1 *~ milli meter


-- | Grid offsets are within 1mm.
sameOffset :: GridOffset -> GridOffset -> Bool
sameOffset go1 go2 = check deltaNorth && check deltaEast && check deltaAltitude
   where check f = f go1 - f go2 < 1 *~ milli meter


-- | Degrees, minutes and seconds into radians. 
dms :: Int -> Int -> Double -> Dimensionless Double
dms d m s = fromIntegral d *~ degree + fromIntegral m *~ arcminute + s *~ arcsecond

-- | Round-trip is identity (approximately)
prop_WGS84_and_back :: Geodetic LocalEllipsoid -> Bool
prop_WGS84_and_back p = samePlace p $ toLocal (ellipsoid p) $ toWGS84 p


-- | Sample points for UK tests. The oracle for these values is the script at 
-- http://www.movable-type.co.uk/scripts/latlong-convert-coords.html, which uses
-- the same Helmert transform as this library. Hence the results should match to within 30 cm.
ukPoints :: [(String, Geodetic WGS84, Geodetic OSGB36)]
ukPoints = [
   ("Greenwich",        Geodetic (dms 51 28 40.86) (dms 0 0 (-5.83)) m0 WGS84, 
                        Geodetic (dms 51 28 39.00) (dms 0 0 0) m0 OSGB36),
   ("Edinburgh Castle", Geodetic (dms 55 56 56.30) (dms (-3) (-12) (-2.73)) m0 WGS84, 
                        Geodetic (dms 55 56 56.51) (dms (-3) (-11) (-57.61)) m0 OSGB36),
   ("Lands End",        Geodetic (dms 50 03 56.68) (dms (-5) (-42) (-51.20)) m0 WGS84,
                        Geodetic (dms 50 03 54.51) (dms (-5) (-42) (-47.87)) m0 OSGB36),
   ("Gt. Yarmouth Pier",Geodetic (dms 52 36 29.33) (dms 1 44 27.79) m0 WGS84,
                        Geodetic (dms 52 36 27.84) (dms 1 44 34.52) m0 OSGB36),
   ("Stanhope",         Geodetic (dms 54 44 49.08) (dms (-2) 0 (-19.89)) m0 WGS84,
                        Geodetic (dms 54 44 48.71) (dms (-2) 0 (-14.41)) m0 OSGB36) ]
   where m0 = 0 *~ meter

-- Convert a named point into a test
pointTest :: (Ellipsoid e2) => (String, Geodetic WGS84, Geodetic e2) -> Test
pointTest (name, wgs84, local) =  testCase name $ HU.assertBool "" $ samePlace wgs84 (toWGS84 local)


-- The negation of the sum of a list of offsets is equal to the sum of the negated items.
prop_offset1 :: [GridOffset] -> Bool
prop_offset1 offsets = sameOffset (offsetNegate $ mconcat offsets) (mconcat $ map offsetNegate offsets)

-- A polar offset multiplied by a scalar is equal to an offset in the same direction with the length multiplied.
prop_offset2 :: Distance -> Bearing -> Scalar -> Bool
prop_offset2 (Distance d) (Bearing h) (Scalar s) = sameOffset go1 go2
   where 
      go1 = offsetScale s $ polarOffset d h
      go2 = polarOffset (d * s) h

-- | A polar offset has the offset distance and bearing of its arguments.
prop_offset3 :: GridOffset -> Bool
prop_offset3 delta = sameOffset delta0 
                                (polarOffset (offsetDistance delta0) (offsetBearing delta))
   where delta0 = delta {deltaAltitude = 0 *~ meter}

-- | Given a grid point and an offset, applying the offset to the point gives a new point which
-- is offset from the first point by the argument offset.
prop_grid1 :: GridPoint (GridTM LocalEllipsoid) -> GridOffset -> Bool
prop_grid1 p d = sameOffset d $ p `gridOffset` (applyOffset d p)


-- | Converting a UK grid reference to a GridPoint and back is a null operation.
prop_ukGrid1 :: GridRef -> Bool
prop_ukGrid1 (GridRef str) = 
   str ==
   (fromJust $ toUkGridReference ((length str P.- 2) `div` 2) $ fst $ fromJust $ fromUkGridReference str)

-- | UK Grid Reference points. The oracle for these points was the 
-- UK Grid Reference Finder (gridreferencefinder.com), retrieved on 26 Jan 2013.
ukSampleGrid :: [(String, GridPoint UkNationalGrid, Geodetic WGS84, String)]
ukSampleGrid = map convert [
 -- Grid Reference, X,      Y,      Latitude,  Longitude,    Description
   ("SW3425625070", 134256, 025070, 50.066230, -5.7148278,   "Lands End"),
   ("TR3302139945", 633021, 139945, 51.111396,  1.3277159,   "Dover Harbour"),
   ("TQ3001980417", 530019, 180417, 51.507736, -0.12793230,  "Nelsons Column"),
   ("TA2542370644", 525423, 470644, 54.116376, -0.082668990, "Flamborough Lighthouse"),
   ("NK1354745166", 413547, 845166, 57.496512, -1.7756310,   "Peterhead harbour"),
   ("ND3804872787", 338048, 972787, 58.638518, -3.0688688,   "John O Groats"),
   ("SC3915875189", 239158, 475189, 54.147275, -4.4641148,   "Douglas Harbour"),
   ("ST1922474591", 319224, 174591, 51.464505, -3.1641741,   "Torchwood HQ"),
   ("SK3520736502", 435207, 336502, 52.924784, -1.4777486,   "Derby Cathedral")]
   where
      convert (grid, x, y, lat, long, desc) = 
         (grid, GridPoint (x *~ meter) (y *~ meter) (0 *~ meter) UkNationalGrid,
          Geodetic (lat *~ degree) (long *~ degree) (0 *~ meter) WGS84, desc)

type GridPointTest = (String, GridPoint UkNationalGrid, Geodetic WGS84, String) -> Test

-- | Check that grid reference to grid point works for sample points.
ukGridTest2 :: GridPointTest
ukGridTest2 (gridRef, gp, _, name) = testCase name $ HU.assertBool "" 
   $ (fst $ fromJust $ fromUkGridReference gridRef) == gp

-- | Check that grid point to grid reference works for sample points.
ukGridTest3 :: GridPointTest
ukGridTest3 (gridRef, gp, _, name) = testCase name $ HU.assertBool "" 
   $ toUkGridReference 5 gp == Just gridRef

-- | Check that grid point to WGS84 works close enough for sample points. 
ukGridTest4 :: GridPointTest
ukGridTest4 (_, gp, geo, name) = testCase name $ HU.assertBool ""
   $ closeEnough geo $ toWGS84 $ fromGrid gp
   
-- | Check that WGS84 to grid point works close enough for sample points.
ukGridTest5 :: GridPointTest
ukGridTest5 (_, gp, geo, name) = testCase name $ HU.assertBool ""
   $ offsetDistance (gridOffset gp $ toGrid UkNationalGrid $ toLocal OSGB36 geo) < 10 *~ meter

-- Some sample values for ad-hoc tests
e :: LocalEllipsoid
e = LocalE {nameLocal = "Local_FQN", majorRadiusLocal = 6378349.238479761 *~ meter, flatRLocal = 297.1206701694604 *~ one, helmertLocal = Helmert {cX = 49.0 *~ meter, cY = 86.0 *~ meter, cZ = 50.0 *~ meter, helmertScale = 4.0 *~ one, rX = 2e-4 *~ one, rY = 2.827464643683094e-4 *~ one, rZ = (-3e-4) *~ one}}

pA = Geodetic (1 *~ degree) (51 *~ degree) (0 *~ meter) e


-- Worked example for UK Geodetic to GridPoint, taken from "A Guide to Coordinate Systems in Great Britain"

ukTest :: Geodetic OSGB36
ukTest = Geodetic (dms 52 39 27.2531) (dms 1 43 4.5177) (0 *~ meter) OSGB36

{- 
   v = 6.3885023333E+06
   rho = 6.3727564399E+06
   eta2 = 2.4708136169E-03
   m = 4.0668829596E+05
   I = 3.0668829596E+05
   II = 1.5404079092E+06
   III = 1.5606875424E+05
   IIIa = -2.0671123011E+04
   IV = 3.8751205749E+06
   V = -1.7000078208E+05
   VI = -1.0134470432E+05
   E = 651409.903 m
   N = 313177.270 m
-}

