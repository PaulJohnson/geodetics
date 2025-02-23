{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE ScopedTypeVariables #-}

module LatLongParser (parserTests) where

import ArbitraryInstances ()
import Geodetics.Geodetic
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck hiding ((===))

-- | The 'Coordinates' that certain test cases must return.
expected0, expected1, expected2, expected3, expected4, expected5 :: Coordinates
expected0 = Coordinates 0.6025447300549828 (-0.8069065539014753)
expected1 = Coordinates 0.6025545620764358 (-0.8069006197820184)
expected2 = Coordinates 0.6024294801467094 (-0.806633002630046)
expected3 = Coordinates 0.6025435083245063 (-0.8069065829902962)
expected4 = Coordinates 0.6025409872933646 (-0.8069044982914674)
expected5 = Coordinates 0.3139836658436815 (-5.259743626357356e-4)

parserTests :: SpecWith ()
parserTests =
  describe "LatLongParser" $ do
    testTree
      "Signed decimal degrees"
      ("34.52327", "-46.23234")
      (commaOrNot (produces expected0))
    testTree
      "Decimal degrees NSEW"
      ("34.52327N", "46.23234W")
      (abOrBa (commaOrNot (produces expected0)))
    testTree
      "Degrees and decimal minutes, with decimal point"
      ("34° 31.43' N", "46° 13.92' W")
      (allVariations expected1)
    testTree
      "Degrees and decimal minutes, without decimal point"
      ("34° 31' N", "46° 13' W")
      (allVariations expected2)
    testTree
      "Degrees, minutes and seconds, with decimal point"
      ("34° 31' 23.52\" N", "46° 13' 56.43\" W")
      (allVariations expected3)
    testTree
      "Degrees, minutes and seconds, without decimal point"
      ("34° 31' 23\" N", "46° 13' 56\" W")
      (allVariations expected4)
    testTree
      "DDDMMSS format, no leading zeros, with decimal point"
      ("343123.52N", "461356.43W")
      (abOrBa (commaOrNot (produces expected3)))
    testTree
      "DDDMMSS format, with leading zeros, with decimal point"
      ("343123.52N", "0461356.43W")
      (abOrBa (commaOrNot (produces expected3)))
    testTree
      "DDDMMSS format, no decimal point"
      ("343123N", "0461356W")
      (abOrBa (commaOrNot (produces expected4)))
    testTree
      "Very small angles in the DDDMMSS format are unambiguous"
      ("175923.78N", "00148.49W")
      (abOrBa (commaOrNot (produces expected5)))
    prop "degrees, minutes, seconds roundtrip" dmsRoundtrip
    prop "signed decimal roundtrip" signedDecimalRoundtrip
    prop "decimal degrees NSEW roundtrip" decimalDegreesNSEWRoundtrip
    prop "DDDMMSS roundtrip" dddmmssRoundtrip


----------------------------------------------------------------------------
-- The Coordinates type

-- | A pair of coordinates: latitude and longitude expressed in radians.
data Coordinates = Coordinates Double Double
  deriving (Show, Eq)

instance Arbitrary Coordinates where
  arbitrary = do
    x :: Geodetic WGS84 <- arbitrary
    return (Coordinates (latitude x) (longitude x))

-- | Convert 'Coordinates' to a 'Geodetic'.
coordinatesToGeodetic :: Coordinates -> Geodetic WGS84
coordinatesToGeodetic (Coordinates lat long) =
  Geodetic
    { latitude = lat,
      longitude = long,
      geoAlt = 0.0,
      ellipsoid = WGS84
    }

-- | Convert a 'Geodetic' to 'Coordinates'.
geodeticToCoordinates :: Geodetic WGS84 -> Coordinates
geodeticToCoordinates x = Coordinates (latitude x) (longitude x)

----------------------------------------------------------------------------
-- The machinery for unit tests

-- | A tree-like structure of test cases.
testTree ::
  -- | Name of the test tree
  String ->
  -- | The input
  input ->
  -- | The test generator
  (input -> SpecWith ()) ->
  -- | The resulting test
  SpecWith ()
testTree testTreeName input generator =
  describe testTreeName (generator input)


-- | Stipulate that parsing the input produces the expected result.
produces ::
  -- | The value that should be produced by the parser
  Coordinates ->
  -- | The input
  String ->
  -- | The resulting collection of tests
  SpecWith (Arg Expectation)
produces expected input =
  it caseName $ Just expected `shouldBe` mcoordinates
  where
    caseName = "Input: \"" ++ input ++ "\" -> " ++ show expected
    mcoordinates = geodeticToCoordinates <$> readGroundPosition WGS84 input


-- | The two 'String's will be examined in various orders.
abOrBa :: ((String, String) -> SpecWith (Arg Expectation)) -> (String, String) -> SpecWith (Arg Expectation)
abOrBa f (a, b) = do
  describe "straight order" (f (a, b))
  describe "reverse order" (f (b, a))


-- | Two two 'String's will be concatenated with a space between them and
-- with a comma and a space between them.
commaOrNot :: (String -> SpecWith (Arg Expectation)) -> (String, String) -> SpecWith (Arg Expectation)
commaOrNot f (a, b) = do
  describe "with comma" (f (a ++ ", " ++ b))
  describe "without comma" (f (a ++ " " ++ b))


-- | Test with and without units both in normal and Unicode version.
variousUnits :: (String -> SpecWith ()) -> String -> SpecWith ()
variousUnits f input = do
  describe "without Unicode" (f input)
  describe "with Unicode" (f unicode)
  describe "without units" (f noUnits)
  where
    unicode = g <$> input
    g '\'' = '′'
    g '\"' = '″'
    g x = x
    noUnits = filter (not . isUnit) input
    isUnit = \case
      '°' -> True
      '\'' -> True
      '\"' -> True
      '′' -> True
      '″' -> True
      _ -> False


-- | Tests with all possible variations.
allVariations :: Coordinates -> (String, String) -> SpecWith (Arg Expectation)
allVariations expected =
  abOrBa (commaOrNot (variousUnits (produces expected)))


----------------------------------------------------------------------------
-- The machinery for roundtrip property tests

-- | A lax version of 'Test.QuickCheck.===' that is forgiving with respect
-- to precision lost due to printing.
infix 4 ===

(===) :: Coordinates -> Coordinates -> Property
Coordinates lata longa === Coordinates latb longb =
  consider "latitude" lata latb
    .&&. consider "longitude" longa longb
  where
    consider str x y =
      counterexample
        (str ++ " " ++ show x ++ " vs " ++ show y)
        (abs (x - y) < epsilon)
    epsilon = 0.0000001

roundtrip :: Coordinates -> String -> Property
roundtrip x rendered =
  counterexample rendered $
    case geodeticToCoordinates <$> readGroundPosition WGS84 rendered of
      Nothing -> property False
      Just y -> y === x

dmsRoundtrip :: Coordinates -> Property
dmsRoundtrip x =
  roundtrip x (showGeodeticLatLong (coordinatesToGeodetic x))

signedDecimalRoundtrip :: Coordinates -> Property
signedDecimalRoundtrip x =
  roundtrip x (showGeodeticSignedDecimal (coordinatesToGeodetic x))

decimalDegreesNSEWRoundtrip :: Coordinates -> Property
decimalDegreesNSEWRoundtrip x =
  roundtrip x (showGeodeticNSEWDecimal (coordinatesToGeodetic x))

dddmmssRoundtrip :: Bool -> Coordinates -> Property
dddmmssRoundtrip useLeadingZeros x =
  roundtrip x (showGeodeticDDDMMSS useLeadingZeros (coordinatesToGeodetic x))
