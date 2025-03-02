{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE NumericUnderscores    #-}
{-# LANGUAGE FlexibleContexts #-}

{- | Universal Transverse Mercator (UTM)

The UTM grid system covers the whole world between 84°N and 80°S. It divides the world into 
grid zones of 6° longitude by 8° latitude. Each zone has a 2 digit number for longitude and
a letter for latitude. This regular system has two exceptions:

* North of Norway the zones 32X, 34X and 36X are not used, with 31X, 33X, 35X and 37X being
  wider instead.

* Zone 32V is widened to cover the south-western end of Norway.

There are two notations for writing UTM grid positions:

* The UTM standard: Zone number, N or S for hemisphere, and then northings and eastings
  relative to the equator.

* The Military Grid Reference System (MGRS): Zone number, latitude band letter, a
  2 letter code for the 100km square within the zone, and then northings and eastings within
  that square.

In this library each UTM longitude zone has two grids, one for the northern hemisphere and
one for the south.

For more details see

* https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system.

* THE UNIVERSAL GRIDS: Universal Transverse Mercator (UTM) and Universal Polar Stereographic (UPS).
  DMA Technical Manual. AD-A226497. https://apps.dtic.mil/sti/tr/pdf/ADA266497.pdf
-}
module Geodetics.UTM (
  UtmHemisphere (..),
  UtmZoneNumber,
  utmZoneNumber,
  UtmZone (utmHemisphere, utmZoneNum, utmProjection),
  utmZone,
  mkUtmZone,
  mkUtmZoneUnsafe,
  fromUtmGridReference,
  toUtmGridReference,
  mgrsBandLetterToLatitude,
  mgrsLatitudeToBandLetter,
  fromMgrsGridReference,
  toMgrsGridReference
) where

import Control.Monad (mplus, guard, void, when, unless)
import Data.Array
import Data.Char
import Data.List
import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid
import Geodetics.TransverseMercator
import Text.Parsec
import Text.Parsec.Error
import Text.Printf
import Text.Read



-- | In UTM the northern and southern hemispheres have different false origins.
data UtmHemisphere = UtmNorth | UtmSouth deriving Eq

instance Show UtmHemisphere where
  show UtmNorth = "N"
  show UtmSouth = "S"


-- | A UTM Zone number. Must be between 1 and 60.
type UtmZoneNumber = Int


-- | A UTM Zone, representing a band of typically 6 degrees of latitude between the equator and one of
-- the poles. The projection *must* match the hemisphere and zone.
data UtmZone = UtmZone {
  utmHemisphere :: UtmHemisphere,
  utmZoneNum :: UtmZoneNumber,
  utmProjection :: GridTM WGS84
} deriving (Show)

instance Eq UtmZone where
  z1 == z2  = utmHemisphere z1 == utmHemisphere z2 && utmZoneNum z1 == utmZoneNum z2

instance GridClass UtmZone WGS84 where
  fromGrid p = fromGrid $ unsafeGridCoerce (utmProjection $ gridBasis p) p
  toGrid grid = unsafeGridCoerce grid . toGrid (utmProjection grid)
  gridEllipsoid _ = WGS84


-- Internal data type representing a "rectangle" of latitude/longitude with an exceptional zone number.
data UtmException = UtmE {
  uteSW :: (Int, Int),  -- South west corner in integer degrees (lat, long), inclusive.
  uteNE :: (Int, Int),  -- North east corner in integer degrees (lat, long), exclusive.
  uteActual :: UtmZoneNumber
} deriving Show


-- | Determine if the integer latitude and longitude are within the exception area.
inException :: Int -> Int -> UtmException -> Bool
inException lat long e =
    inR lat  (fst $ uteSW e) (fst $ uteNE e) &&
    inR long (snd $ uteSW e) (snd $ uteNE e)
  where
    inR v v1 v2 = v1 <= v && v < v2


-- The UTM zone that encloses a given geodetic position. For most of the world this is based on
-- @longitude/6@, but there are exceptions around Norway and Svalbard.
utmZoneNumber :: Geodetic a -> Maybe UtmZoneNumber
utmZoneNumber geo = do
    guard $ lat1 >= (-80) && lat1 < 84
    return $ maybe zone1 uteActual exception
  where
    lat1 = floor $ latitude geo / degree
    long1 = floor $ longitude geo / degree
    zone1 = (long1 `div` 6 + 30) `mod` 60 + 1
    exception = find (inException lat1 long1)
      [
        UtmE (56,03) (64,12) 32,  -- Southwestern end of Norway around Bergen.
        UtmE (72,00) (84,09) 31,
        UtmE (72,09) (84,21) 33,  -- Svalbard.
        UtmE (72,21) (84,33) 35,
        UtmE (72,33) (84,42) 37
      ]


-- | The UTM Zone for the given location, if it exists.
utmZone :: Geodetic a -> Maybe UtmZone
utmZone geo = do
  let hemi = if latitude geo >= 0 then UtmNorth else UtmSouth
  zn <- utmZoneNumber geo
  mkUtmZone hemi zn


-- | Construct a UTM Zone value. Returns @Nothing@ if the zone number is out of range.
mkUtmZone :: UtmHemisphere -> UtmZoneNumber -> Maybe UtmZone
mkUtmZone h n = do
    guard $ n >= 1 && n <= 60
    return $ mkUtmZoneUnsafe h n


-- | Construct a UTM Zone value without checking whether the zone number is valid.
mkUtmZoneUnsafe :: UtmHemisphere -> UtmZoneNumber -> UtmZone
mkUtmZoneUnsafe h n = UtmZone h n $ mkGridTM trueO falseO scale
  where
    trueO = Geodetic 0 (degree * fromIntegral (n * 6 - 183)) 0 WGS84
    falseO = case h of
      UtmNorth -> GridOffset (-500_000) 0 0
      UtmSouth -> GridOffset (-500_000) (-10_000_000) 0
    scale = 0.999_6


-- | Units for UTM grid coordinates.
data UtmGridUnit = UtmMeters | UtmKilometers deriving (Eq, Show)


-- | Convert a grid reference to a position, if the reference is valid.
--
-- The northings and eastings cannot contain more than 20 digits each,
-- including an optional decimal point. Negative values are not permitted.
--
-- Northings and eastings can each be followed by an optional unit. The unit
-- must be either \"m\" or \"km\". The units for both
-- must be the same because otherwise its probably an error. The default is meters.
--
-- Northings may be followed by an \"N\" and Eastings may be followed by an \"E\".
--
-- If the argument cannot be parsed then one or more error messages are returned.
fromUtmGridReference :: String -> Either [String] (GridPoint UtmZone)
fromUtmGridReference str = case parse gridP str str of
    Left err -> Left $ lines $ showErrorMessages
      "or" "unknown parse error" "expecting" "unexpected" "end of input"
      (errorMessages err)
    Right r -> Right r
  where
    gridP = do
      spaces1
      zone <- readZone <?> "Zone number"
      hemi <- readHemi <?> "Hemisphere (N or S)"
      spaces1
      (eastings1, eastUnit) <- readDistance
      spaces
      optional (oneOf "Ee" <?> "E")
      spaces1
      (northings1, northUnit) <- readDistance
      unless (eastUnit == northUnit) $ fail "Northings and Eastings units don't match."
      spaces1
      optional (oneOf "Nn" <?> "N")
      spaces1
      eof
      return $ GridPoint eastings1 northings1 0 $ mkUtmZoneUnsafe hemi zone
    readZone :: Parsec String () UtmZoneNumber
    readZone = do
      ds <- many1 digit
      case readMaybe ds of
        Nothing -> fail "Zone number not found."
        Just n ->
          if n < 1 || n > 60
            then fail $ "Zone number " <> show n <> " out of range."
            else return n
    readHemi :: Parsec String () UtmHemisphere
    readHemi = do
      h <- oneOf "NSns"
      case toUpper h of
        'N' -> return UtmNorth
        'S' -> return UtmSouth
        _ -> fail $ "Invalid hemisphere: " <> (h : ". Must be N or S.")
    readDistance :: Parsec String () (Double, UtmGridUnit)  -- (Distance, unit)
    readDistance = do
      digits <- many1 (digit <|> char '.' <?> "number")
      spaces1
      when (length digits > 20) $ fail "Too many digits."
      (multiplier, unit) <- do
        unit <- option "m" (string1' "m" <|> string1' "km" <?> "units (m or km)")
        if unit == "km" then return (1000, UtmKilometers) else return (1, UtmMeters)
      case readMaybe digits of
        Just d -> return (d * multiplier, unit)
        Nothing -> fail $ "Cannot read number: " <> digits
    string1' target = try $ do  -- Case-insensitive version of string'
      cs <- count (length target) anyToken
      if map toLower target == map toLower cs then return cs else unexpected cs
    spaces1 = void $ many (char ' ' <?> "space")  -- Other white space not permitted.


-- | Convert a grid point to a UTM grid reference.
-- The northings and eastings are rounded down to the resolution, so the result is the south-west
-- corner of the grid square enclosing the grid point.
toUtmGridReference ::
  Maybe UtmGridUnit  -- ^ Include explicit units in the output. @Nothing@ means meters without units.
  -> Bool -- ^ Include \"E\" and \"N\" in the output.
  -> Int  -- ^ Digits of resolution. 0 = 1m resolution, 1 = 10m, 2 = 100m etc. (-2) = 1cm.
  -> GridPoint UtmZone
  -> String
toUtmGridReference unit letters res gp =
    zoneStr <> " " <>
    dist (eastings gp)  <> (if letters then "E " else " ") <>
    dist (northings gp) <> (if letters then "N" else "")
  where
    res1 :: Double
    res1 = 10 ** fromIntegral res   -- Resolution in meters.
    floorRes :: Double -> Double
    floorRes d = res1 * fromIntegral (floor (d/res1) :: Integer)
    b = gridBasis gp
    zoneStr = printf "%02d" (utmZoneNum b) <> show (utmHemisphere b)
    dist d = case unit of
      Nothing            -> printf "%.*f" (-res) $ floorRes d
      Just UtmMeters     -> printf "%.*fm" (-res) $ floorRes d
      Just UtmKilometers -> printf "%.*fkm" (3-res) $ floorRes d / 1000



-- | The MGRS latitude band code letters, excluding A and B used for Antarctica (south of -80 degrees)
-- and Y and Z used for the Arctic (north of 84 degrees).
mgrsBandLetters :: [Char]
mgrsBandLetters = "CDEFGHJKLMNPQRSTUVWX"


-- | Find the southern boundary of a latitude band letter.
mgrsBandLetterToLatitude :: Char -> Maybe Double
mgrsBandLetterToLatitude band = do
    n1 <- ix band
    return $ degree * fromIntegral (-80 + n1 * 8)
  where
    indexMap :: Array Char (Maybe Int)
    indexMap = accumArray mplus Nothing ('A', 'Z') [(c1, Just n) | (c1, n) <- zip mgrsBandLetters [0..]]
    ix c1 = if inRange (bounds indexMap) c1 then indexMap ! c1 else Nothing


-- | Find the band letter for a latitude, if it is in the range (-80, 84) degrees.
-- (Argument in radians)
mgrsLatitudeToBandLetter :: Double -> Maybe Char
mgrsLatitudeToBandLetter lat = do
    guard $ -80 <= dlat && dlat <= 84
    return $ indexMap ! latIdx
  where
    dlat = lat / degree
    ilat :: Int
    ilat = floor dlat
    latIdx = min 19 $ (ilat + 80) `div` 8  -- Band 19 (X) extends an extra 4 degrees.
    indexMap = listArray (0,19) mgrsBandLetters



-- | Letters A-Z except for I and O.
mgrsEastingsLetters :: [Char]
mgrsEastingsLetters = "ABCDEFGHJKLMNPQRSTUVWXYZ"

-- | If zone number is in range and the letter is one of the valid Eastings letters for that zone
-- then return the UTM easting in meters.
--
-- Zone 1 starts with \'A\'. Each zone is 8 characters wide. Hence the letters repeat every 3 zones.
mgrsLetterToEasting :: UtmZoneNumber -> Char -> Maybe Double
mgrsLetterToEasting zn c = do
    guard $ 1 <= zn && zn <= 60
    n1 <- ix c
    let n2 = n1 - base
    guard $ 0 <= n2 && n2 <= 7
    return $ fromIntegral (n2+1) * 100_000
  where
    indexMap = accumArray mplus Nothing ('A', 'Z') [(c1, Just n) | (c1, n) <- zip mgrsEastingsLetters [0..]]
    base = ((zn-1) `mod` 3) * 8
    ix c1 = if inRange (bounds indexMap) c1 then indexMap ! c1 else Nothing


-- | If the zone number is in range and the eastings are between 100,000 and 900,000 then
-- return the Eastings letter.
mgrsEastingToLetter :: UtmZoneNumber -> Double -> Maybe Char
mgrsEastingToLetter zn east = do
    guard $ 1 <= zn && zn <= 60
    guard $ 100 * kilometer <= east && east < 900 * kilometer
    return $ indexMap ! ix
  where
    indexMap = listArray (0,23) mgrsEastingsLetters
    base = ((zn-1) `mod` 3) * 8
    square = max 0 $ min 7 $ floor $ (east - 100 * kilometer)/(100 * kilometer)  -- Clamped in range (0,7).
    ix = base + square  -- Must be in range (0,23)


-- | Letters A-V except for I and O.
mgrsNorthingsLetters :: [Char]
mgrsNorthingsLetters = "ABCDEFGHJKLMNPQRSTUV"


-- | MGRS Northings letters have rather complex relationship to the latitude bands. The 20 letters
-- repeat every 2,000km going north and south from the equator, so the latitude band is needed
-- to disambiguate which repetition of the Northings letter is meant.

-- Unfortunately this repetition of the letters does not neatly coincide with
-- the latitude band boundaries, which are based on degrees of latitude.

-- The base letter just north of the equator is A in odd-numbered zones and F in even numbered zones.
--
-- This uses the latitude band to estimate the range of northings that would be valid there,
-- and hence determine which possible grid band is meant by the northings letter.
-- The algorithm used is approximate and deliberately very forgiving: it will accept some grid squares
-- which are north or south of the band given.
mgrsLetterToNorthings ::
  UtmZoneNumber
  -> Char  -- ^ Latitude band letter (@C@ - @X@ excluding @I@ and @O@).
  -> Char  -- ^ MGRS Northings letter (@A@ - @V@ excluding @I@ and @O@).
  -> Maybe Double
mgrsLetterToNorthings zone bandC northingsC = do
  guard $ 1 <= zone && zone <= 60
  band <- (/degree) <$> mgrsBandLetterToLatitude bandC
  northings0 <- (baseNorthingsOffset +) <$> ix northingsC
  let bandDist = band * metersPerDegree -- Approx dist from equator to southern edge of band.
      bandGridLower = floor $ bandDist / 100_000 - 2  -- Lower limit of band in 100km units
      bandGridUpper = ceiling $ if band > 71 * degree  -- Upper limit of band in 100km units.
        then (bandDist + 12 * metersPerDegree) / 100_000 + 1  -- Band X.
        else (bandDist + 8 * metersPerDegree) / 100_000 + 2  -- Other bands.
      rep = (bandGridLower - northings0 - 1) `div` 20  -- Lower limit in 2,000,000km units.
      grid = (rep+1)*20 + northings0
  guard $ grid >= bandGridLower
  guard $ grid <= bandGridUpper
  return $ fromIntegral $ grid * 100_000
  where
    metersPerDegree = 10_002_000 / 90  -- Equator to north pole.
    baseNorthingsOffset :: Int
    baseNorthingsOffset = if odd zone then 0 else -5
    indexMap :: Array Char (Maybe Int)
    indexMap = accumArray mplus Nothing ('A', 'Z') [(c, Just n) | (c, n) <- zip mgrsEastingsLetters [0..]]
    ix c1 = if inRange (bounds indexMap) c1 then indexMap ! c1 else Nothing


-- | Find the northings letter of the 100km square containing the given Northings.
--
-- The input is not range checked. It just assumes that the northings letters repeat forever.
mgrsNorthingToLetter :: UtmZoneNumber -> Double -> Char
mgrsNorthingToLetter zone northings1 =
  letters ! ((gridNum + baseNorthingsOffset) `mod` 20)
  where
    gridNum :: Int
    gridNum = floor $ northings1 / (100 * kilometer)
    baseNorthingsOffset = if odd zone then 0 else 5
    letters = listArray (0,19) mgrsNorthingsLetters


-- | Convert an MGRS grid reference to a UTM @GridPoint@, if the reference is valid.
-- E.g. \"30U XC 99304 10208\" is the grid reference for Nelson's Column in London.
-- 
-- If the input contains spaces then these are used to delimit the fields. Any or all spaces
-- may be omitted. Multiple spaces are treated as a single space.
--
-- If the reference is valid this returns the position of the south-west corner of the
-- nominated grid square and an offset to its centre. Altitude is set to zero.
fromMgrsGridReference :: String -> Either [String] (GridPoint UtmZone, GridOffset)
fromMgrsGridReference str = case parse mgrsP str str of
    Left err -> Left $ filter (not . null) $ lines $ showErrorMessages
      "or" "unknown parse error" "expecting" "unexpected" "end of input"
      (errorMessages err)
    Right r -> Right r
  where
    mgrsP = do
      zoneNum <- read <$> many1 digit  -- Safe because we can only read digits here.
      spaces
      band <- upper <?> "latitude band letter"
      let (hemi, falseNorthing) = if band >= 'N'
            then (UtmNorth, 0)
            else (UtmSouth, -10_000_000)
      zone <- maybe (fail "Invalid zone") return $ mkUtmZone hemi zoneNum
      spaces
      squareEast <- upper <?> "eastings letter"
      eastingBase <- maybe (fail "Invalid eastings letter") return $
        mgrsLetterToEasting zoneNum squareEast
      spaces
      squareNorth <- upper <?> "northings letter"
      northingBase <- maybe (fail "Invalid northings letter") return $
        mgrsLetterToNorthings zoneNum band squareNorth
      spaces
      (eastingChars, northingChars) <- try spaced <|> try unspaced <|> noDigits
      when (length eastingChars /= length northingChars) $
        fail "Northings and Eastings must be the same length."
      if northingChars == "" && eastingChars == ""
        then -- No digits, just return the outer 100km grid square
          return (GridPoint eastingBase (northingBase - falseNorthing) 0 zone,
                  GridOffset (100 * kilometer / 2) (100 * kilometer / 2) 0)
        else do
          (northing, offset) <- maybe (fail "Invalid northing digits") return $
            fromGridDigits (100 * kilometer) northingChars
          (easting, _) <- maybe (fail "Invalid easting digits") return $
            fromGridDigits (100 * kilometer) eastingChars
          return (GridPoint (eastingBase + easting) (northingBase + northing - falseNorthing) 0 zone,
                  GridOffset (offset/2) (offset/2) 0)
    spaced = do
      e <- many1 digit
      skipMany1 space  -- A space is mandatory here.
      n <- many1 digit
      return (e,n)
    unspaced = do
      digits <- many1 digit
      let c = length digits
      when (odd c) $ fail "Northings and Eastings must be the same length."
      return (splitAt (c `div` 2) digits)
    noDigits = do
      eof
      return ("", "")

-- | Convert UTM @GridPoint@ to an MGRS grid reference.
toMgrsGridReference ::
  Bool  -- ^ Include spaces in the output. The standard says no spaces, but they make 
        -- the output easier to read.
  -> Int  -- ^ Number of digits of precision in the easting and northing. Must be 0-5.
  -> GridPoint UtmZone
  -> Maybe String
toMgrsGridReference withSpaces precision gp = do
  guard $ precision >= 0 && precision <= 5
  band <- mgrsLatitudeToBandLetter $ latitude $ fromGrid gp
  let
    zoneNum = utmZoneNum $ gridBasis gp
    northLetter = mgrsNorthingToLetter zoneNum $ northings gp
  eastLetter <- mgrsEastingToLetter zoneNum $ eastings gp
  (_, northDigits) <- toGridDigits (100 * kilometer) precision $ northings gp
  (_, eastDigits) <- toGridDigits (100 * kilometer) precision $ eastings gp
  let part1 = printf "%02d" zoneNum <> [band]
      part2 = [eastLetter, northLetter]
  return $ if withSpaces
    then part1 <> " " <> part2 <> " " <> eastDigits <> " " <> northDigits
    else part1 <> part2 <> eastDigits <> northDigits

