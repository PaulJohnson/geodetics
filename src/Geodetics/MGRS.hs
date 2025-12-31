{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE NumericUnderscores #-}

{- | The Military Grid Reference System (MGRS)

In MGRS there are two syntaxes for grid references:

  1. Between 80 South and 84 North a grid reference has a zone number, latitude band letter, a
    2 letter code for the 100km square within the zone, and then northings and eastings within
    that square.

  2. In the polar regions a grid reference has a latitude band letter (A or B for South, Y or Z 
    for North), a 2 letter code for the 100km square within the polar region, and then northings
    and eastings within that square. There is no zone number in the polar regions.
-}
module Geodetics.MGRS (
  -- * MGRS Grid Basis
  MgrsGrid (..),
  mgrsGrid,
  utmToMgrsPoint,
  upsToMgrsPoint,
  fromMgrsPoint,
  -- * Textual Representation
  mgrsBandLetterToLatitude,
  mgrsLatitudeToBandLetter,
  fromMgrsGridReference,
  toMgrsGridReference
) where

import Control.Monad
import Data.Array
import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid
import Geodetics.PolarStereographic
import Geodetics.UTM
import Text.Parsec
import Text.Parsec.Error
import Text.Printf


-- | MGRS grid references can be anywhere on Earth. Hence the position can be either on a UTM or
-- a UPS grid.
data MgrsGrid =
  MgrsUtm UtmZone
  | MgrsUps UpsGrid
  deriving (Eq, Show)

instance GridClass MgrsGrid WGS84 where
  fromGrid p = case gridBasis p of
    MgrsUtm zone -> fromGrid $ unsafeGridCoerce zone p
    MgrsUps grid -> fromGrid $ unsafeGridCoerce grid p

  toGrid (MgrsUtm zone) = unsafeGridCoerce (MgrsUtm zone) . toGrid zone
  toGrid (MgrsUps grid) = unsafeGridCoerce (MgrsUps grid) . toGrid grid

  gridEllipsoid _ = WGS84


-- | Find the most appropriate grid for the given geodetic position
mgrsGrid :: Geodetic WGS84 -> MgrsGrid
mgrsGrid geo =
  case utmZone geo of
    Just zone -> MgrsUtm zone
    Nothing -> MgrsUps $ if latitude geo > 0 then upsNorth else upsSouth


-- | Generalise from UTM to MGRS
utmToMgrsPoint :: GridPoint UtmZone -> GridPoint MgrsGrid
utmToMgrsPoint gp = unsafeGridCoerce (MgrsUtm $ gridBasis gp) gp


-- | Generalise from UPS to MGRS
upsToMgrsPoint :: GridPoint UpsGrid -> GridPoint MgrsGrid
upsToMgrsPoint gp = unsafeGridCoerce (MgrsUps $ gridBasis gp) gp


-- | Convert an MGRS grid point to either UTM or UPS, depending on its basis.
fromMgrsPoint :: GridPoint MgrsGrid -> Either (GridPoint UtmZone) (GridPoint UpsGrid)
fromMgrsPoint gp = case gridBasis gp of
  MgrsUtm zone -> Left $ unsafeGridCoerce zone gp
  MgrsUps grid -> Right $ unsafeGridCoerce grid gp


-- | Convert a list of letters into an inverse lookup function.
letterTable :: [Char] -> Char -> Maybe Int
letterTable cs = \c -> if inRange (bounds arr) c then arr ! c else Nothing
  where arr = accumArray mplus Nothing ('A', 'Z') [(c1, Just n) | (c1, n) <- zip cs [0..]]
  -- Second argument in a lambda so that partial application returns a function with a CAF
  -- rather than computing the array for every call.


-- | The MGRS latitude band code letters, excluding A and B used for Antarctica (south of -80 degrees)
-- and Y and Z used for the Arctic (north of 84 degrees).
mgrsBandLetters :: [Char]
mgrsBandLetters = "CDEFGHJKLMNPQRSTUVWX"

mgrsBandIx :: Char -> Maybe Int
mgrsBandIx = letterTable mgrsBandLetters

-- | Polar regions have special bands A, B, Y and Z.
mgrsBandLetterToPole :: Char -> Maybe Pole
mgrsBandLetterToPole 'A' = Just SouthPole
mgrsBandLetterToPole 'B' = Just SouthPole
mgrsBandLetterToPole 'Y' = Just NorthPole
mgrsBandLetterToPole 'Z' = Just NorthPole
mgrsBandLetterToPole _ = Nothing

polarBandLetters :: [Char]
polarBandLetters = "ABYZ"

-- | Polar regions use different northings and eastings letters to the rest of the world.
polarEastingLetters :: [Char]
polarEastingLetters = "ABCFGHJKLPQRSTUXYZ"

polarEastingIx :: Char -> Maybe Int
polarEastingIx = letterTable polarEastingLetters

polarNorthingLetters :: [Char]
polarNorthingLetters = "ABCDEFGHJKLMNPQRSTUVWXYZ"

polarNorthingIx :: Char -> Maybe Int
polarNorthingIx = letterTable polarNorthingLetters


-- | Convert the three letters of an MGRS grid reference to the south-west corner of a 100km square.
--
-- If any of the letters are invalid for the positions they are in then an error message is returned.
mgrsLettersToPolar ::
  Char   -- ^ Latitude band: A, B, Y or Z.
  -> Char  -- ^ Eastings letter.
  -> Char  -- ^ Northings letter.
  -> Parsec s u (Pole, Double, Double)
mgrsLettersToPolar bandC eastC northC =
  case (pole, baseEasting, baseNorthing) of
    (Just pole1, Just easting, Just northing) -> pure (pole1, easting, northing + polarNorthingLetterOrigin pole1)
    _ -> fail "Invalid polar grid letters"
  where
    pole = mgrsBandLetterToPole bandC
    baseEasting = (bandEasting +) . (100_000 *) . fromIntegral <$> polarEastingIx eastC
    baseNorthing = (100_000 *) . fromIntegral <$> polarNorthingIx northC
    bandEasting = if bandC `elem` "AY" then 200*kilometer else 2000*kilometer
      -- A and Y denote the "eastern" halves of the south and north polar regions respectively.


-- | Find the southern boundary of a latitude band letter (excluding poles).
mgrsBandLetterToLatitude :: Char -> Maybe Double
mgrsBandLetterToLatitude band = do
    n1 <- mgrsBandIx band
    return $ degree * fromIntegral (-80 + n1 * 8)


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

mgrsEastingIx :: Char -> Maybe Int
mgrsEastingIx = letterTable mgrsEastingsLetters

-- | If zone number is in range and the letter is one of the valid Eastings letters for that zone
-- then return the UTM easting in meters.
--
-- Zone 1 starts with \'A\'. Each zone is 8 characters wide. Hence the letters repeat every 3 zones.
mgrsLetterToEasting :: UtmZoneNumber -> Char -> Maybe Double
mgrsLetterToEasting zn c = do
    guard $ 1 <= zn && zn <= 60
    n1 <- mgrsEastingIx c
    let n2 = n1 - base
    guard $ 0 <= n2 && n2 <= 7
    return $ fromIntegral (n2+1) * 100_000
  where
    base = ((zn-1) `mod` 3) * 8


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

mgrsNorthingIx :: Char -> Maybe Int
mgrsNorthingIx = letterTable mgrsEastingsLetters


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
  northings0 <- (baseNorthingsOffset +) <$> mgrsNorthingIx northingsC
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


polarEastingsToLetter :: Double -> Char
polarEastingsToLetter eastings1 =
  letters ! ((floor (eastings1 / (100 * kilometer)) - 20) `mod` 18 :: Int)
  where
    letters = listArray (0,17) polarEastingLetters


-- | The southern edge of the @A@ northings letter differs between north and south poles.
polarNorthingLetterOrigin :: Pole -> Double
polarNorthingLetterOrigin NorthPole = 1_300 * kilometer
polarNorthingLetterOrigin SouthPole =   800 * kilometer


polarNorthingsToLetter :: Pole -> Double -> Maybe Char
polarNorthingsToLetter pole northings1 = do
  let i :: Int
      i = floor ((northings1 - polarNorthingLetterOrigin pole) / (100 * kilometer)) `mod` 24
  guard $ i >= 0 && i < 24
  pure $ letters ! i
  where
    letters = listArray (0,23) polarNorthingLetters

-- | Convert an MGRS grid reference to a UTM @GridPoint@, if the reference is valid.
-- E.g. \"30U XC 99304 10208\" is the grid reference for Nelson's Column in London.
-- 
-- If the input contains spaces then these are used to delimit the fields. Multiple
-- spaces are treated as a single space.
--
-- If the reference is valid this returns the position of the south-west corner of the
-- nominated grid square and an offset to its centre. Altitude is set to zero.
fromMgrsGridReference :: String -> Either [String] (GridPoint MgrsGrid, GridOffset)
fromMgrsGridReference str = case parse parseMgrsGridReference "" str of
    Left err -> Left $ filter (not . null) $ lines $ showErrorMessages
      "or" "unknown parse error" "expecting" "unexpected" "end of input"
      (errorMessages err)
    Right r -> Right r

parseMgrsGridReference :: Parsec String u (GridPoint MgrsGrid, GridOffset)
parseMgrsGridReference = mgrsUtmParser <|> mgrsUpsParser
  where
    mgrsUtmParser = do
      zoneNum <- read <$> many1 digit <?> "UTM zone number"
          -- Safe because we can only read digits
      spaces
      band <- oneOf mgrsBandLetters <?> "latitude band letter"
      let (hemi, falseNorthing) = if band >= 'N'
            then (UtmNorth, 0)
            else (UtmSouth, -10_000_000)
      zone <- maybe (fail "Invalid zone") (pure . MgrsUtm) $ mkUtmZone hemi zoneNum
      spaces
      squareEast <- oneOf mgrsEastingsLetters <?> "eastings letter"
      eastingBase <- maybe (fail "Invalid eastings letter") pure $
        mgrsLetterToEasting zoneNum squareEast
      spaces
      squareNorth <- oneOf mgrsNorthingsLetters <?> "northings letter"
      northingBase <- maybe (fail "Invalid northings letter") pure $
        mgrsLetterToNorthings zoneNum band squareNorth
      spaces
      (eastingChars, northingChars) <- try spaced <|> try unspaced <|> noDigits
      when (length eastingChars /= length northingChars) $
        fail "Northings and Eastings must be the same length."
      if northingChars == "" && eastingChars == ""
        then -- No digits, just return the outer 100km grid square
          pure (GridPoint eastingBase (northingBase - falseNorthing) 0 zone,
                  GridOffset (100 * kilometer / 2) (100 * kilometer / 2) 0)
        else do
          (northing, offset) <- maybe (fail "Invalid northing digits") pure $
            fromGridDigits (100 * kilometer) northingChars
          (easting, _) <- maybe (fail "Invalid easting digits") pure $
            fromGridDigits (100 * kilometer) eastingChars
          pure (GridPoint (eastingBase + easting) (northingBase + northing - falseNorthing) 0 zone,
                  GridOffset (offset/2) (offset/2) 0)

    mgrsUpsParser = do
      spaces
      band <- oneOf polarBandLetters <?> "polar band letter"
      spaces
      squareEast <- oneOf polarEastingLetters <?> "eastings letter"
      spaces
      squareNorth <- oneOf polarNorthingLetters <?> "northings letter"
      spaces
      (pole, eastingBase, northingBase) <- mgrsLettersToPolar band squareEast squareNorth
      let grid = case pole of
            NorthPole -> MgrsUps upsNorth
            SouthPole -> MgrsUps upsSouth
      (eastingChars, northingChars) <- try spaced <|> try unspaced <|> noDigits
      when (length eastingChars /= length northingChars) $
        fail "Northings and Eastings must be the same length."
      if northingChars == "" && eastingChars == ""
        then -- No digits, just return the outer 100km grid square
          pure (GridPoint eastingBase northingBase 0 grid,
                  GridOffset (100 * kilometer / 2) (100 * kilometer / 2) 0)
        else do
          (northing, offset) <- maybe (fail "Invalid northing digits") pure $
            fromGridDigits (100 * kilometer) northingChars
          (easting, _) <- maybe (fail "Invalid easting digits") pure $
            fromGridDigits (100 * kilometer) eastingChars
          pure (GridPoint (eastingBase + easting) (northingBase + northing) 0 grid,
                  GridOffset (offset/2) (offset/2) 0)
    spaced = do
      e <- many1 digit
      skipMany1 space  -- A space is mandatory here.
      n <- many1 digit
      pure (e,n)
    unspaced = do
      digits <- many1 digit
      let c = length digits
      when (odd c) $ fail "Northings and Eastings must be the same length."
      pure (splitAt (c `div` 2) digits)
    noDigits = do
      eof
      pure ("", "")

-- | Convert a UTM or UPS @GridPoint@ to an MGRS grid reference.
toMgrsGridReference ::
  Bool  -- ^ Include spaces in the output. The standard says no spaces, but they make 
        -- the output easier to read.
  -> Int  -- ^ Number of digits of precision in the easting and northing. Must be 0-8.
  -> GridPoint MgrsGrid
  -> Maybe String
toMgrsGridReference withSpaces precision gp = do
  guard $ precision >= 0 && precision <= 8
  case gridBasis gp of
    MgrsUtm zone -> do
      band <- mgrsLatitudeToBandLetter $ latitude $ fromGrid gp
      let
        zoneNum = utmZoneNum zone
        northLetter = mgrsNorthingToLetter zoneNum $ northings gp
      eastLetter <- mgrsEastingToLetter zoneNum $ eastings gp
      (_, northDigits) <- toGridDigits (100 * kilometer) precision $ northings gp
      (_, eastDigits) <- toGridDigits (100 * kilometer) precision $ eastings gp
      let part1 = printf "%02d" zoneNum <> [band]
          part2 = [eastLetter, northLetter]
      pure $ if withSpaces
        then part1 <> " " <> part2 <> " " <> eastDigits <> " " <> northDigits
        else part1 <> part2 <> eastDigits <> northDigits
    MgrsUps grid -> do
      let zoneLetter = case trueOrigin grid of
            NorthPole -> if eastings gp < 2_000 * kilometer then 'Y' else 'Z'
            SouthPole -> if eastings gp < 2_000 * kilometer then 'A' else 'B'
          eastLetter = polarEastingsToLetter $ eastings gp
      northLetter <- polarNorthingsToLetter (trueOrigin grid) $ northings gp
      (_, northDigits) <- toGridDigits (100 * kilometer) precision $ northings gp
      (_, eastDigits) <- toGridDigits (100 * kilometer) precision $ eastings gp
      pure $ if withSpaces
        then zoneLetter : ' ' : eastLetter : northLetter : " " <> eastDigits <> " " <> northDigits
        else zoneLetter : eastLetter : northLetter : eastDigits <> northDigits
