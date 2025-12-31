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
  parseUtmGridReference,
  toUtmGridReference
) where

import Control.Monad (guard, void, when, unless)
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
fromUtmGridReference str = case parse parseUtmGridReference str str of
    Left err -> Left $ lines $ showErrorMessages
      "or" "unknown parse error" "expecting" "unexpected" "end of input"
      (errorMessages err)
    Right r -> Right r


parseUtmGridReference :: Stream s m Char => ParsecT s u m (GridPoint UtmZone)
parseUtmGridReference = do
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
      pure $ GridPoint eastings1 northings1 0 $ mkUtmZoneUnsafe hemi zone
  where
    readZone :: Stream s m Char => ParsecT s u m UtmZoneNumber
    readZone = do
      ds <- many1 digit
      case readMaybe ds of
        Nothing -> fail "Zone number not found."
        Just n ->
          if n < 1 || n > 60
            then fail $ "Zone number " <> show n <> " out of range."
            else pure n
    readHemi :: Stream s m Char => ParsecT s u m UtmHemisphere
    readHemi = do
      h <- oneOf "NSns"
      case toUpper h of
        'N' -> pure UtmNorth
        'S' -> pure UtmSouth
        _ -> fail $ "Invalid hemisphere: " <> (h : ". Must be N or S.")
    readDistance :: Stream s m Char => ParsecT s u m (Double, Maybe GridUnit)  -- (Distance, unit)
    readDistance = do
      digits <- many1 (digit <|> char '.' <?> "number")
      spaces1
      when (length digits > 20) $ fail "Too many digits."
      (multiplier, unit) <- do
        unit <- optionMaybe (string1' "m" <|> string1' "km" <?> "units (m or km)")
        case unit of
          Just "km" -> pure (1000, Just GridKilometers)
          Just _    -> pure (1, Just GridMeters)
          Nothing   -> pure (1, Nothing)
      case readMaybe digits of
        Just d -> pure (d * multiplier, unit)
        Nothing -> fail $ "Cannot read number: " <> digits
    string1' target = try $ do  -- Case-insensitive version of string'
      cs <- count (length target) anyToken
      if map toLower target == map toLower cs then pure cs else unexpected cs
    spaces1 :: Stream s m Char => ParsecT s u m ()
    spaces1 = void $ many (char ' ' <?> "space")  -- Other white space not permitted.


-- | Convert a grid point to a UTM grid reference.
-- The northings and eastings are rounded down to the resolution, so the result is the south-west
-- corner of the grid square enclosing the grid point.
toUtmGridReference ::
  Maybe GridUnit  -- ^ Include explicit units in the output. @Nothing@ means meters without units.
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
      Just GridMeters     -> printf "%.*fm" (-res) $ floorRes d
      Just GridKilometers -> printf "%.*fkm" (3-res) $ floorRes d / 1000
