geodetics
=========

Haskell library of data types and calculations for positions on planet Earth

This library provides "geodetic" positions. That is, latitude, longitude and altitude on a 
specified Terrestrial Reference Frame (TRF). The basic TRF is the WGS84, which is the
one used by GPS and Google Earth. Others can be added by describing the underlying ellipsoid
and the difference in angle and centre with WGS84, and a position in one TRF can be
transformed into another. Given two points in the same TRF you can find the shortest distance 
between them and the bearing from one to the other.

Once you have a geodetic position defined you can project it onto a flat plane, or Grid.
At present Transverse Mercator and Oblique Stereographic grids are provided. More can be
added by defining new instances of the Grid typeclass: see "AddingProjections.txt" for
detais.

The Paths module defines a path as a parametric function of distance that returns a
position and a bearing. Given two paths you can find their intersection using a fast
iterative algorithm.

Release Notes
-------------

Version 0.0.2: Tided up cabal file and removed spurious dependency on Parsec.

Version 0.0.3: Updated for Haskell Platform 2014.2.0.0 and GHC 7.8.3. Fixed
   some minor documentation issues.

Version 0.0.4: Updated for Dimensional 1.0.
