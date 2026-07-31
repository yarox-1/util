// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Yannik Schnell <yannik.schnell@email.uni-freiburg.de>

#include "util/Test.h"
#include "util/log/Log.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/geo/RTree.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/tests/GeoTest.h"

using namespace util;
using namespace util::geo;

// _____________________________________________________________________________
void GeoTest::testCRS() {
  {
    // Test with IRI in front of WKT.

    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> MULTIPOINT(0 0, 1 1)", 0), ==,
         WKTType::MULTIPOINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84>MULTIPOINT(0 0, 1 1)", 0), ==,
         WKTType::MULTIPOINT);
    TEST(util::geo::getWKTType(" <http://www.opengis.net/def/crs/OGC/1.3/CRS84> MULTIPOINT(0 0, 1 1)", 0), ==,
         WKTType::MULTIPOINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> MMULTIPOINT(0 0, 1 1)", 0), ==,
         WKTType::MULTIPOINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84>MMULTIPOINT(0 0, 1 1)", 0), ==,
         WKTType::MULTIPOINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84>  MMULTIPOINT (0 0, 1 1)", 0), ==,
         WKTType::MULTIPOINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84>  MMULTIPOINT Z(0 0, 1 1)", 0), ==,
         WKTType::MULTIPOINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84>  MULTIPOINT Z(0 0, 1 1)", 0), ==,
         WKTType::MULTIPOINT);

    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT(1 1)", 0), ==, WKTType::POINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> POIN (1 1)", 0), ==, WKTType::NONE);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT (1 1)", 0), ==, WKTType::POINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> Point (1 1)", 0), ==, WKTType::POINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT Z(1 1)", 0), ==, WKTType::POINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326>   POINT Z (1 1)", 0), ==, WKTType::POINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT ZM(1 1)", 0), ==, WKTType::POINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> Point MZ (1 1)", 0), ==, WKTType::POINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> mPoint MZ (1 1)", 0), ==, WKTType::POINT);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326> oint MZ (1 1)", 0), ==, WKTType::NONE);

    // Invalid IRI.
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326 Point(1 1)", 0), ==, WKTType::NONE);
    TEST(util::geo::getWKTType("<http://www.opengis.net/def/crs/EPSG/0/4326>> Point(1 1)", 0), ==, WKTType::NONE);
  }

  {
    // Test CRS Type
    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> MULTIPOINT(0 0, 1 1)", 0), ==,
         CRSType::CRS84);
    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/OGC/1.3/CRS84>MULTIPOINT(0 0, 1 1)", 0), ==,
         CRSType::CRS84);

    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT(1 1)", 0), ==, CRSType::WGS84);
    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT (1 1)", 0), ==, CRSType::WGS84);

    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/EPSG/0/3857>LINESTRING(0 0, 1 1)", 0), ==, CRSType::WEB_MERCATOR);
    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/EPSG/0/3857>MLINESTRING(0 0, 1 1)", 0), ==, CRSType::WEB_MERCATOR);

    TEST(util::geo::getCRSType("MULTIPOINT(0 0, 1 1)", 0), ==, CRSType::CRS84);
    TEST(util::geo::getCRSType("MMULTIPOINT(0 0, 1 1)", 0), ==, CRSType::CRS84);

    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/EPSG/0/28992> LINESTRING(0 0, 1 1)", 0), ==, CRSType::UNSUPPORTED);
    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/EPSG/0/28992> MLINESTRING(0 0, 1 1)", 0), ==, CRSType::UNSUPPORTED);

    // Invalid IRI.
    TEST(util::geo::getCRSType("<http://www.opengis.net/def/crs/EPSG/0/4326 Point(1 1)", 0), ==, CRSType::UNSUPPORTED);
  }

  {
    // Test transformations between different CRSs.
    util::geo::Point<double> lngLatPoint = Point<double>(4.0, 5.0);
    TEST(util::geo::lngLatToLatLng(lngLatPoint).getX(), ==, 5.0);
    TEST(util::geo::lngLatToLatLng(lngLatPoint).getY(), ==, 4.0);

    util::geo::Point<double> latLngPoint = Point<double>(5.0, 4.0);
    TEST(util::geo::latLngToLngLat(latLngPoint).getX(), ==, 4.0);
    TEST(util::geo::latLngToLngLat(latLngPoint).getY(), ==, 5.0);

    // Test convert function.
    util::geo::Point<double> webMercPoint = Point<double>(445277.96317309426, 557305.2572745768);
    // CRS84 to CRS84
    TEST(util::geo::projectToCRS(lngLatPoint, util::geo::CRS84, util::geo::CRS84).getX(), ==, 4.0);
    TEST(util::geo::projectToCRS(lngLatPoint, util::geo::CRS84, util::geo::CRS84).getY(), ==, 5.0);
    // CRS84 to WGS84
    TEST(util::geo::projectToCRS(lngLatPoint, util::geo::CRS84, util::geo::WGS84).getX(), ==, 5.0);
    TEST(util::geo::projectToCRS(lngLatPoint, util::geo::CRS84, util::geo::WGS84).getY(), ==, 4.0);
    // CRS84 to WEB_MERCATOR
    TEST(util::geo::projectToCRS(lngLatPoint, util::geo::CRS84, util::geo::WEB_MERCATOR).getX(), ==, approx(445277.96317309426));
    TEST(util::geo::projectToCRS(lngLatPoint, util::geo::CRS84, util::geo::WEB_MERCATOR).getY(), ==, approx(557305.2572745768));
    // WGS84 to CRS84
    TEST(util::geo::projectToCRS(latLngPoint, util::geo::WGS84, util::geo::CRS84).getX(), ==, 4.0);
    TEST(util::geo::projectToCRS(latLngPoint, util::geo::WGS84, util::geo::CRS84).getY(), ==, 5.0);
    // WGS84 to WGS84
    TEST(util::geo::projectToCRS(latLngPoint, util::geo::WGS84, util::geo::WGS84).getX(), ==, 5.0);
    TEST(util::geo::projectToCRS(latLngPoint, util::geo::WGS84, util::geo::WGS84).getY(), ==, 4.0);
    // WGS84 to WEB_MERCATOR
    TEST(util::geo::projectToCRS(latLngPoint, util::geo::WGS84, util::geo::WEB_MERCATOR).getX(), ==, approx(445277.96317309426));
    TEST(util::geo::projectToCRS(latLngPoint, util::geo::WGS84, util::geo::WEB_MERCATOR).getY(), ==, approx(557305.2572745768));
    // WEB_MERCATOR to CRS84
    TEST(util::geo::projectToCRS(webMercPoint, util::geo::WEB_MERCATOR, util::geo::CRS84).getX(), ==, approx(4.0));
    TEST(util::geo::projectToCRS(webMercPoint, util::geo::WEB_MERCATOR, util::geo::CRS84).getY(), ==, approx(5.0));
    // WEB_MERCATOR to WGS84
    TEST(util::geo::projectToCRS(webMercPoint, util::geo::WEB_MERCATOR, util::geo::WGS84).getX(), ==, approx(5.0));
    TEST(util::geo::projectToCRS(webMercPoint, util::geo::WEB_MERCATOR, util::geo::WGS84).getY(), ==, approx(4.0));
    // WEB_MERCATOR to WEB_MERCATOR
    TEST(util::geo::projectToCRS(webMercPoint, util::geo::WEB_MERCATOR, util::geo::WEB_MERCATOR).getX(), ==, approx(445277.96317309426));
    TEST(util::geo::projectToCRS(webMercPoint, util::geo::WEB_MERCATOR, util::geo::WEB_MERCATOR).getY(), ==, approx(557305.2572745768));
  }
  {
    // Test already specified convert functions
    util::geo::Point<double> lngLatPoint = Point<double>(4.0, 5.0);
    util::geo::Point<double> latLngPoint = Point<double>(5.0, 4.0);
    util::geo::Point<double> webMercPoint = Point<double>(445277.96317309426, 557305.2572745768);

    // CRS84 to CRS84
    TEST(util::geo::projectToCRS84(lngLatPoint, util::geo::CRS84).getX(), ==, 4.0);
    TEST(util::geo::projectToCRS84(lngLatPoint, util::geo::CRS84).getY(), ==, 5.0);

    // WGS84 to CRS84
    TEST(util::geo::projectToCRS84(latLngPoint, util::geo::WGS84).getX(), ==, 4.0);
    TEST(util::geo::projectToCRS84(latLngPoint, util::geo::WGS84).getY(), ==, 5.0);
    // WEB_MERCATOR to CRS84
    TEST(util::geo::projectToCRS84(webMercPoint, util::geo::WEB_MERCATOR).getX(), ==, approx(4.0));
    TEST(util::geo::projectToCRS84(webMercPoint, util::geo::WEB_MERCATOR).getY(), ==, approx(5.0));

    // CRS84 to WGS84
    TEST(util::geo::projectToWGS84(lngLatPoint, util::geo::CRS84).getX(), ==, 5.0);
    TEST(util::geo::projectToWGS84(lngLatPoint, util::geo::CRS84).getY(), ==, 4.0);
    // WGS84 to WGS84
    TEST(util::geo::projectToWGS84(latLngPoint, util::geo::WGS84).getX(), ==, 5.0);
    TEST(util::geo::projectToWGS84(latLngPoint, util::geo::WGS84).getY(), ==, 4.0);
    //  WEB_MERCATOR to WGS84
    TEST(util::geo::projectToWGS84(webMercPoint, util::geo::WEB_MERCATOR).getX(), ==, approx(5.0));
    TEST(util::geo::projectToWGS84(webMercPoint, util::geo::WEB_MERCATOR).getY(), ==, approx(4.0));

    // CRS84 to WEB_MERCATOR
    TEST(util::geo::projectToWebMerc(lngLatPoint, util::geo::CRS84).getX(), ==, approx(445277.96317309426));
    TEST(util::geo::projectToWebMerc(lngLatPoint, util::geo::CRS84).getY(), ==, approx(557305.2572745768));
    // WGS84 to WEB_MERCATOR
    TEST(util::geo::projectToWebMerc(latLngPoint, util::geo::WGS84).getX(), ==, approx(445277.96317309426));
    TEST(util::geo::projectToWebMerc(latLngPoint, util::geo::WGS84).getY(), ==, approx(557305.2572745768));
    // WEB_MERCATOR to WEB_MERCATOR
    TEST(util::geo::projectToWebMerc(webMercPoint, util::geo::WEB_MERCATOR).getX(), ==, approx(445277.96317309426));
    TEST(util::geo::projectToWebMerc(webMercPoint, util::geo::WEB_MERCATOR).getY(), ==, approx(557305.2572745768));
  }
  {
    // Test PointProjection

    // TO CRS84
    auto projCRS84 = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> POINT(2 3)", util::geo::projectToCRS84<double>);
    TEST(projCRS84.getX(), ==, 2.0);
    TEST(projCRS84.getY(), ==, 3.0);

    auto projWGS84 = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT(3 2)", util::geo::projectToCRS84<double>);
    TEST(projWGS84.getX(), ==, 2.0);
    TEST(projWGS84.getY(), ==, 3.0);

    auto projWebMerc = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> POINT(222638.98158654 334111.17140195)", util::geo::projectToCRS84<double>);
    TEST(projWebMerc.getX(), ==, approx(2.0));
    TEST(projWebMerc.getY(), ==, approx(3.0));

    // to WGS84
    auto projCRS842 = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> POINT(2 3)", util::geo::projectToWGS84<double>);
    TEST(projCRS842.getX(), ==, 3.0);
    TEST(projCRS842.getY(), ==, 2.0);

    auto projWGS842 = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT(3 2)", util::geo::projectToWGS84<double>);
    TEST(projWGS842.getX(), ==, 3.0);
    TEST(projWGS842.getY(), ==, 2.0);

    auto projWebMerc2 = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> POINT(222638.98158654 334111.17140195)", util::geo::projectToWGS84<double>);
    TEST(projWebMerc2.getX(), ==, approx(3.0));
    TEST(projWebMerc2.getY(), ==, approx(2.0));

    // to WebMerc
    auto projCRS843 = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> POINT(2 3)", util::geo::projectToWebMerc<double>);
    TEST(projCRS843.getX(), ==, approx(222638.98158654));
    TEST(projCRS843.getY(), ==, approx(334111.17140195));

    auto projWGS843 = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT(3 2)", util::geo::projectToWebMerc<double>);
    TEST(projWGS843.getX(), ==, approx(222638.98158654));
    TEST(projWGS843.getY(), ==, approx(334111.17140195));

    auto projWebMerc3 = util::geo::pointFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> POINT(222638.98158654 334111.17140195)", util::geo::projectToWebMerc<double>);
    TEST(projWebMerc3.getX(), ==, approx(222638.98158654));
    TEST(projWebMerc3.getY(), ==, approx(334111.17140195));


    // Without projection function.
    auto projCRS844 = util::geo::pointFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> POINT(2 3)");
    TEST(projCRS844.getX(), ==, 2.0);
    TEST(projCRS844.getY(), ==, 3.0);

    auto projWGS844 = util::geo::pointFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> POINT(3 2)");
    TEST(projWGS844.getX(), ==, 2.0);
    TEST(projWGS844.getY(), ==, 3.0);

    auto projWebMerc4 = util::geo::pointFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> POINT(222638.98158654 334111.17140195)");
    TEST(projWebMerc4.getX(), ==, approx(2.0));
    TEST(projWebMerc4.getY(), ==, approx(3.0));
  }
  {
    // Test LineProjection

    // TO CRS84
    auto projCRS84 = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> LINESTRING(7 47, 7 48)", util::geo::projectToCRS84<double>);
    TEST(projCRS84.size(), ==, 2);
    TEST(projCRS84[0].getX(), ==, 7.0);
    TEST(projCRS84[0].getY(), ==, 47.0);
    TEST(projCRS84[1].getX(), ==, 7.0);
    TEST(projCRS84[1].getY(), ==, 48.0);

    auto projWGS84 = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> LINESTRING(47 7, 48 7)", util::geo::projectToCRS84<double>);
    TEST(projWGS84.size(), ==, 2);
    TEST(projWGS84[0].getX(), ==, 7.0);
    TEST(projWGS84[0].getY(), ==, 47.0);
    TEST(projWGS84[1].getX(), ==, 7.0);
    TEST(projWGS84[1].getY(), ==, 48.0);

    auto projWebMerc = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> LINESTRING(779236.435552915 5942074.072431109, 779236.435552915 6106854.834885074)", util::geo::projectToCRS84<double>);
    TEST(projWebMerc.size(), ==, 2);
    TEST(projWebMerc[0].getX(), ==, approx(7.0));
    TEST(projWebMerc[0].getY(), ==, approx(47.0));
    TEST(projWebMerc[1].getX(), ==, approx(7.0));
    TEST(projWebMerc[1].getY(), ==, approx(48.0));

    // to WGS84
    auto projCRS842 = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> LINESTRING(7 47, 7 48)", util::geo::projectToWGS84<double>);
    TEST(projCRS842.size(), ==, 2);
    TEST(projCRS842[0].getX(), ==, approx(47.0));
    TEST(projCRS842[0].getY(), ==, approx(7.0));
    TEST(projCRS842[1].getX(), ==, approx(48.0));
    TEST(projCRS842[1].getY(), ==, approx(7.0));

    auto projWGS842 = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> LINESTRING(47 7, 48 7)", util::geo::projectToWGS84<double>);
    TEST(projWGS842.size(), ==, 2);
    TEST(projWGS842[0].getX(), ==, approx(47.0));
    TEST(projWGS842[0].getY(), ==, approx(7.0));
    TEST(projWGS842[1].getX(), ==, approx(48.0));
    TEST(projWGS842[1].getY(), ==, approx(7.0));

    auto projWebMerc2 = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> LINESTRING(779236.435552915 5942074.072431109, 779236.435552915 6106854.834885074)", util::geo::projectToWGS84<double>);
    TEST(projWebMerc2.size(), ==, 2);
    TEST(projWebMerc2[0].getX(), ==, approx(47.0));
    TEST(projWebMerc2[0].getY(), ==, approx(7.0));
    TEST(projWebMerc2[1].getX(), ==, approx(48.0));
    TEST(projWebMerc2[1].getY(), ==, approx(7.0));

    // to WebMerc
    auto projCRS843 = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> LINESTRING(7 47, 7 48)", util::geo::projectToWebMerc<double>);
    TEST(projCRS843.size(), ==, 2);
    TEST(projCRS843[0].getX(), ==, approx(779236.435552915));
    TEST(projCRS843[0].getY(), ==, approx(5942074.072431109));
    TEST(projCRS843[1].getX(), ==, approx(779236.435552915));
    TEST(projCRS843[1].getY(), ==, approx(6106854.834885074));

    auto projWGS843 = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> LINESTRING(47 7, 48 7)", util::geo::projectToWebMerc<double>);
    TEST(projWGS843.size(), ==, 2);
    TEST(projWGS843[0].getX(), ==, approx(779236.435552915));
    TEST(projWGS843[0].getY(), ==, approx(5942074.072431109));
    TEST(projWGS843[1].getX(), ==, approx(779236.435552915));
    TEST(projWGS843[1].getY(), ==, approx(6106854.834885074));

    auto projWebMerc3 = util::geo::lineFromWKTProj<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> LINESTRING(779236.435552915 5942074.072431109, 779236.435552915 6106854.834885074)", util::geo::projectToWebMerc<double>);
    TEST(projWebMerc3.size(), ==, 2);
    TEST(projWebMerc3[0].getX(), ==, approx(779236.435552915));
    TEST(projWebMerc3[0].getY(), ==, approx(5942074.072431109));
    TEST(projWebMerc3[1].getX(), ==, approx(779236.435552915));
    TEST(projWebMerc3[1].getY(), ==, approx(6106854.834885074));


    // Without projection function.
    auto projCRS844 = util::geo::lineFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> LINESTRING(7 47, 7 48)");
    TEST(projCRS844.size(), ==, 2);
    TEST(projCRS844[0].getX(), ==, 7.0);
    TEST(projCRS844[0].getY(), ==, 47.0);
    TEST(projCRS844[1].getX(), ==, 7.0);
    TEST(projCRS844[1].getY(), ==, 48.0);

    auto projWGS844 = util::geo::lineFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> LINESTRING(47 7, 48 7)");
    TEST(projWGS844.size(), ==, 2);
    TEST(projWGS844[0].getX(), ==, 7.0);
    TEST(projWGS844[0].getY(), ==, 47.0);
    TEST(projWGS844[1].getX(), ==, 7.0);
    TEST(projWGS844[1].getY(), ==, 48.0);

    auto projWebMerc4 = util::geo::lineFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> LINESTRING(779236.435552915 5942074.072431109, 779236.435552915 6106854.834885074)");
    TEST(projWebMerc4.size(), ==, 2);
    TEST(projWebMerc4[0].getX(), ==, approx(7.0));
    TEST(projWebMerc4[0].getY(), ==, approx(47.0));
    TEST(projWebMerc4[1].getX(), ==, approx(7.0));
    TEST(projWebMerc4[1].getY(), ==, approx(48.0));
  }
  {
    // Test PolygonProjection

    // Without projection function.
    auto projCRS84 = util::geo::polygonFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> POLYGON((7 47, 8 47, 8 48, 7 48, 7 47))");
    TEST(getWKT(projCRS84), ==, "POLYGON((7 47,8 47,8 48,7 48,7 47))");

    auto projWGS84 = util::geo::polygonFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> POLYGON((47 7, 47 8, 48 8, 48 7, 47 7))");
    TEST(getWKT(projWGS84), ==, "POLYGON((7 47,8 47,8 48,7 48,7 47))");

    auto projWebMerc = util::geo::polygonFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> POLYGON((779236.435552915 5942074.072431109, 890555.9263461885 5942074.072431109, 890555.9263461885 6106854.834885074, 779236.435552915 6106854.834885074, 779236.435552915 5942074.072431109))");
    TEST(getWKT(projWebMerc), ==, "POLYGON((7 47,8 47,8 48,7 48,7 47))");
  }
  {
    // Test MultiLineProjection

    // Without projection function.
    auto projCRS84 = util::geo::multiLineFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> MULTILINESTRING((7 47, 8 47), (8 48, 7 48))");
    TEST(getWKT(projCRS84), ==, "MULTILINESTRING((7 47,8 47),(8 48,7 48))");

    auto projWGS84 = util::geo::multiLineFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> MULTILINESTRING((47 7, 47 8), (48 8, 48 7))");
    TEST(getWKT(projWGS84), ==, "MULTILINESTRING((7 47,8 47),(8 48,7 48))");

    auto projWebMerc = util::geo::multiLineFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> MULTILINESTRING((779236.435552915 5942074.072431109, 890555.9263461885 5942074.072431109), (890555.9263461885 6106854.834885074, 779236.435552915 6106854.834885074))");
    TEST(getWKT(projWebMerc), ==, "MULTILINESTRING((7 47,8 47),(8 48,7 48))");
  }
  {
    // Test MultiPointProjection

    // Without projection function. MULTIPOINT((1 1), (2 2)) syntax
    auto projCRS84 = util::geo::multiPointFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> MULTIPOINT((7 47), (8 48))");
    TEST(getWKT(projCRS84), ==, "MULTIPOINT(7 47,8 48)");

    auto projWGS84 = util::geo::multiPointFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> MULTIPOINT((47 7), (48 8))");
    TEST(getWKT(projWGS84), ==, "MULTIPOINT(7 47,8 48)");

    auto projWebMerc = util::geo::multiPointFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> MULTIPOINT((779236.435552915 5942074.072431109), (890555.9263461885 6106854.834885074))");
    TEST(getWKT(projWebMerc), ==, "MULTIPOINT(7 47,8 48)");

    // MULTIPOINT(1 1, 2 2) syntax
    auto projCRS842 = util::geo::multiPointFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> MULTIPOINT(7 47, 8 48)");
    TEST(getWKT(projCRS842), ==, "MULTIPOINT(7 47,8 48)");

    auto projWGS842 = util::geo::multiPointFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> MULTIPOINT(47 7, 48 8)");
    TEST(getWKT(projWGS842), ==, "MULTIPOINT(7 47,8 48)");

    auto projWebMerc2 = util::geo::multiPointFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> MULTIPOINT(779236.435552915 5942074.072431109, 890555.9263461885 6106854.834885074)");
    TEST(getWKT(projWebMerc2), ==, "MULTIPOINT(7 47,8 48)");
  }
  {
    // Test MultiPolygonProjection

    // Without projection function.
    auto projCRS84 = util::geo::multiPolygonFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> MULTIPOLYGON(((7 47, 8 47, 8 48, 7 48, 7 47)), ((8 48, 9 48, 9 49, 8 49, 8 48)))");
    TEST(getWKT(projCRS84), ==, "MULTIPOLYGON(((7 47,8 47,8 48,7 48,7 47)),((8 48,9 48,9 49,8 49,8 48)))");

    auto projWGS84 = util::geo::multiPolygonFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> MULTIPOLYGON(((47 7, 47 8, 48 8, 48 7, 47 7)), ((48 8, 48 9, 49 9, 49 8, 48 8)))");
    TEST(getWKT(projWGS84), ==, "MULTIPOLYGON(((7 47,8 47,8 48,7 48,7 47)),((8 48,9 48,9 49,8 49,8 48)))");

    auto projWebMerc = util::geo::multiPolygonFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> MULTIPOLYGON(((779236.435552915 5942074.072431109, 890555.9263461885 5942074.072431109, 890555.9263461885 6106854.834885074, 779236.435552915 6106854.834885074, 779236.435552915 5942074.072431109)), ((890555.9263461885 6106854.834885074, 1001875.4171394621 6106854.834885074, 1001875.4171394621 6274861.394006576, 890555.9263461885 6274861.394006576, 890555.9263461885 6106854.834885074)))");
    TEST(getWKT(projWebMerc), ==, "MULTIPOLYGON(((7 47,8 47,8 48,7 48,7 47)),((8 48,9 48,9 49,8 49,8 48)))");
  }
  {
    // Test CollectionProjection

    // Without projection function.
    auto projCRS84 = util::geo::collectionFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> GEOMETRYCOLLECTION(POINT(7 47), LINESTRING(7 47, 8 48), POLYGON((7 47, 8 47, 8 48, 7 48, 7 47)))");
    TEST(getWKT(projCRS84), ==, "GEOMETRYCOLLECTION(POINT(7 47),LINESTRING(7 47,8 48),POLYGON((7 47,8 47,8 48,7 48,7 47)))");

    auto projWGS84 = util::geo::collectionFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> GEOMETRYCOLLECTION(POINT(47 7), LINESTRING(47 7, 48 8), POLYGON((47 7, 47 8, 48 8, 48 7, 47 7)))");
    TEST(getWKT(projWGS84), ==, "GEOMETRYCOLLECTION(POINT(7 47),LINESTRING(7 47,8 48),POLYGON((7 47,8 47,8 48,7 48,7 47)))");

    auto projWebMerc = util::geo::collectionFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> GEOMETRYCOLLECTION(POINT(779236.435552915 5942074.072431109), LINESTRING(779236.435552915 5942074.072431109, 890555.9263461885 6106854.834885074), POLYGON((779236.435552915 5942074.072431109, 890555.9263461885 5942074.072431109, 890555.9263461885 6106854.834885074, 779236.435552915 6106854.834885074, 779236.435552915 5942074.072431109)))");
    TEST(getWKT(projWebMerc), ==, "GEOMETRYCOLLECTION(POINT(7 47),LINESTRING(7 47,8 48),POLYGON((7 47,8 47,8 48,7 48,7 47)))");

    // Test collection of other geometries.
    auto projCRS842 = util::geo::collectionFromWKT<double>("<http://www.opengis.net/def/crs/OGC/1.3/CRS84> GEOMETRYCOLLECTION(MULTIPOINT(7 47, 8 48), MULTIPOLYGON(((7 47, 8 47, 8 48, 7 48, 7 47)), ((8 48, 9 48, 9 49, 8 49, 8 48))), MULTILINESTRING((7 47, 8 47), (8 48, 7 48)))");
    TEST(getWKT(projCRS842), ==, "GEOMETRYCOLLECTION(MULTIPOINT(7 47,8 48),MULTIPOLYGON(((7 47,8 47,8 48,7 48,7 47)),((8 48,9 48,9 49,8 49,8 48))),MULTILINESTRING((7 47,8 47),(8 48,7 48)))");

    auto projWGS842 = util::geo::collectionFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/4326> GEOMETRYCOLLECTION(MULTIPOINT(47 7, 48 8), MULTIPOLYGON(((47 7, 47 8, 48 8, 48 7, 47 7)), ((48 8, 48 9, 49 9, 49 8, 48 8))), MULTILINESTRING((47 7, 47 8), (48 8, 48 7)))");
    TEST(getWKT(projWGS842), ==, "GEOMETRYCOLLECTION(MULTIPOINT(7 47,8 48),MULTIPOLYGON(((7 47,8 47,8 48,7 48,7 47)),((8 48,9 48,9 49,8 49,8 48))),MULTILINESTRING((7 47,8 47),(8 48,7 48)))");

    auto projWebMerc2 = util::geo::collectionFromWKT<double>("<http://www.opengis.net/def/crs/EPSG/0/3857> GEOMETRYCOLLECTION(MULTIPOINT(779236.435552915 5942074.072431109, 890555.9263461885 6106854.834885074), MULTIPOLYGON(((779236.435552915 5942074.072431109, 890555.9263461885 5942074.072431109, 890555.9263461885 6106854.834885074, 779236.435552915 6106854.834885074, 779236.435552915 5942074.072431109)), ((890555.9263461885 6106854.834885074, 1001875.4171394621 6106854.834885074, 1001875.4171394621 6274861.394006576, 890555.9263461885 6274861.394006576, 890555.9263461885 6106854.834885074))), MULTILINESTRING((779236.435552915 5942074.072431109, 890555.9263461885 5942074.072431109), (890555.9263461885 6106854.834885074, 779236.435552915 6106854.834885074)))");
    TEST(getWKT(projWebMerc2), ==, "GEOMETRYCOLLECTION(MULTIPOINT(7 47,8 48),MULTIPOLYGON(((7 47,8 47,8 48,7 48,7 47)),((8 48,9 48,9 49,8 49,8 48))),MULTILINESTRING((7 47,8 47),(8 48,7 48)))");
  }
}
