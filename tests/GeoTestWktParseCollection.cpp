// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

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
void GeoTest::testWktParseCollection() {

  {
    TEST(util::geo::getWKTType("GEOMETRYCOLLECTION(MULTIPOLYGON(((1 1,3 3,1 "
                               "1),(0 0,1 1,0 0)),((1 3,3 1,1 3))))"),
         ==, WKTType::COLLECTION);
    TEST(util::geo::getWKTType("GEOMETRYCOLLECTIO (MULTIPOLYGON(((1 1,3 3,1 "
                               "1),(0 0,1 1,0 0)),((1 3,3 1,1 3))))"),
         ==, WKTType::NONE);
    TEST(util::geo::getWKTType(" GEOMETRYCOLLECTION (MULTIPOLYGON(((1 1,3 3,1 "
                               "1),(0 0,1 1,0 0)),((1 3,3 1,1 3))))"),
         ==, WKTType::COLLECTION);
    TEST(util::geo::getWKTType(" GeometryCollection (MULTIPOLYGON(((1 1,3 3,1 "
                               "1),(0 0,1 1,0 0)),((1 3,3 1,1 3))))"),
         ==, WKTType::COLLECTION);
    TEST(util::geo::getWKTType("\t\tGeometryCollection (MULTIPOLYGON(((1 1,3 "
                               "3,1 1),(0 0,1 1,0 0)),((1 3,3 1,1 3))))"),
         ==, WKTType::COLLECTION);
    TEST(util::geo::getWKTType("\t\tmGeometryCollection (MULTIPOLYGON(((1 1,3 "
                               "3,1 1),(0 0,1 1,0 0)),((1 3,3 1,1 3))))"),
         ==, WKTType::COLLECTION);
    TEST(util::geo::getWKTType("\t\tmGeometryCollection M(MULTIPOLYGON(((1 1,3 "
                               "3,1 1),(0 0,1 1,0 0)),((1 3,3 1,1 3))))"),
         ==, WKTType::COLLECTION);
    TEST(util::geo::getWKTType("\t\tmGeometryCollection ZM(MULTIPOLYGON(((1 "
                               "1,3 3,1 1),(0 0,1 1,0 0)),((1 3,3 1,1 3))))"),
         ==, WKTType::COLLECTION);

    TEST(util::geo::getWKT(collectionFromWKT<int>(
             "GEOMETRYCOLLECTION(MULTIPOLYGON(((1 1,3 3,1 1), (0 0,1 1,0 "
             "0)),((1 3,3 1, 1 3))))")),
         ==,
         "GEOMETRYCOLLECTION(MULTIPOLYGON(((1 1,3 3,1 1),(0 0,1 1,0 0)),((1 "
         "3,3 1,1 3))))"); 

    TEST(util::geo::getWKT(collectionFromWKT<int>(
             "GEOMETRYCOLLECTION (MULTIPOLYGON (((1 1,3 3,1 1), (0 0,1 1,0 "
             "0)),((1 3,3 1, 1 3))))")),
         ==,
         "GEOMETRYCOLLECTION(MULTIPOLYGON(((1 1,3 3,1 1),(0 0,1 1,0 0)),((1 "
         "3,3 1,1 3))))");

    TEST(util::geo::getWKT(collectionFromWKT<int>(
             "GEOMETRYCOLLECTION (  MULTIPOLYGON Z (((1 1,3 3,1 1), (0 0,1 1,0 "
             "0)),((1 3,3 1, 1 3))))")),
         ==,
         "GEOMETRYCOLLECTION(MULTIPOLYGON(((1 1,3 3,1 1),(0 0,1 1,0 0)),((1 "
         "3,3 1,1 3))))");

    TEST(util::geo::getWKT(collectionFromWKT<int>(
             "GEOMETRYCOLLECTION (  MMULTIPOLYGON Z (((1 1,3 3,1 1), (0 0,1 "
             "1,0 0)),((1 3,3 1, 1 3))))")),
         ==,
         "GEOMETRYCOLLECTION(MULTIPOLYGON(((1 1,3 3,1 1),(0 0,1 1,0 0)),((1 "
         "3,3 1,1 3))))"); 
  }
}
