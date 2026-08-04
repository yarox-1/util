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
void GeoTest::run() {
  testWktParseLine();
  testWktParseMultiPoint();
  testWktParsePoint();
  testWktParsePolygon();
  testWktParseMultiLine();
  testWktParseMultiPolygon();
  testWktParseCollection();
  testSignedRingArea();
  testAngles();
  testDE9IM();
  testLineLinePredicates();
  testLinePolygonPredicates();
  testPolygonPolygonPredicates();
  testSegmentOrdering();
  testPolygonInternals();
  testPolygonFill();
  testXSortedPolygons();
  testFrechet();
  testGrid();
  testRTree();
  testLineOps();
  testWKTParsing();
  testGeometryBasics();
  testWktSerialize();
  testDist();
  testCRS();
}
