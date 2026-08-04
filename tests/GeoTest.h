// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_TEST_GEOTEST_H_
#define UTIL_TEST_GEOTEST_H_

class GeoTest {
  public:
    void run();

  private:
    // run() dispatches to these; the tests were split out of one giant
    // function so each stays small enough to keep -O optimization (esp. GCC's
    // points-to analysis) cheap to compile.
    void testWktParseLine();
    void testWktParseMultiPoint();
    void testWktParsePoint();
    void testWktParsePolygon();
    void testWktParseMultiLine();
    void testWktParseMultiPolygon();
    void testWktParseCollection();
    void testSignedRingArea();
    void testAngles();
    void testDE9IM();
    void testLineLinePredicates();
    void testLinePolygonPredicates();
    void testPolygonPolygonPredicates();
    void testSegmentOrdering();
    void testPolygonInternals();
    void testPolygonFill();
    void testXSortedPolygons();
    void testFrechet();
    void testGrid();
    void testRTree();
    void testLineOps();
    void testWKTParsing();
    void testGeometryBasics();
    void testWktSerialize();
    void testDist();
    void testCRS();
};

#endif
