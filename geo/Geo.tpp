// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

namespace util {
namespace geo {

const static size_t EST_CHECKS_THRESHOLD_XSORTED = 1000;
const static size_t EST_MULTI_CHECKS_THRESHOLD_XSORTED = 16;

// _____________________________________________________________________________
template <typename T>
Box<T> pad(const Box<T>& box, double xPadding, double yPadding) {
  return Box<T>(Point<T>(boundedSub(box.getLowerLeft().getX(), xPadding),
                         boundedSub(box.getLowerLeft().getY(), yPadding)),
                Point<T>(boundedAdd(box.getUpperRight().getX(), xPadding),
                         boundedAdd(box.getUpperRight().getY(), yPadding)));
}

// _____________________________________________________________________________
template <typename T>
Box<T> pad(const Box<T>& box, double padding) {
  return pad(box, padding, padding);
}

// _____________________________________________________________________________
template <typename T>
RotatedBox<T> pad(const RotatedBox<T>& box, double padding) {
  return {pad(box.getBox(), padding), box.getDegree(), box.getCenter()};
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const Point<T> p) {
  return p;
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const LineSegment<T>& ls) {
  return Point<T>((1.0 * ls.first.getX() + ls.second.getX()) / T(2),
                  (1.0 * ls.first.getY() + ls.second.getY()) / T(2));
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const Line<T>& ls) {
  if (ls.size() == 0) return {0, 0};  // undefined behavior
  if (ls.size() == 1) return ls[0];

  double x = 0, y = 0, sum = 0;
  for (size_t i = 1; i < ls.size(); i++) {
    double l = len(LineSegment<T>{ls[i - 1], ls[i]});
    sum += l;
    x += l * ((ls[i - 1].getX() + ls[i].getX()) / 2);
    y += l * ((ls[i - 1].getY() + ls[i].getY()) / 2);
  }

  if (sum == 0) return ls[0];

  return Point<T>(x / sum, y / sum);
}

// _____________________________________________________________________________
template <typename T>
Point<T> ringCentroid(const Line<T>& ls) {
  if (ls.size() == 0) return {0, 0};  // undefined behavior
  if (ls.size() == 1) return ls[0];

  double xOff = ls[0].getX();
  double yOff = ls[0].getY();

  double x = 0, y = 0, ta = 0;

  size_t j = ls.size() - 1;
  for (size_t i = 0; i < ls.size(); i++) {
    double a = 1.0 * (ls[j].getX() - xOff) * (ls[i].getY() - yOff) -
               1.0 * (ls[i].getX() - xOff) * (ls[j].getY() - yOff);
    ta += a;
    x += a * ((ls[j].getX() - xOff) + (ls[i].getX() - xOff));
    y += a * ((ls[j].getY() - yOff) + (ls[i].getY() - yOff));
    j = i;
  }

  if (ta == 0) return centroid(ls);

  return Point<T>((x / (3.0 * ta)) + xOff, (y / (3.0 * ta)) + yOff);
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const Polygon<T>& o) {
  if (o.getOuter().size() == 0) return {0, 0};  // undefined behavior
  double sumA = 0, x = 0, y = 0;

  double outerArea = ringArea(o.getOuter());
  if (outerArea == 0) outerArea = len(o.getOuter());  // fallback for lower dim

  // if len and area are 0, return single point
  if (outerArea == 0) return o.getOuter().front();

  Point<T> outerCentroid = ringCentroid(o.getOuter());

  sumA += outerArea;
  x += outerCentroid.getX() * outerArea;
  y += outerCentroid.getY() * outerArea;

  for (const auto& inner : o.getInners()) {
    double area = -ringArea(inner);
    if (area == 0) continue;
    Point<T> c = ringCentroid(inner);
    sumA += area;

    x += c.getX() * area;
    y += c.getY() * area;
  }

  return Point<T>(x / sumA, y / sumA);
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const Box<T>& box) {
  return centroid(LineSegment<T>(box.getLowerLeft(), box.getUpperRight()));
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const std::vector<Point<T>>& multigeo) {
  double sumA = 0, x = 0, y = 0;

  for (const auto& g : multigeo) {
    Point<T> c = centroid(g);
    sumA += 1;
    x += c.getX();
    y += c.getY();
  }

  return Point<T>(x / sumA, y / sumA);
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const std::vector<Line<T>>& multigeo) {
  double sumA = 0, x = 0, y = 0;

  for (const auto& g : multigeo) {
    double w = len(g);
    Point<T> c = centroid(g);
    sumA += w;
    x += c.getX() * w;
    y += c.getY() * w;
  }

  return Point<T>(x / sumA, y / sumA);
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const std::vector<LineSegment<T>>& multigeo) {
  double sumA = 0, x = 0, y = 0;

  for (const auto& g : multigeo) {
    double w = len(g);
    Point<T> c = centroid(g);
    sumA += w;
    x += c.getX() * w;
    y += c.getY() * w;
  }

  return Point<T>(x / sumA, y / sumA);
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const std::vector<Polygon<T>>& multigeo) {
  double sumA = 0, x = 0, y = 0;

  if (area(multigeo) > 0) {
    // dim 3
    for (const auto& g : multigeo) {
      double w = area(g);
      Point<T> c = centroid(g);
      sumA += w;
      x += c.getX() * w;
      y += c.getY() * w;
    }
  } else {
    // dim 2
    for (const auto& g : multigeo) {
      double w = len(g);
      Point<T> c = centroid(g);
      sumA += w;
      x += c.getX() * w;
      y += c.getY() * w;
    }
  }

  if (sumA == 0) return centroid(multigeo[0]);

  return Point<T>(x / sumA, y / sumA);
}

// _____________________________________________________________________________
template <typename T>
int8_t dimension(const Collection<T>& col) {
  int8_t dim = 0;
  for (const auto& g : col) {
    if (g.getType() == 0 && dim < 1) dim = 1;
    if (g.getType() == 1 && dim < 2) dim = 2;
    if (g.getType() == 2 && dim < 3) dim = 3;
    if (g.getType() == 3 && dim < 2) dim = 2;
    if (g.getType() == 4 && dim < 3) dim = 3;
    if (g.getType() == 5) {
      int8_t ldim = dimension(g.getCollection());
      if (dim < ldim) dim = ldim;
    }
  }
  return dim;
}

// _____________________________________________________________________________
template <typename T>
Point<T> centroid(const Collection<T>& col) {
  int8_t dim = dimension(col);
  double sum = 0, x = 0, y = 0;

  if (dim == 0) return Point<T>{0, 0};  // undefined behavior

  for (const auto& g : col) {
    double w = 1;
    if (dim == 2) {
      if (g.getType() == 0) w = len(g.getPoint());
      if (g.getType() == 1) w = len(g.getLine());
      if (g.getType() == 2) w = len(g.getPolygon());
      if (g.getType() == 3) w = len(g.getMultiLine());
      if (g.getType() == 4) w = len(g.getMultiPolygon());
      if (g.getType() == 5) w = len(g.getCollection());
    }
    if (dim == 3) {
      if (g.getType() == 0) w = area(g.getPoint());
      if (g.getType() == 1) w = area(g.getLine());
      if (g.getType() == 2) w = area(g.getPolygon());
      if (g.getType() == 3) w = area(g.getMultiLine());
      if (g.getType() == 4) w = area(g.getMultiPolygon());
      if (g.getType() == 5) w = area(g.getCollection());
    }

    sum += w;
    Point<T> cen;
    if (g.getType() == 0) cen = centroid(g.getPoint());
    if (g.getType() == 1) cen = centroid(g.getLine());
    if (g.getType() == 2) cen = centroid(g.getPolygon());
    if (g.getType() == 3) cen = centroid(g.getMultiLine());
    if (g.getType() == 4) cen = centroid(g.getMultiPolygon());
    if (g.getType() == 5) cen = centroid(g.getCollection());

    x += cen.getX() * w;
    y += cen.getY() * w;
  }

  return Point<T>{x / sum, y / sum};
}

// _____________________________________________________________________________
template <typename T>
Point<T> rotate(const Point<T>& p, double) {
  return p;
}

// _____________________________________________________________________________
template <typename T>
Point<T> rotateRAD(const Point<T>& p, double) {
  return p;
}

// _____________________________________________________________________________
template <typename T>
Point<T> rotateSinCos(Point<T> p, double si, double co, const Point<T>& c) {
  p = p - c;

  return Point<T>(p.getX() * co - p.getY() * si,
                  p.getX() * si + p.getY() * co) +
         c;
}

// _____________________________________________________________________________
template <typename T>
Point<T> rotateRAD(Point<T> p, double rad, const Point<T>& c) {
  return rotateSinCos(p, sin(rad), cos(rad), c);
}

// _____________________________________________________________________________
template <typename T>
Point<T> rotate(Point<T> p, double deg, const Point<T>& c) {
  return rotateRAD(p, deg * -RAD, c);
}

// _____________________________________________________________________________
template <typename T>
LineSegment<T> rotate(LineSegment<T> geo, double deg, const Point<T>& c) {
  geo.first = rotate(geo.first, deg, c);
  geo.second = rotate(geo.second, deg, c);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
LineSegment<T> rotateRAD(LineSegment<T> geo, double deg, const Point<T>& c) {
  geo.first = rotateRAD(geo.first, deg, c);
  geo.second = rotateRAD(geo.second, deg, c);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
LineSegment<T> rotate(LineSegment<T> geo, double deg) {
  return rotate(geo, deg, centroid(geo));
}

// _____________________________________________________________________________
template <typename T>
LineSegment<T> rotateRAD(LineSegment<T> geo, double deg) {
  return rotateRAD(geo, deg, centroid(geo));
}

// _____________________________________________________________________________
template <typename T>
Line<T> rotateSinCos(Line<T> geo, double si, double co, const Point<T>& c) {
  for (size_t i = 0; i < geo.size(); i++)
    geo[i] = rotateSinCos(geo[i], si, co, c);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
Line<T> rotateRAD(Line<T> geo, double deg, const Point<T>& c) {
  const double si = sin(deg);
  const double co = cos(deg);
  return rotateSinCos(geo, si, co, c);
}

// _____________________________________________________________________________
template <typename T>
Line<T> rotate(Line<T> geo, double deg, const Point<T>& c) {
  return rotateRAD(geo, deg * -RAD, c);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> rotateSinCos(Polygon<T> geo, double si, double co,
                        const Point<T>& c) {
  for (size_t i = 0; i < geo.getOuter().size(); i++)
    geo.getOuter()[i] = rotateSinCos(geo.getOuter()[i], si, co, c);
  for (size_t i = 0; i < geo.getInners().size(); i++)
    for (size_t j = 0; j < geo.getInners()[i].size(); j++)
      geo.getInners()[i][j] = rotateSinCos(geo.getInners()[i][j], si, co, c);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> rotateRAD(Polygon<T> geo, double deg, const Point<T>& c) {
  const double si = sin(deg);
  const double co = cos(deg);

  return rotateSinCos(geo, si, co, c);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> rotate(Polygon<T> geo, double deg, const Point<T>& c) {
  return rotateRAD(geo, deg * -RAD, c);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
std::vector<Geometry<T>> rotate(std::vector<Geometry<T>> multigeo, double deg,
                                const Point<T>& c) {
  for (size_t i = 0; i < multigeo.size(); i++)
    multigeo[i] = rotate(multigeo[i], deg, c);
  return multigeo;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
std::vector<Geometry<T>> rotateRAD(std::vector<Geometry<T>> multigeo,
                                   double deg, const Point<T>& c) {
  for (size_t i = 0; i < multigeo.size(); i++)
    multigeo[i] = rotateRAD(multigeo[i], deg, c);
  return multigeo;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
std::vector<Geometry<T>> rotate(std::vector<Geometry<T>> multigeo, double deg) {
  auto c = centroid(multigeo);
  for (size_t i = 0; i < multigeo.size(); i++)
    multigeo[i] = rotate(multigeo[i], deg, c);
  return multigeo;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
std::vector<Geometry<T>> rotateRAD(std::vector<Geometry<T>> multigeo,
                                   double deg) {
  auto c = centroid(multigeo);
  for (size_t i = 0; i < multigeo.size(); i++)
    multigeo[i] = rotateRAD(multigeo[i], deg, c);
  return multigeo;
}

// _____________________________________________________________________________
template <typename T>
Point<T> move(const Point<T>& geo, double x, double y) {
  return Point<T>(geo.getX() + x, geo.getY() + y);
}

// _____________________________________________________________________________
template <typename T>
Line<T> move(Line<T> geo, double x, double y) {
  for (size_t i = 0; i < geo.size(); i++) geo[i] = move(geo[i], x, y);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
LineSegment<T> move(LineSegment<T> geo, double x, double y) {
  geo.first = move(geo.first, x, y);
  geo.second = move(geo.second, x, y);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> move(Polygon<T> geo, double x, double y) {
  for (size_t i = 0; i < geo.getOuter().size(); i++)
    geo.getOuter()[i] = move(geo.getOuter()[i], x, y);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
Box<T> move(Box<T> geo, double x, double y) {
  geo.setLowerLeft(move(geo.getLowerLeft(), x, y));
  geo.setUpperRight(move(geo.getUpperRight(), x, y));
  return geo;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
std::vector<Geometry<T>> move(std::vector<Geometry<T>> multigeo, double x,
                              double y) {
  for (size_t i = 0; i < multigeo.size(); i++)
    multigeo[i] = move(multigeo[i], x, y);
  return multigeo;
}

// _____________________________________________________________________________
template <typename T>
Box<T> minbox() {
  return Box<T>();
}

// _____________________________________________________________________________
template <typename T>
RotatedBox<T> shrink(const RotatedBox<T>& b, double d) {
  double xd =
      b.getBox().getUpperRight().getX() - b.getBox().getLowerLeft().getX();
  double yd =
      b.getBox().getUpperRight().getY() - b.getBox().getLowerLeft().getY();

  if (xd <= 2 * d) d = xd / 2 - 1;
  if (yd <= 2 * d) d = yd / 2 - 1;

  Box<T> r(Point<T>(b.getBox().getLowerLeft().getX() + d,
                    b.getBox().getLowerLeft().getY() + d),
           Point<T>(b.getBox().getUpperRight().getX() - d,
                    b.getBox().getUpperRight().getY() - d));

  return RotatedBox<T>(r, b.getDegree(), b.getCenter());
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Point<T>& p, uint16_t prec) {
  std::string ret = "POINT(";
  ret.reserve(6 + prec + 3 + prec + 3 + 1);
  ret.append(formatFloat(p.getX(), prec));
  ret.push_back(' ');
  ret.append(formatFloat(p.getY(), prec));
  ret.push_back(')');
  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Point<T>& p) {
  return getWKT(p, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const std::vector<Point<T>>& p, uint16_t prec) {
  std::string ret = "MULTIPOINT(";
  ret.reserve(10 + 1 + p.size() * (prec + 3) * 2 + 1);
  for (size_t i = 0; i < p.size(); i++) {
    if (i) ret.push_back(',');
    ret.append(formatFloat(p[i].getX(), prec));
    ret.push_back(' ');
    ret.append(formatFloat(p[i].getY(), prec));
  }
  ret.push_back(')');
  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const std::vector<Point<T>>& p) {
  return getWKT(p, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Line<T>& l, uint16_t prec) {
  std::string ret = "LINESTRING(";
  ret.reserve(10 + 1 + l.size() * (prec + 3) * 2 + 1);
  for (size_t i = 0; i < l.size(); i++) {
    if (i) ret.push_back(',');
    ret.append(formatFloat(l[i].getX(), prec));
    ret.push_back(' ');
    ret.append(formatFloat(l[i].getY(), prec));
  }
  ret.push_back(')');
  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Line<T>& l) {
  return getWKT(l, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const std::vector<Line<T>>& ls, uint16_t prec) {
  std::string ret = "MULTILINESTRING(";

  if (ls.size()) ret.reserve(15 + 2 + ls[0].size() * (prec + 3) * 2 + 2);

  for (size_t j = 0; j < ls.size(); j++) {
    if (j) ret.push_back(',');
    ret.push_back('(');
    for (size_t i = 0; i < ls[j].size(); i++) {
      if (i) ret.push_back(',');
      ret.append(formatFloat(ls[j][i].getX(), prec));
      ret.push_back(' ');
      ret.append(formatFloat(ls[j][i].getY(), prec));
    }
    ret.push_back(')');
  }

  ret.push_back(')');
  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const std::vector<Line<T>>& ls) {
  return getWKT(ls, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const XSortedPolygon<T>& ls, uint16_t prec) {
  std::stringstream ss;
  ss << "MULTILINESTRING(";

  for (size_t j = 0; j < ls.getOuter().rawRing().size(); j++) {
    if (ls.getOuter().rawRing()[j].out()) continue;
    if (j) ss << ",";
    ss << "(";
    ss << formatFloat(ls.getOuter().rawRing()[j].seg().first.getX() * 1.0 / 10,
                      prec)
       << " "
       << formatFloat(ls.getOuter().rawRing()[j].seg().first.getY() * 1.0 / 10,
                      prec)
       << ",";
    ss << formatFloat(ls.getOuter().rawRing()[j].seg().second.getX() * 1.0 / 10,
                      prec)
       << " "
       << formatFloat(ls.getOuter().rawRing()[j].seg().second.getY() * 1.0 / 10,
                      prec);
    ss << ")";
  }

  ss << ")";
  return ss.str();
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const XSortedPolygon<T>& ls) {
  return getWKT(ls, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const LineSegment<T>& l, uint16_t prec) {
  return getWKT(Line<T>{l.first, l.second}, prec);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const LineSegment<T>& l) {
  return getWKT(l, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Box<T>& l, uint16_t prec) {
  std::string ret = "POLYGON((";
  ret.reserve(7 + 2 + 4 * (prec + 3) * 2 + 2);

  ret.append(formatFloat(l.getLowerLeft().getX(), prec));
  ret.push_back(' ');
  ret.append(formatFloat(l.getLowerLeft().getY(), prec));
  ret.push_back(',');
  ret.append(formatFloat(l.getUpperRight().getX(), prec));
  ret.push_back(' ');
  ret.append(formatFloat(l.getLowerLeft().getY(), prec));
  ret.push_back(',');
  ret.append(formatFloat(l.getUpperRight().getX(), prec));
  ret.push_back(' ');
  ret.append(formatFloat(l.getUpperRight().getY(), prec));
  ret.push_back(',');
  ret.append(formatFloat(l.getLowerLeft().getX(), prec));
  ret.push_back(' ');
  ret.append(formatFloat(l.getUpperRight().getY(), prec));
  ret.push_back(',');
  ret.append(formatFloat(l.getLowerLeft().getX(), prec));
  ret.push_back(' ');
  ret.append(formatFloat(l.getLowerLeft().getY(), prec));
  ret.append("))");
  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Box<T>& l) {
  return getWKT(l, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Polygon<T>& p, uint16_t prec) {
  if (p.getOuter().size() == 0) return "POLYGON()";
  std::string ret = "POLYGON((";
  ret.reserve(7 + 2 + p.getOuter().size() * (prec + 3) * 2 + 2);

  for (size_t i = 0; i < p.getOuter().size(); i++) {
    if (i > 0) ret.push_back(',');
    ret.append(formatFloat(p.getOuter()[i].getX(), prec));
    ret.push_back(' ');
    ret.append(formatFloat(p.getOuter()[i].getY(), prec));
  }

  if (p.getOuter().front() != p.getOuter().back()) {
    ret.push_back(',');
    ret.append(formatFloat(p.getOuter().front().getX(), prec));
    ret.push_back(' ');
    ret.append(formatFloat(p.getOuter().front().getY(), prec));
  }
  ret.push_back(')');

  for (const auto& inner : p.getInners()) {
    ret.append(",(");
    for (size_t i = 0; i < inner.size(); i++) {
      if (i > 0) ret.push_back(',');
      ret.append(formatFloat(inner[i].getX(), prec));
      ret.push_back(' ');
      ret.append(formatFloat(inner[i].getY(), prec));
    }

    if (inner.front() != inner.back()) {
      ret.push_back(',');
      ret.append(formatFloat(inner.front().getX(), prec));
      ret.push_back(' ');
      ret.append(formatFloat(inner.front().getY(), prec));
    }
    ret.push_back(')');
  }
  ret.push_back(')');
  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Polygon<T>& p) {
  return getWKT(p, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const std::vector<Polygon<T>>& ls, uint16_t prec) {
  std::string ret = "MULTIPOLYGON(";
  if (ls.size())
    ret.reserve(12 + 2 + ls[0].getOuter().size() * (prec + 3) * 2 + 2);

  for (size_t j = 0; j < ls.size(); j++) {
    if (j) ret.push_back(',');
    ret.push_back('(');
    ret.push_back('(');
    for (size_t i = 0; i < ls[j].getOuter().size(); i++) {
      if (i > 0) ret.push_back(',');
      ret.append(formatFloat(ls[j].getOuter()[i].getX(), prec));
      ret.push_back(' ');
      ret.append(formatFloat(ls[j].getOuter()[i].getY(), prec));
    }

    if (ls[j].getOuter().front() != ls[j].getOuter().back()) {
      ret.push_back(',');
      ret.append(formatFloat(ls[j].getOuter().front().getX(), prec));
      ret.push_back(' ');
      ret.append(formatFloat(ls[j].getOuter().front().getY(), prec));
    }
    ret.push_back(')');

    for (const auto& inner : ls[j].getInners()) {
      ret.push_back(',');
      ret.push_back('(');
      for (size_t i = 0; i < inner.size(); i++) {
        if (i > 0) ret.push_back(',');
        ret.append(formatFloat(inner[i].getX(), prec));
        ret.push_back(' ');
        ret.append(formatFloat(inner[i].getY(), prec));
      }
      if (inner.front() != inner.back()) {
        ret.push_back(',');
        ret.append(formatFloat(inner.front().getX(), prec));
        ret.push_back(' ');
        ret.append(formatFloat(inner.front().getY(), prec));
      }
      ret.push_back(')');
    }
    ret.push_back(')');
  }

  ret.push_back(')');
  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const std::vector<Polygon<T>>& ls) {
  return getWKT(ls, 6);
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Collection<T>& coll, uint16_t prec) {
  std::string ret = "GEOMETRYCOLLECTION(";

  std::string delim = "";

  for (const auto& g : coll) {
    ret += delim;
    delim = ",";
    if (g.getType() == 0) ret += util::geo::getWKT(g.getPoint(), prec);
    if (g.getType() == 1) ret += util::geo::getWKT(g.getLine(), prec);
    if (g.getType() == 2) ret += util::geo::getWKT(g.getPolygon(), prec);
    if (g.getType() == 3) ret += util::geo::getWKT(g.getMultiLine(), prec);
    if (g.getType() == 4) ret += util::geo::getWKT(g.getMultiPolygon(), prec);
    if (g.getType() == 5) ret += util::geo::getWKT(g.getCollection(), prec);
  }

  return ret + ")";
}

// _____________________________________________________________________________
template <typename T>
std::string getWKT(const Collection<T>& coll) {
  return getWKT(coll, 6);
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Point<T>& p, const Box<T>& box) {
  // check if point lies in box
  return (fabs(p.getX() - box.getLowerLeft().getX()) < EPSILON ||
          p.getX() > box.getLowerLeft().getX()) &&
         (fabs(p.getX() - box.getUpperRight().getX()) < EPSILON ||
          p.getX() < box.getUpperRight().getX()) &&
         (fabs(p.getY() - box.getLowerLeft().getY()) < EPSILON ||
          p.getY() > box.getLowerLeft().getY()) &&
         (fabs(p.getY() - box.getUpperRight().getY()) < EPSILON ||
          p.getY() < box.getUpperRight().getY());
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Line<T>& l, const Box<T>& box) {
  // check if line lies in box
  for (const auto& p : l)
    if (!contains(p, box)) return false;
  return true;
}

// _____________________________________________________________________________
template <typename T>
bool contains(const LineSegment<T>& l, const Box<T>& box) {
  // check if line segment lies in box
  return contains(l.first, box) && contains(l.second, box);
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Box<T>& b, const Box<T>& box) {
  // check if box b lies in box
  return contains(b.getLowerLeft(), box) && contains(b.getUpperRight(), box);
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Point<T>& p, const LineSegment<T>& ls) {
  // check if point p lies in (on) line segment ls
  return fabs(crossProd(p, ls)) < EPSILON && contains(p, getBoundingBox(ls));
}

// _____________________________________________________________________________
template <typename T>
bool contains(const LineSegment<T>& a, const LineSegment<T>& b) {
  // check if line segment a is contained in line segment b
  return contains(a.first, b) && contains(a.second, b);
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Point<T>& p, const Line<T>& l) {
  // check if point p lies in line l
  for (size_t i = 1; i < l.size(); i++) {
    if (contains(p, LineSegment<T>(l[i - 1], l[i]))) return true;
  }
  return false;
}

// _____________________________________________________________________________
template <typename T>
bool ringContains(const Point<T>& p, const Ring<T>& ph) {
  // check if point p lies in polygon hull

  // see https://de.wikipedia.org/wiki/Punkt-in-Polygon-Test_nach_Jordan
  int8_t c = -1;

  for (size_t i = 1; i < ph.size(); i++) {
    c *= polyContCheck(p, ph[i - 1], ph[i]);
    if (c == 0) return true;
  }

  c *= polyContCheck(p, ph.back(), ph[0]);

  return c >= 0;
}

// _____________________________________________________________________________
template <typename T>
int8_t polyContCheck(const Point<T>& a, Point<T> b, Point<T> c) {
  if (a.getY() == b.getY() && a.getY() == c.getY())
    return (!((b.getX() <= a.getX() && a.getX() <= c.getX()) ||
              (c.getX() <= a.getX() && a.getX() <= b.getX())));

  if (a.getX() == b.getX() && a.getX() == c.getX()) {
    return (!((b.getY() <= a.getY() && a.getY() <= c.getY()) ||
              (c.getY() <= a.getY() && a.getY() <= b.getY())));
  }

  if (fabs(a.getY() - b.getY()) < EPSILON &&
      fabs(a.getX() - b.getX()) < EPSILON)
    return 0;

  if (b.getX() > c.getX()) {
    Point<T> tmp = b;
    b = c;
    c = tmp;
  }

  if (a.getX() <= b.getX() || a.getX() > c.getX()) {
    return 1;
  }

  double d =
      (1.0 * b.getY() - 1.0 * a.getY()) * (1.0 * c.getX() - 1.0 * a.getX()) -
      (1.0 * b.getX() - 1.0 * a.getX()) * (1.0 * c.getY() - 1.0 * a.getY());
  if (d > 0) return -1;
  if (d < 0) return 1;
  return 0;
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Point<T>& p, const Polygon<T>& poly) {
  if (!ringContains(p, poly.getOuter())) return false;

  for (const auto& inner : poly.getInners()) {
    if (ringContains(p, inner)) return false;
  }

  return true;
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
void _withinDistCheck(const XSortedCollection<T>& a, util::geo::BoxVal<T> ab,
                      const XSortedCollection<T>& b, util::geo::BoxVal<T> bb,
                      double& minDist, PF&& paddingFunc,
                      double& maxEuclideanDist, DF&& distFunc) {
  if (util::geo::dist(a.boundingBox(ab.type, ab.id),
                      b.boundingBox(bb.type, bb.id)) > maxEuclideanDist)
    return;

  // padding function which keeps track of the upper euclidean distance bound
  // found so far, done this way here so that we do not have to change signature
  // of all individual withinDist() function to also optionally return the
  // euclidean distance upper bound
  std::function<double(double, double, const Box<T>&, const Box<T>&)>
      capturingPaddingFunc = [&](double euclideanDist, double dist,
                                 const Box<T>& ba, const Box<T>& bb) -> double {
    maxEuclideanDist = std::min(maxEuclideanDist, euclideanDist);
    return paddingFunc(euclideanDist, dist, ba, bb);
  };

  if (ab.type == SWEEP_POLYGON && bb.type == SWEEP_POLYGON) {
    minDist = std::min(
        minDist, withinDist(a.getPolygon(ab.id), b.getPolygon(bb.id), minDist,
                            capturingPaddingFunc, maxEuclideanDist, distFunc));
  }
  if (ab.type == SWEEP_LINESTRING && bb.type == SWEEP_LINESTRING) {
    minDist = std::min(
        minDist, withinDist(a.getLine(ab.id), b.getLine(bb.id), minDist,
                            capturingPaddingFunc, maxEuclideanDist, distFunc));
  }
  if (ab.type == SWEEP_POINT && bb.type == SWEEP_POINT) {
    minDist =
        std::min(minDist, dist(a.getPoint(ab.id), b.getPoint(bb.id), distFunc));
  }
  if (ab.type == SWEEP_POINT && bb.type == SWEEP_LINESTRING) {
    minDist = std::min(
        minDist, withinDist(a.getPoint(ab.id), b.getLine(bb.id), minDist,
                            capturingPaddingFunc, maxEuclideanDist, distFunc));
  }
  if (ab.type == SWEEP_LINESTRING && bb.type == SWEEP_POINT) {
    minDist = std::min(
        minDist, withinDist(a.getLine(ab.id), b.getPoint(bb.id), minDist,
                            capturingPaddingFunc, maxEuclideanDist, distFunc));
  }
  if (ab.type == SWEEP_POLYGON && bb.type == SWEEP_LINESTRING) {
    minDist = std::min(
        minDist, withinDist(a.getPolygon(ab.id), b.getLine(bb.id), minDist,
                            capturingPaddingFunc, maxEuclideanDist, distFunc));
  }
  if (ab.type == SWEEP_LINESTRING && bb.type == SWEEP_POLYGON) {
    minDist = std::min(
        minDist, withinDist(a.getLine(ab.id), b.getPolygon(bb.id), minDist,
                            capturingPaddingFunc, maxEuclideanDist, distFunc));
  }
  if (ab.type == SWEEP_POINT && bb.type == SWEEP_POLYGON) {
    minDist = std::min(
        minDist, withinDist(a.getPoint(ab.id), b.getPolygon(bb.id), minDist,
                            capturingPaddingFunc, maxEuclideanDist, distFunc));
  }
  if (ab.type == SWEEP_POLYGON && bb.type == SWEEP_POINT) {
    minDist = std::min(
        minDist, withinDist(a.getPolygon(ab.id), b.getPoint(bb.id), minDist,
                            capturingPaddingFunc, maxEuclideanDist, distFunc));
  }
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
void _withinDistCheck(const XSortedCollection<T>& a,
                      const XSortedCollection<T>& b, double& minDist,
                      PF&& paddingFunc, double& maxEuclideanDist, DF&& distFunc,
                      const BoxVal<T>& boxVal,
                      const util::geo::IntervalIdx<T, BoxVal<T>>& actives) {
  if (util::geo::dist(a.boundingBox(boxVal.type, boxVal.id), b.boundingBox()) >
      maxEuclideanDist)
    return;
  for (const auto& other : actives.overlap_find_all({boxVal.loY, boxVal.upY})) {
    _withinDistCheck(a, boxVal, b, other.v, minDist, paddingFunc,
                     maxEuclideanDist, distFunc);
  }
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF,
          template <typename> class XSORTED>
double withinDist(const XSortedCollection<T>& a, const XSORTED<T>& b,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  // TODO: might be faster to avoid the Collection construction here
  return withinDist(a, XSortedCollection<T>(b), maxDist, paddingFunc,
                    maxEuclideanDist, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const XSortedCollection<T>& a, const XSortedCollection<T>& b,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  double minDist = nextafter(maxDist, std::numeric_limits<double>::infinity());

  if (a.size() == 0 || b.size() == 0)
    return std::numeric_limits<double>::infinity();

  const auto& boxA = a.boundingBox();
  const auto& boxB = b.boundingBox();

  if (util::geo::dist(boxA, boxB) > maxEuclideanDist) {
    return std::numeric_limits<double>::max();
  }

  const auto& ea = a.sweepEvents();
  auto eb = b.sweepEvents();

  // simple probing between extreme points in sorted lists, also ensures that
  // neither maxDist nor maxEuclideanDist is infinite
  for (size_t i = 0; i < ea.size(); i += ea.size() - 1) {
    for (size_t j = 0; j < eb.size(); j += eb.size() - 1) {
      _withinDistCheck(a, ea[i], b, eb[j], minDist, paddingFunc,
                       maxEuclideanDist, distFunc);

      // early abort
      if (minDist == 0) return 0;
    }
  }

  if (ea.size() * eb.size() == 1) {
    // probing already found distance
    return minDist;
  }

  double euclideanDistUpperBound = maxEuclideanDist;
  double distanceUpperBound = minDist;

  // minimum euclidean X dist based on bounding boxes for better y padding
  double euclideanXDistLowerBound =
      dist(LineSegment<T>{Point<T>{boxA.getLowerLeft().getX(), 0},
                          Point<T>{boxA.getUpperRight().getX(), 0}},
           LineSegment<T>{Point<T>{boxB.getLowerLeft().getX(), 0},
                          Point<T>{boxB.getUpperRight().getX(), 0}});
  double euclideanYDistLowerBound =
      dist(LineSegment<T>{Point<T>{0, boxA.getLowerLeft().getY()},
                          Point<T>{0, boxA.getUpperRight().getY()}},
           LineSegment<T>{Point<T>{0, boxB.getLowerLeft().getY()},
                          Point<T>{0, boxB.getUpperRight().getY()}});

  T padding = std::min(
      std::numeric_limits<T>::max() * 1.0,
      paddingFunc(euclideanDistUpperBound, distanceUpperBound, boxA, boxB) *
          1.0);

  T xPadding =
      std::min(std::numeric_limits<T>::max() * 1.0,
               (sqrt(std::max(0.0, padding * 1.0 * padding * 1.0 -
                                       euclideanYDistLowerBound * 1.0 *
                                           euclideanYDistLowerBound * 1.0))));
  T yPadding =
      std::min(std::numeric_limits<T>::max() * 1.0,
               (sqrt(std::max(0.0, padding * 1.0 * padding * 1.0 -
                                       euclideanXDistLowerBound * 1.0 *
                                           euclideanXDistLowerBound * 1.0))));

  // pad the bounding boxes of B, ensure that we stay within the limits of T
  for (auto& b : eb) {
    if (!b.out) {
      b.x = boundedSub(b.x, xPadding);
    } else {
      b.x = boundedAdd(b.x, xPadding);
    }
    b.loY = boundedSub(b.loY, yPadding);
    b.upY = boundedAdd(b.upY, yPadding);
  }

  std::sort(eb.begin(), eb.end());

  util::geo::IntervalIdx<T, BoxVal<T>> activesA;
  util::geo::IntervalIdx<T, BoxVal<T>> activesB;

  size_t i = 0;
  size_t j = 0;

  while (i < ea.size() && j < eb.size()) {
    if (ea[i] < eb[j]) {
      if (!ea[i].out) {
        activesA.insert({ea[i].loY, ea[i].upY}, ea[i]);
      } else {
        activesA.erase({ea[i].loY, ea[i].upY}, ea[i]);
        _withinDistCheck(a, b, minDist, paddingFunc, maxEuclideanDist, distFunc,
                         ea[i], activesB);
      }

      i++;
    } else if (eb[j] < ea[i]) {
      if (!eb[j].out) {
        activesB.insert({eb[j].loY, eb[j].upY}, eb[j]);
      } else {
        activesB.erase({eb[j].loY, eb[j].upY}, eb[j]);
        _withinDistCheck(b, a, minDist, paddingFunc, maxEuclideanDist, distFunc,
                         eb[j], activesA);
      }

      j++;
    } else {
      // equal
      if (!ea[i].out) {
        activesA.insert({ea[i].loY, ea[i].upY}, ea[i]);
      } else {
        activesA.erase({ea[i].loY, ea[i].upY}, ea[i]);
        _withinDistCheck(a, b, minDist, paddingFunc, maxEuclideanDist, distFunc,
                         ea[i], activesB);
      }

      if (!eb[j].out) {
        activesB.insert({eb[j].loY, eb[j].upY}, eb[j]);
      } else {
        activesB.erase({eb[j].loY, eb[j].upY}, eb[j]);
        _withinDistCheck(b, a, minDist, paddingFunc, maxEuclideanDist, distFunc,
                         eb[j], activesA);
      }

      i++;
      j++;
    }

    // early stop
    if (minDist == 0) return 0;
  }

  return minDist;
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF,
          template <typename> class Geometry>
double withinDist(const Collection<T>& a, const Geometry<T>& b, double maxDist,
                  PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  return withinDist(b, XSortedCollection<T>(a), maxDist, paddingFunc,
                    maxEuclideanDist, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF,
          template <typename> class Geometry>
double withinDist(const Collection<T>& a, const std::vector<Geometry<T>>& b,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  return withinDist(XSortedCollection<T>(b), XSortedCollection<T>(a), maxDist,
                    paddingFunc, maxEuclideanDist, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF,
          template <typename> class Geometry>
double withinDist(const std::vector<Geometry<T>>& b, const Collection<T>& a,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  return withinDist(XSortedCollection<T>(b), XSortedCollection<T>(a), maxDist,
                    paddingFunc, maxEuclideanDist, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const Collection<T>& a, const Collection<T>& b,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  return withinDist(XSortedCollection<T>(a), XSortedCollection<T>(b), maxDist,
                    paddingFunc, maxEuclideanDist, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
std::pair<double, bool> withinDist(const Point<T>& p, const XSortedRing<T>& ph,
                                   double maxDist, PF&& paddingFunc,
                                   double maxEuclideanDist, DF&& distFunc) {
  // check if point p is within distance "maxDist" to ring, while also only
  // checking elements that are within the euclidean distance "maxEuclideanDist"
  //
  // returns {distance, contained}

  if (ph.rawRing().size() == 0)
    return {std::numeric_limits<double>::max(), false};

  int8_t c = -1;

  size_t i = 0;

  // do some distance probing to get a lower upper bound
  double euclideanDistUpperBound =
      std::min(maxEuclideanDist, std::min(dist(p, ph.rawRing().back().other),
                                          dist(p, ph.rawRing().front().p)));
  double minDist = std::min(
      distFunc(p, ph.rawRing().back().other,
               std::numeric_limits<double>::max()),
      distFunc(p, ph.rawRing().front().p, std::numeric_limits<double>::max()));
  T padding = std::min(std::numeric_limits<T>::max() * 1.0,
                       paddingFunc(euclideanDistUpperBound, minDist,
                                   getBoundingBox(p), ph.boundingBox()) *
                           1.0);

  // skip irrelevant parts in poly
  if (ph.getMaxSegLen() < std::numeric_limits<T>::max() &&
      i < ph.rawRing().size() &&
      ph.rawRing()[i].p.getX() <
          boundedSub(boundedSub(p.getX(), ph.getMaxSegLen()), padding)) {
    i = std::lower_bound(
            ph.rawRing().begin() + i, ph.rawRing().end(),
            XSortedTuple<T>{
                {boundedSub(boundedSub(p.getX(), ph.getMaxSegLen()), padding),
                 0},
                false}) -
        ph.rawRing().begin();
  }

  while (i < ph.rawRing().size() &&
         ph.rawRing()[i].seg().second.getX() < boundedSub(p.getX(), padding))
    i++;

  for (; i < ph.rawRing().size(); i++) {
    if (ph.rawRing()[i].out()) continue;
    // there won't be coming any more lines intersecting a straight north/south
    // line through p
    if (boundedSub(ph.rawRing()[i].seg().first.getX(), padding) > p.getX())
      break;
    c *= polyContCheck(p, ph.rawRing()[i].seg().first,
                       ph.rawRing()[i].seg().second);
    if (c == 0) return {0, false};

    double euclideanDist = dist(ph.rawRing()[i].seg(), p);
    if (euclideanDist <= padding) {
      auto p2 = projectOn(ph.rawRing()[i].seg().first, p,
                          ph.rawRing()[i].seg().second);

      double dist = distFunc(p, p2, std::min(maxDist, minDist));
      if (dist <= maxDist && dist < minDist) minDist = dist;
      if (minDist == 0) return {0, false};
    }
    if (euclideanDist < euclideanDistUpperBound) {
      euclideanDistUpperBound = euclideanDist;
      padding = paddingFunc(euclideanDistUpperBound, minDist, getBoundingBox(p),
                            ph.boundingBox());
    }
  }

  return {minDist, c >= 0};
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const Point<T>& p, const XSortedPolygon<T>& poly,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  // first do the check in the outer boundary
  auto r = withinDist(p, poly.getOuter(), maxDist, paddingFunc,
                      maxEuclideanDist, distFunc);

  // if we are not included in the outer ring,  we can abort here
  if (!r.second) return r.first;

  // if we are included in the outer ring, we have to check the inner rings too
  if (r.second && poly.getInners().size()) {
    size_t i = 0;

    // skip irrelevant inner rings by their bounding box
    // NOTE: we are only interested in inner polygons which contain
    // the point, so no need to add the maxEuclideanDist here to anything
    if (poly.getInnerMaxSegLen() < std::numeric_limits<double>::infinity()) {
      i = std::lower_bound(
              poly.getInnerBoxIdx().begin(), poly.getInnerBoxIdx().end(),
              std::pair<T, size_t>{
                  boundedSub(p.getX(), poly.getInnerMaxSegLen()), 0}) -
          poly.getInnerBoxIdx().begin();
    }

    for (; i < poly.getInners().size(); i++) {
      // again skip any inner rings if p is not contained in their boxes
      if (poly.getInnerBoxes()[i].getLowerLeft().getX() > p.getX()) break;
      if (!util::geo::contains(p, poly.getInnerBoxes()[i])) continue;

      auto r = withinDist(p, poly.getInners()[i], maxDist, paddingFunc,
                          maxEuclideanDist, distFunc);

      // if we are contained in the inner ring, directly return the distance to
      // its border. We can safely abort has as we assume that inner rings never
      // intersect
      if (r.second) return r.first;
    }
  }

  return r.second ? 0 : r.first;
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const XSortedPolygon<T>& poly, const Point<T>& p,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  return withinDist(p, poly, maxDist, paddingFunc, maxEuclideanDist, distFunc);
}

// _____________________________________________________________________________
template <typename T>
std::pair<bool, bool> ringContains(const Point<T>& p, const XSortedRing<T>& ph,
                                   size_t i) {
  if (!contains(p, ph.boundingBox())) return {false, false};
  // returns {contains, covers}
  // check if point p lies in polygon

  // see https://de.wikipedia.org/wiki/Punkt-in-Polygon-Test_nach_Jordan
  int8_t c = -1;

  // skip irrelevant parts in poly
  if (ph.getMaxSegLen() < std::numeric_limits<double>::infinity() &&
      i < ph.rawRing().size() &&
      ph.rawRing()[i].p.getX() < boundedSub(p.getX(), ph.getMaxSegLen())) {
    i = std::lower_bound(
            ph.rawRing().begin() + i, ph.rawRing().end(),
            XSortedTuple<T>{{boundedSub(p.getX(), ph.getMaxSegLen()), 0},
                            false}) -
        ph.rawRing().begin();
  }

  while (i < ph.rawRing().size() &&
         ph.rawRing()[i].seg().second.getX() < p.getX())
    i++;

  for (; i < ph.rawRing().size(); i++) {
    if (ph.rawRing()[i].out()) continue;
    // there won't be coming any more lines intersecting a straight north/south
    // line through p
    if (ph.rawRing()[i].seg().first.getX() > p.getX()) break;
    c *= polyContCheck(p, ph.rawRing()[i].seg().first,
                       ph.rawRing()[i].seg().second);
    if (c == 0) return {0, 1};
  }

  return {c >= 0, c >= 0};
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const Point<T>& p, const XSortedLine<T>& line, double maxDist,
                  PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  // check if point p lies on line

  if (line.rawLine().size() == 0) return std::numeric_limits<double>::max();

  size_t i = 0;

  // do some distance probing to get a lower upper bound
  double euclideanDistUpperBound =
      std::min(maxEuclideanDist, std::min(dist(p, line.rawLine().back().other),
                                          dist(p, line.rawLine().front().p)));
  double minDist = std::min(
      std::min(maxDist, distFunc(p, line.rawLine().back().other, maxDist)),
      distFunc(p, line.rawLine().front().p, maxDist));

  auto padding = paddingFunc(euclideanDistUpperBound, minDist,
                             getBoundingBox(p), line.boundingBox());

  // skip irrelevant parts
  if (line.getMaxSegLen() < std::numeric_limits<double>::infinity() &&
      i < line.rawLine().size()) {
    i = std::lower_bound(
            line.rawLine().begin() + i, line.rawLine().end(),
            XSortedTuple<T>{
                {boundedSub(boundedSub(p.getX(), line.getMaxSegLen()), padding),
                 0},
                false}) -
        line.rawLine().begin();
  }

  while (i < line.rawLine().size() &&
         line.rawLine()[i].seg().second.getX() < p.getX() - padding)
    i++;

  for (; i < line.rawLine().size(); i++) {
    const auto& cur = line.rawLine()[i];
    if (cur.out()) continue;
    if (cur.seg().first.getX() - padding > p.getX()) break;

    double euclideanDist = dist(cur.seg(), p);
    if (euclideanDist <= padding) {
      auto p2 = projectOn(cur.seg().first, p, cur.seg().second);
      double dist = distFunc(p, p2, std::min(maxDist, minDist));
      if (dist <= maxDist && dist < minDist) minDist = dist;
      if (minDist == 0) return 0;
    }

    if (euclideanDist < euclideanDistUpperBound) {
      euclideanDistUpperBound = euclideanDist;
      padding = paddingFunc(euclideanDistUpperBound, minDist, getBoundingBox(p),
                            line.boundingBox());
    }
  }

  return minDist;
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const XSortedLine<T>& line, const Point<T>& p, double maxDist,
                  PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  return withinDist(p, line, maxDist, paddingFunc, maxEuclideanDist, distFunc);
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool> intersectsContains(const Point<T>& p,
                                          const XSortedLine<T>& line,
                                          size_t i) {
  // check if point p lies on line
  // returns {intersects, contains}

  // skip irrelevant parts
  if (line.getMaxSegLen() < std::numeric_limits<T>::max() &&
      i < line.rawLine().size()) {
    i = std::lower_bound(
            line.rawLine().begin() + i, line.rawLine().end(),
            XSortedTuple<T>{{boundedSub(p.getX(), line.getMaxSegLen()), 0},
                            false}) -
        line.rawLine().begin();
  }

  for (; i < line.rawLine().size(); i++) {
    const auto& cur = line.rawLine()[i];
    if (cur.out()) continue;
    if (cur.seg().first.getX() > p.getX()) break;
    auto c = polyContCheck(p, cur.seg().first, cur.seg().second);
    if (c == 0)
      return {1, !((p == cur.origSeg().first && cur.rawPrevAng() == 32767) ||
                   (p == cur.origSeg().second && cur.rawNextAng() == 32767))};
  }

  return {0, 0};
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool> intersectsContains(const Point<T>& p,
                                          const XSortedLine<T>& line) {
  return intersectsContains(p, line, 0);
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const Point<T>& p, const XSortedLine<T>& line, size_t i) {
  auto res = intersectsContains(p, line, i);
  if (std::get<1>(res)) return M0FFFFF102;
  if (std::get<0>(res)) return MF0FFFF102;
  return MFF0FFF102;
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const Point<T>& p, const XSortedLine<T>& line) {
  return DE9IM(p, line, 0);
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const XSortedLine<T>& line, size_t i, const Point<T>& p) {
  auto res = intersectsContains(p, line, i);
  if (res.first) return M0F1FF0FF2;
  if (res.second) return MFF1F00FF2;
  return MFF1FF00F2;
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const XSortedLine<T>& line, const Point<T>& p) {
  return DE9IM(line, 0, p);
}

// _____________________________________________________________________________
template <typename T>
std::pair<bool, bool> containsCovers(const Point<T>& p,
                                     const XSortedPolygon<T>& poly, size_t i) {
  // returns {contains, covers}

  auto r = ringContains(p, poly.getOuter(), i);

  if (!r.second) return {0, 0};

  bool onInnerBorder = false;

  if (poly.getInners().size()) {
    size_t i = 0;

    if (poly.getInnerMaxSegLen() < std::numeric_limits<T>::max()) {
      i = std::lower_bound(
              poly.getInnerBoxIdx().begin(), poly.getInnerBoxIdx().end(),
              std::pair<T, size_t>{
                  boundedSub(p.getX(), poly.getInnerMaxSegLen()), 0}) -
          poly.getInnerBoxIdx().begin();
    }

    for (; i < poly.getInners().size(); i++) {
      if (poly.getInnerBoxes()[i].getLowerLeft().getX() > p.getX()) break;
      if (!util::geo::contains(p, poly.getInnerBoxes()[i])) continue;

      auto r = ringContains(p, poly.getInners()[i], 0);
      if (r.first) return {0, 0};  // surely contained in inner
      if (r.second) onInnerBorder = true;
    }
  }

  return {r.first && !onInnerBorder, true};
}

// _____________________________________________________________________________
template <typename T>
std::pair<bool, bool> containsCovers(const Point<T>& p,
                                     const XSortedPolygon<T>& poly) {
  return containsCovers(p, poly, 0);
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const Point<T>& p, const XSortedPolygon<T>& poly, size_t i) {
  auto res = containsCovers(p, poly, i);
  if (std::get<0>(res)) return M0FFFFF212;
  if (std::get<1>(res)) return MF0FFFF212;
  return MFF0FFF212;
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const Point<T>& p, const XSortedPolygon<T>& poly) {
  return DE9IM(p, poly, 0);
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const XSortedPolygon<T>& poly, size_t i, const Point<T>& p) {
  auto res = containsCovers(p, poly, i);
  if (res.first) return M0F2FF1FF2;
  if (res.second) return MFF20F1FF2;
  return MFF2FF10F2;
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const XSortedPolygon<T>& poly, const Point<T>& p) {
  return DE9IM(poly, 0, p);
}

// _____________________________________________________________________________
template <typename T>
bool ringContains(const Ring<T>& a, const Ring<T>& b) {
  for (size_t i = 1; i < a.size(); i++) {
    if (!ringContains(LineSegment<T>(a[i - 1], a[i]), b)) return false;
  }

  // also check the last hop
  if (!ringContains(LineSegment<T>(a.back(), a.front()), b)) return false;

  return true;
}

// _____________________________________________________________________________
template <typename T>
bool ringIntersects(const Polygon<T>& polyC, const Ring<T>& h) {
  // check if polygon polyC intersects the ring h

  // if h is contained in an inner ring of polyC, we do not intersect
  for (const auto& inner : polyC.getInners()) {
    if (ringContains(h, inner)) return false;
  }

  for (size_t i = 1; i < polyC.getOuter().size(); i++) {
    if (ringIntersects(
            LineSegment<T>(polyC.getOuter()[i - 1], polyC.getOuter()[i]), h))
      return true;
  }

  // also check the last hop
  if (ringIntersects(
          LineSegment<T>(polyC.getOuter().back(), polyC.getOuter().front()), h))
    return true;

  for (size_t i = 1; i < h.size(); i++) {
    if (ringIntersects(LineSegment<T>(h[i - 1], h[i]), polyC.getOuter()))
      return true;
  }

  // also check the last hop
  if (ringIntersects(LineSegment<T>(h.back(), h.front()), polyC.getOuter()))
    return true;

  return false;
}

// _____________________________________________________________________________
template <typename T>
bool ringContains(const Polygon<T>& polyC, const Ring<T>& poly) {
  return ringContains(polyC.getOuter(), poly);
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Polygon<T>& polyC, const Polygon<T>& poly) {
  // check if polygon polyC lies in polygon poly

  // check outer
  if (!ringContains(polyC, poly.getOuter())) return false;

  // check inners
  for (const auto& innerRing : poly.getInners()) {
    if (ringIntersects(polyC, innerRing)) return false;
  }

  return true;
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Polygon<T>& polyA, const Polygon<T>& polyB) {
  // check if polygon polyC intersects polygon poly

  // check outer
  if (!ringIntersects(polyA, polyB.getOuter()) &&
      !ringIntersects(polyB, polyA.getOuter()))
    return false;

  // check inners
  for (const auto& innerRing : polyB.getInners()) {
    if (ringContains(polyA, innerRing)) return false;
  }

  for (const auto& innerRing : polyA.getInners()) {
    if (ringContains(polyB, innerRing)) return false;
  }

  return true;
}

// _____________________________________________________________________________
template <typename T>
bool ringContains(const LineSegment<T>& ls, const Ring<T>& ph) {
  // check if linesegment ls lies in polygon poly

  // if one of the endpoints lies outside, abort
  if (!ringContains(ls.first, ph)) return false;
  if (!ringContains(ls.second, ph)) return false;

  for (size_t i = 1; i < ph.size(); i++) {
    auto seg = LineSegment<T>(ph[i - 1], ph[i]);
    if (!(contains(ls.first, seg) || contains(ls.second, seg)) &&
        intersects(seg, ls)) {
      return false;
    }
  }

  auto seg = LineSegment<T>(ph.back(), ph.front());
  if (!(contains(ls.first, seg) || contains(ls.second, seg)) &&
      intersects(seg, ls)) {
    return false;
  }

  return true;
}

// _____________________________________________________________________________
template <typename T>
bool contains(const LineSegment<T>& ls, const Polygon<T>& p) {
  if (!ringContains(ls, p.getOuter())) return false;

  for (const auto& inner : p.getInners()) {
    if (ringContains(ls, inner)) return false;
  }

  return true;
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Line<T>& l, const Polygon<T>& poly) {
  for (size_t i = 1; i < l.size(); i++) {
    if (!contains(LineSegment<T>(l[i - 1], l[i]), poly)) {
      return false;
    }
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Line<T>& l, const Line<T>& other) {
  for (const auto& p : l) {
    if (!contains(p, other)) return false;
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Box<T>& b, const Polygon<T>& poly) {
  return contains(convexHull(b), poly);
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Polygon<T>& poly, const Box<T>& b) {
  // check of poly lies in box
  for (const auto& p : poly.getOuter()) {
    if (!contains(p, b)) return false;
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
bool contains(const Polygon<T>& poly, const Line<T>& l) {
  for (const auto& p : poly.getOuter()) {
    if (!contains(p, l)) return false;
  }
  return true;
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T>
bool contains(const std::vector<GeometryA<T>>& multigeo,
              const GeometryB<T>& geo) {
  for (const auto& g : multigeo)
    if (!contains(g, geo)) return false;
  return true;
}

// _____________________________________________________________________________
template <typename T>
bool intersectsNaive(const std::vector<XSortedTuple<T>>& ls1,
                     const std::vector<XSortedTuple<T>>& ls2) {
  for (size_t i = 0; i < ls1.size(); i++) {
    if (ls1[i].out()) continue;
    for (size_t j = 0; j < ls2.size(); j++) {
      if (ls2[j].out()) continue;
      if (intersects(ls1[i].seg(), ls2[j].seg())) return true;
    }
  }

  return false;
}

// _____________________________________________________________________________
template <typename T, template <typename> class C>
uint8_t intersectsHelper(const std::vector<XSortedTuple<T>>& ls1,
                         const std::vector<XSortedTuple<T>>& ls2, T maxSegLenA,
                         T maxSegLenB, const Box<T>& boxA, const Box<T>& boxB,
                         size_t* firstRelIn1, size_t* firstRelIn2) {
  // returns {intersects, strict intersects, inside, overlaps}

  // ls2 is assumed to be a polygon ring
  if (ls1.size() == 0 || ls2.size() == 0) return 0;

  // shortcuts
  if (ls1.front().p.getX() > ls2.back().p.getX()) return 0;
  if (ls2.front().p.getX() > ls1.back().p.getX()) return 0;

  size_t i = 0;  // position in ls1
  size_t j = 0;  // position in ls2

  if (firstRelIn1) i = *firstRelIn1;
  if (firstRelIn2) j = *firstRelIn2;

  if (i >= ls1.size()) return 0;
  if (j >= ls2.size()) return 0;

  // skip irrelevant parts in ls1
  if (maxSegLenA < std::numeric_limits<T>::max() &&
      ls1[i].p.getX() < ls2[j].p.getX() - maxSegLenA) {
    i = std::lower_bound(
            ls1.begin() + i, ls1.end(),
            XSortedTuple<T>{{static_cast<T>(ls2[j].p.getX() - maxSegLenA), 0},
                            false}) -
        ls1.begin();
  }

  // skip irrelevant parts in ls2
  if (maxSegLenB < std::numeric_limits<T>::max() &&
      ls2[j].p.getX() < ls1[i].p.getX() - maxSegLenB) {
    j = std::lower_bound(
            ls2.begin() + j, ls2.end(),
            XSortedTuple<T>{{static_cast<T>(ls1[i].p.getX() - maxSegLenB), 0},
                            false}) -
        ls2.begin();
  }

  while (i < ls1.size() &&
         ls1[i].seg().second.getX() < boxB.getLowerLeft().getX())
    i++;
  while (j < ls2.size() &&
         ls2[j].seg().second.getX() < boxA.getLowerLeft().getX())
    j++;

  if (firstRelIn1) *firstRelIn1 = i;
  if (firstRelIn2) *firstRelIn2 = j;

  uint8_t ret = 0;

  std::set<AngledLineSegment<T>> active1, active2;

  while (i < ls1.size() && j < ls2.size()) {
    if (ls1[i].p.getX() < ls2[j].p.getX() ||
        (ls1[i].p.getX() == ls2[j].p.getX() &&
         (!ls1[i].out() || ls2[j].out()))) {
      // advance ls1

      // we are past ls2
      if (ls1[i].p.getX() > ls2.back().p.getX()) break;

      // ignore segments out of the X range
      if (ls1[i].seg().second.getX() < boxB.getLowerLeft().getX()) {
        i++;
        continue;
      }

      // ignore segments out of the Y range
      if (ls1[i].seg().first.getY() < boxB.getLowerLeft().getY() &&
          ls1[i].seg().second.getY() < boxB.getLowerLeft().getY()) {
        i++;
        continue;
      }

      if (ls1[i].seg().first.getY() > boxB.getUpperRight().getY() &&
          ls1[i].seg().second.getY() > boxB.getUpperRight().getY()) {
        i++;
        continue;
      }

      if (!ls1[i].out()) {
        auto above = active2.lower_bound(ls1[i].origSegAng());

        if (above != active2.end()) {
          auto a = above;

          do {
            auto r = C<T>::check(
                above->seg, above->prevAng, above->firstIsBoundary,
                above->nextAng, above->secondIsBoundary, ls1[i].origSeg(),
                ls1[i].rawPrevAng(), ls1[i].firstIsBoundary(),
                ls1[i].rawNextAng(), ls1[i].secondIsBoundary());
            ret |= r;

            if (!r) break;

            a = std::next(a);
          } while (a != active2.end());
        }
        if (above != active2.begin()) {
          auto aa = std::prev(above);
          auto a = aa;
          do {
            auto r =
                C<T>::check(a->seg, a->prevAng, a->firstIsBoundary, a->nextAng,
                            a->secondIsBoundary, ls1[i].origSeg(),
                            ls1[i].rawPrevAng(), ls1[i].firstIsBoundary(),
                            ls1[i].rawNextAng(), ls1[i].secondIsBoundary());
            ret |= r;

            if (!r || a == active2.begin()) break;
            a = std::prev(a);
          } while (1);
        }

        active1.insert(ls1[i].origSegAng());
      } else {
        auto toDel = active1.find(ls1[i].origSegAng());
        if (toDel != active1.end()) {
          auto above = active2.lower_bound(ls1[i].origSegAng());
          auto a = toDel;
          if (above != active2.end() && a != active1.begin()) {
            auto b = std::prev(a);
            do {
              bool found = false;
              auto aa = above;
              do {
                auto r = C<T>::check(aa->seg, aa->prevAng, aa->firstIsBoundary,
                                     aa->nextAng, aa->secondIsBoundary, b->seg,
                                     b->prevAng, b->firstIsBoundary, b->nextAng,
                                     b->secondIsBoundary);
                ret |= r;

                if (!r) break;
                found = true;

                aa = std::next(aa);
              } while (aa != active2.end());

              if (!found || b == active1.begin()) break;
              b = std::prev(b);
            } while (1);
          }

          toDel++;
          if (toDel != active1.end() && above != active2.begin()) {
            auto aa = std::prev(above);
            auto b = toDel;
            do {
              auto a = aa;
              bool found = false;
              do {
                auto r = C<T>::check(a->seg, a->prevAng, a->firstIsBoundary,
                                     a->nextAng, a->secondIsBoundary, b->seg,
                                     b->prevAng, b->firstIsBoundary, b->nextAng,
                                     b->secondIsBoundary);
                ret |= r;
                if (!r) break;
                found = true;
                if (a == active2.begin()) break;
                a = std::prev(a);
              } while (1);
              if (!found) break;
              b = std::next(b);
            } while (b != active1.end());
          }

          active1.erase(std::prev(toDel));
        }
      }

      i++;
    } else if (ls2[j].p.getX() < ls1[i].p.getX() ||
               (ls2[j].p.getX() == ls1[i].p.getX() && !ls2[j].out())) {
      // advance ls2

      // we are past ls1, so simply return (active1 is guaranteed to be
      // empty!)
      if (ls2[j].p.getX() > ls1.back().p.getX()) break;

      // ignore segments out of the X range
      if (ls2[j].seg().second.getX() < boxA.getLowerLeft().getX()) {
        j++;
        continue;
      }

      // ignore segments out of the Y range
      if (ls2[j].seg().first.getY() < boxA.getLowerLeft().getY() &&
          ls2[j].seg().second.getY() < boxA.getLowerLeft().getY()) {
        j++;
        continue;
      }
      if (ls2[j].seg().first.getY() > boxA.getUpperRight().getY() &&
          ls2[j].seg().second.getY() > boxA.getUpperRight().getY()) {
        j++;
        continue;
      }

      if (!ls2[j].out()) {
        auto above = active1.lower_bound(ls2[j].origSegAng());
        if (above != active1.end()) {
          auto a = above;
          do {
            auto r = C<T>::check(ls2[j].origSeg(), ls2[j].rawPrevAng(),
                                 ls2[j].firstIsBoundary(), ls2[j].rawNextAng(),
                                 ls2[j].secondIsBoundary(), a->seg, a->prevAng,
                                 a->firstIsBoundary, a->nextAng,
                                 a->secondIsBoundary);
            ret |= r;
            if (!r) break;
            a = std::next(a);
          } while (a != active1.end());
        }

        if (above != active1.begin()) {
          auto aa = std::prev(above);
          auto a = aa;
          do {
            auto r = C<T>::check(ls2[j].origSeg(), ls2[j].rawPrevAng(),
                                 ls2[j].firstIsBoundary(), ls2[j].rawNextAng(),
                                 ls2[j].secondIsBoundary(), a->seg, a->prevAng,
                                 a->firstIsBoundary, a->nextAng,
                                 a->secondIsBoundary);
            ret |= r;
            if (!r || a == active1.begin()) break;
            a = std::prev(a);
          } while (1);
        }

        active2.insert(ls2[j].origSegAng());
      } else {
        auto toDel = active2.find(ls2[j].origSegAng());
        if (toDel != active2.end()) {
          auto above = active1.lower_bound(ls2[j].origSegAng());
          auto a = toDel;
          if (above != active1.end() && a != active2.begin()) {
            auto b = std::prev(a);
            do {
              auto aa = above;
              bool found = false;
              do {
                auto r = C<T>::check(b->seg, b->prevAng, b->firstIsBoundary,
                                     b->nextAng, b->secondIsBoundary, aa->seg,
                                     aa->prevAng, aa->firstIsBoundary,
                                     aa->nextAng, aa->secondIsBoundary);
                ret |= r;
                if (!r) break;
                found = true;

                aa = std::next(aa);
              } while (aa != active1.end());

              if (!found || b == active2.begin()) break;
              b = std::prev(b);
            } while (1);
          }

          toDel++;
          if (toDel != active2.end() && above != active1.begin()) {
            auto aa = std::prev(above);
            auto b = toDel;
            do {
              auto a = aa;
              bool found = false;
              do {
                auto r = C<T>::check(b->seg, b->prevAng, b->firstIsBoundary,
                                     b->nextAng, b->secondIsBoundary, a->seg,
                                     a->prevAng, a->firstIsBoundary, a->nextAng,
                                     a->secondIsBoundary);
                ret |= r;
                if (!r) break;
                found = true;
                if (a == active1.begin()) break;
                a = std::prev(a);
              } while (1);
              if (!found) break;
              b = std::next(b);
            } while (b != active2.end());
          }

          active2.erase(std::prev(toDel));
        }
      }

      j++;
    }
  }

  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool, bool, bool, bool> intersectsPoly(
    const std::vector<XSortedTuple<T>>& ls1,
    const std::vector<XSortedTuple<T>>& ls2, T maxSegLenA, T maxSegLenB,
    const Box<T>& boxA, const Box<T>& boxB, size_t* firstRelIn1,
    size_t* firstRelIn2) {
  uint8_t ret = intersectsHelper<T, IntersectorPoly>(
      ls1, ls2, maxSegLenA, maxSegLenB, boxA, boxB, firstRelIn1, firstRelIn2);

  return {(ret >> 0) & 1, (ret >> 1) & 1, (ret >> 2) & 1, (ret >> 3) & 1,
          (ret >> 4) & 1, (ret >> 5) & 1, (ret >> 6) & 1, (ret >> 7) & 1};
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool> intersectsLinePoly(
    const std::vector<XSortedTuple<T>>& ls1,
    const std::vector<XSortedTuple<T>>& ls2, T maxSegLenA, T maxSegLenB,
    const Box<T>& boxA, const Box<T>& boxB, size_t* firstRelIn1,
    size_t* firstRelIn2) {
  uint8_t ret = intersectsHelper<T, IntersectorLinePoly>(
      ls1, ls2, maxSegLenA, maxSegLenB, boxA, boxB, firstRelIn1, firstRelIn2);

  return {(ret >> 0) & 1, (ret >> 1) & 1, (ret >> 2) & 1, (ret >> 3) & 1,
          (ret >> 4) & 1};
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::XSortedLine<T>& a,
                 const util::geo::XSortedLine<T>& b, size_t* firstRelIn1,
                 size_t* firstRelIn2) {
  uint8_t ret = intersectsHelper<T, IntersectorLine>(
      a.rawLine(), b.rawLine(), a.getMaxSegLen(), b.getMaxSegLen(),
      a.boundingBox(), b.boundingBox(), firstRelIn1, firstRelIn2);

  const bool weakIntersect = ret;

  const bool strictIntersect = (ret >> 0) & 1;
  const bool overlaps = (ret >> 1) & 1;
  const bool crosses = (ret >> 4) & 1;
  const bool strictIntersect2 = (ret >> 5) & 1;
  const bool bFirstInA = (ret >> 2) & 1;
  const bool bLastInA = (ret >> 3) & 1;
  const bool aFirstInB = (ret >> 6) & 1;
  const bool aLastInB = (ret >> 7) & 1;

  const bool aInB = !crosses && !strictIntersect && weakIntersect;
  const bool bInA = !crosses && !strictIntersect2 && weakIntersect;

  bool aFirstOnBoundary =
      a.firstPoint() == b.firstPoint() || a.firstPoint() == b.lastPoint();
  bool aLastOnBoundary =
      a.lastPoint() == b.firstPoint() || a.lastPoint() == b.lastPoint();

  bool bFirstOnBoundary =
      b.firstPoint() == a.firstPoint() || b.firstPoint() == a.lastPoint();
  bool bLastOnBoundary =
      b.lastPoint() == a.firstPoint() || b.lastPoint() == a.lastPoint();

  uint16_t ii = overlaps ? D1 : (crosses ? D0 : F);
  uint16_t ib = (bFirstInA || bLastInA) ? D0 : F;
  uint16_t ie = aInB ? F : D1;
  uint16_t bi =
      ((aFirstInB && !aFirstOnBoundary) || (aLastInB && !aLastOnBoundary)) ? D0
                                                                           : F;
  uint16_t bb =
      (a.firstPoint() == b.firstPoint() || a.lastPoint() == b.firstPoint() ||
       a.lastPoint() == b.lastPoint() || a.firstPoint() == b.lastPoint())
          ? D0
          : F;
  uint16_t be =
      ((aFirstInB || aFirstOnBoundary) && (aLastInB || aLastOnBoundary)) ? F
                                                                         : D0;
  uint16_t ei = bInA ? F : D1;
  uint16_t eb =
      ((bFirstInA || bFirstOnBoundary) && (bLastInA || bLastOnBoundary)) ? F
                                                                         : D0;

  return (ii << 0) | (ib << 2) | (ie << 4) | (bi << 6) | (bb << 8) |
         (be << 10) | (ei << 12) | (eb << 14);
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::XSortedLine<T>& a,
                 const util::geo::XSortedLine<T>& b) {
  return util::geo::DE9IM(a, b, 0, 0);
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool> intersectsLinePoly(
    const XSortedRing<T>& ls1, const XSortedRing<T>& ls2, T maxSegLenA,
    T maxSegLenB, size_t* firstRelIn1, size_t* firstRelIn2) {
  return intersectsLinePoly(ls1.rawRing(), ls2.rawRing(), maxSegLenA,
                            maxSegLenB, ls1.boundingBox(), ls2.boundingBox(),
                            firstRelIn1, firstRelIn2);
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool, bool, bool, bool> intersectsPoly(
    const XSortedRing<T>& ls1, const XSortedRing<T>& ls2, T maxSegLenA,
    T maxSegLenB, size_t* firstRelIn1, size_t* firstRelIn2) {
  return intersectsPoly(ls1.rawRing(), ls2.rawRing(), maxSegLenA, maxSegLenB,
                        ls1.boundingBox(), ls2.boundingBox(), firstRelIn1,
                        firstRelIn2);
}

// _____________________________________________________________________________
template <typename T>
uint8_t IntersectorLine<T>::check(const LineSegment<T>& ls1, int16_t prevLs1Ang,
                                  bool ls1FirstIsBoundary, int16_t nextLs1Ang,
                                  bool ls1SecondIsBoundary,
                                  const LineSegment<T>& ls2, int16_t prevLs2Ang,
                                  bool ls2FirstIsBoundary, int16_t nextLs2Ang,
                                  bool ls2SecondIsBoundary) {
  // {strictly intersects, overlaps, ls1.first inside ls2, ls1.last inside ls2,
  // crosses, strictly intersects ls2/ls1, ls2.first inside ls1, ls2.last inside
  // ls1}

  // trivial case: no intersect
  if (!intersects(getBoundingBox(ls1), getBoundingBox(ls2))) return 0;

  const bool ls2FirstInLs1 = contains(ls2.first, ls1);
  const bool ls2SecondInLs1 = contains(ls2.second, ls1);
  const bool ls1FirstInLs2 = contains(ls1.first, ls2);
  const bool ls1SecondInLs2 = contains(ls1.second, ls2);

  // ls1 and ls2 are equivalent
  if (ls2FirstInLs1 && ls2SecondInLs1 && ls1FirstInLs2 && ls1SecondInLs2) {
    return boolArrToInt8(
        {0, 1,
         ls1FirstIsBoundary &&
             ((ls1.first == ls2.first && !ls2FirstIsBoundary) ||
              (ls1.first == ls2.second && !ls2SecondIsBoundary)),
         ls1SecondIsBoundary &&
             ((ls1.second == ls2.first && !ls2FirstIsBoundary) ||
              (ls1.second == ls2.second && !ls2SecondIsBoundary)),
         0, 0,
         ls2FirstIsBoundary &&
             ((ls2.first == ls1.first && !ls1FirstIsBoundary) ||
              (ls2.first == ls1.second && !ls1SecondIsBoundary)),
         ls2SecondIsBoundary &&
             ((ls2.second == ls1.first && !ls1FirstIsBoundary) ||
              (ls2.second == ls1.second && !ls1SecondIsBoundary))});
  }

  // ls2 is completely in ls1
  if (ls2FirstInLs1 && ls2SecondInLs1) {
    int32_t ang1 = ((angBetween(ls2.first, ls2.second) / M_PI) * 32766);
    if (nextLs2Ang == 32767 &&
        (ls2.second == ls1.first || ls2.second == ls1.second))
      ang1 = 32767;

    int32_t ang2 = ((angBetween(ls2.second, ls2.first) / M_PI) * 32766);
    if (prevLs2Ang == 32767 &&
        (ls2.first == ls1.first || ls2.first == ls1.second))
      ang2 = 32767;

    return boolArrToInt8({0, 1, 0, 0, 0,
                          !(prevLs2Ang == ang2 && nextLs2Ang == ang1),
                          ls2FirstIsBoundary, ls2SecondIsBoundary});
  }

  // ls1 is completely in ls2
  if (ls1FirstInLs2 && ls1SecondInLs2) {
    int32_t ang1 = ((angBetween(ls1.first, ls1.second) / M_PI) * 32766);
    if (nextLs1Ang == 32767 &&
        (ls1.second == ls2.first || ls1.second == ls2.second))
      ang1 = 32767;

    int32_t ang2 = ((angBetween(ls1.second, ls1.first) / M_PI) * 32766);
    if (prevLs1Ang == 32767 &&
        (ls1.first == ls2.first || ls1.first == ls2.second))
      ang2 = 32767;

    return boolArrToInt8(
        {!(prevLs1Ang == ang2 && nextLs1Ang == ang1), 1,
         ls1FirstIsBoundary &&
             !(ls2FirstIsBoundary && ls1.first == ls2.first) &&
             !(ls2SecondIsBoundary && ls1.first == ls2.second),
         ls1SecondIsBoundary &&
             !(ls2FirstIsBoundary && ls1.second == ls2.first) &&
             !(ls2SecondIsBoundary && ls1.second == ls2.second),
         0, 0, 0, 0});
  }

  if (ls1.first == ls2.first && !ls1SecondInLs2 && !ls2SecondInLs1) {
    int32_t ang1 = ((angBetween(ls2.first, ls2.second) / M_PI) * 32766);
    int32_t ang2 = ((angBetween(ls1.first, ls1.second) / M_PI) * 32766);

    return boolArrToInt8(
        {prevLs1Ang != ang1, prevLs1Ang == ang1 || prevLs2Ang == ang2,
         ls1FirstIsBoundary && !ls2FirstIsBoundary, 0,
         prevLs2Ang != ang2 && prevLs1Ang != ang1 && prevLs1Ang != 32767 &&
             prevLs2Ang != 32767 && prevLs2Ang != prevLs1Ang,
         prevLs2Ang != ang2, ls2FirstIsBoundary && !ls1FirstIsBoundary, 0});
  }

  if (ls1.first == ls2.second && !ls1SecondInLs2 && !ls2FirstInLs1) {
    int32_t ang1 = ((angBetween(ls2.second, ls2.first) / M_PI) * 32766);
    int32_t ang2 = ((angBetween(ls1.first, ls1.second) / M_PI) * 32766);

    return boolArrToInt8(
        {prevLs1Ang != ang1, prevLs1Ang == ang1 || nextLs2Ang == ang2,
         ls1FirstIsBoundary && !ls2SecondIsBoundary, 0,
         nextLs2Ang != ang2 && prevLs1Ang != ang1 && prevLs1Ang != 32767 &&
             nextLs2Ang != 32767 && prevLs1Ang != nextLs2Ang,
         nextLs2Ang != ang2, 0, ls2SecondIsBoundary && !ls1FirstIsBoundary});
  }

  if (ls1.second == ls2.first && !ls1FirstInLs2 && !ls2SecondInLs1) {
    int32_t ang1 = ((angBetween(ls2.first, ls2.second) / M_PI) * 32766);
    int32_t ang2 = ((angBetween(ls1.second, ls1.first) / M_PI) * 32766);

    return boolArrToInt8(
        {nextLs1Ang != ang1, nextLs1Ang == ang1 || prevLs2Ang == ang2, 0,
         ls1SecondIsBoundary && !ls2FirstIsBoundary,
         prevLs2Ang != ang2 && nextLs1Ang != ang1 && nextLs1Ang != 32767 &&
             prevLs2Ang != 32767 && prevLs2Ang != nextLs1Ang,
         prevLs2Ang != ang2, ls2FirstIsBoundary && !ls1SecondIsBoundary, 0});
  }

  if (ls1.second == ls2.second && !ls1FirstInLs2 && !ls2FirstInLs1) {
    int32_t ang1 = ((angBetween(ls2.second, ls2.first) / M_PI) * 32766);
    int32_t ang2 = ((angBetween(ls1.second, ls1.first) / M_PI) * 32766);

    return boolArrToInt8(
        {ang1 != nextLs1Ang, ang1 == nextLs1Ang || ang2 == nextLs2Ang, 0,
         ls1SecondIsBoundary && !ls2SecondIsBoundary,
         nextLs2Ang != ang2 && nextLs1Ang != ang1 && nextLs1Ang != 32767 &&
             nextLs2Ang != 32767 && nextLs1Ang != nextLs2Ang,
         nextLs2Ang != ang2, 0, ls2SecondIsBoundary && !ls1SecondIsBoundary});
  }

  if (ls2FirstInLs1 && ls1SecondInLs2) {
    int32_t ang1 = ((angBetween(ls1.first, ls1.second) / M_PI) * 32766);
    int32_t ang2 = ((angBetween(ls2.second, ls2.first) / M_PI) * 32766);
    return boolArrToInt8({nextLs1Ang != ang1, 1, 0, ls1SecondIsBoundary, 0,
                          prevLs2Ang != ang2, ls2FirstIsBoundary, 0});
  }

  if (ls2SecondInLs1 && ls1SecondInLs2) {
    int32_t ang1 = ((angBetween(ls1.first, ls1.second) / M_PI) * 32766);
    int32_t ang2 = ((angBetween(ls2.first, ls2.second) / M_PI) * 32766);
    return boolArrToInt8({nextLs1Ang != ang1, 1, 0, ls1SecondIsBoundary, 0,
                          nextLs2Ang != ang2, 0, ls2SecondIsBoundary});
  }

  if (ls2FirstInLs1 && ls1FirstInLs2) {
    int32_t ang1 = ((angBetween(ls1.second, ls1.first) / M_PI) * 32766);
    int32_t ang2 = ((angBetween(ls2.second, ls2.first) / M_PI) * 32766);
    return boolArrToInt8({prevLs1Ang != ang1, 1, ls1FirstIsBoundary, 0, 0,
                          prevLs2Ang != ang2, ls2FirstIsBoundary, 0});
  }

  if (ls2SecondInLs1 && ls1FirstInLs2) {
    int16_t ang1 = ((angBetween(ls1.second, ls1.first) / M_PI) * 32766);
    int16_t ang2 = ((angBetween(ls2.first, ls2.second) / M_PI) * 32766);
    return boolArrToInt8({prevLs1Ang != ang1, 1, ls1FirstIsBoundary, 0, 0,
                          nextLs2Ang != ang2, 0, ls2SecondIsBoundary});
  }

  if (contains(ls1.second, ls2) && !ls1FirstInLs2 && !ls2FirstInLs1 &&
      !ls2SecondInLs1) {
    int16_t ang1 = (angBetween(ls2.first, ls2.second) / M_PI) * 32766;
    int16_t ang2 = (angBetween(ls2.second, ls2.first) / M_PI) * 32766;

    return boolArrToInt8(
        {1, nextLs1Ang == ang1 || nextLs1Ang == ang2, 0, ls1SecondIsBoundary,
         nextLs1Ang != 32767 && nextLs1Ang != ang1 && nextLs1Ang != ang2, 1, 0,
         0});
  }

  if (ls1FirstInLs2 && !ls1SecondInLs2 && !ls2FirstInLs1 && !ls2SecondInLs1) {
    int16_t ang1 = (angBetween(ls2.first, ls2.second) / M_PI) * 32766;
    int16_t ang2 = (angBetween(ls2.second, ls2.first) / M_PI) * 32766;

    return boolArrToInt8(
        {1, prevLs1Ang == ang1 || prevLs1Ang == ang2, ls1FirstIsBoundary, 0,
         prevLs1Ang != 32767 && prevLs1Ang != ang1 && prevLs1Ang != ang2, 1, 0,
         0});
  }

  if (ls2FirstInLs1 && !ls1SecondInLs2 && !ls1FirstInLs2 && !ls2SecondInLs1) {
    int16_t ang1 = (angBetween(ls1.first, ls1.second) / M_PI) * 32766;
    int16_t ang2 = (angBetween(ls1.second, ls1.first) / M_PI) * 32766;

    return boolArrToInt8(
        {1, prevLs2Ang == ang1 || prevLs2Ang == ang2, 0, 0,
         prevLs2Ang != 32767 && prevLs2Ang != ang1 && prevLs2Ang != ang2, 1,
         ls2FirstIsBoundary, 0});
  }

  if (ls2SecondInLs1 && !ls2FirstInLs1 && !ls1FirstInLs2 && !ls1SecondInLs2) {
    int16_t ang1 = (angBetween(ls1.first, ls1.second) / M_PI) * 32766;
    int16_t ang2 = (angBetween(ls1.second, ls1.first) / M_PI) * 32766;

    return boolArrToInt8(
        {1, nextLs2Ang == ang1 || nextLs2Ang == ang2, 0, 0,
         nextLs2Ang != 32767 && nextLs2Ang != ang1 && nextLs2Ang != ang2, 1, 0,
         ls2SecondIsBoundary});
  }

  // the line segments strictly intersect
  if (((crossProd(ls1.first, ls2) < 0) ^ (crossProd(ls1.second, ls2) < 0)) &&
      ((crossProd(ls2.first, ls1) < 0) ^ (crossProd(ls2.second, ls1) < 0))) {
    return bit16(00110001);
  }

  return 0;
}

// _____________________________________________________________________________
template <typename T>
uint8_t IntersectorLinePoly<T>::check(const LineSegment<T>& ls1,
                                      int16_t prevLs1Ang, bool,
                                      int16_t nextLs1Ang, bool,
                                      const LineSegment<T>& ls2,
                                      int16_t prevLs2Ang, bool,
                                      int16_t nextLs2Ang, bool) {
  // trivial case: no intersect
  if (!intersects(getBoundingBox(ls1), getBoundingBox(ls2))) return 0;

  int8_t ret =
      IntersectorPoly<T>::checkOneDir(ls1, prevLs1Ang, nextLs1Ang, ls2);

  if (prevLs2Ang == 32767 && contains(ls2.first, ls1)) ret |= bit16(10000);
  if (nextLs2Ang == 32767 && contains(ls2.second, ls1)) ret |= bit16(10000);

  return ret;
}

// _____________________________________________________________________________
template <typename T>
uint8_t IntersectorPoly<T>::checkOneDir(const LineSegment<T>& ls1,
                                        int16_t prevLs1Ang, int16_t nextLs1Ang,
                                        const LineSegment<T>& ls2) {
  const bool ls2FirstInLs1 = contains(ls2.first, ls1);
  const bool ls2SecondInLs1 = contains(ls2.second, ls1);

  // ls2 is completely in ls1
  if (ls2FirstInLs1 && ls2SecondInLs1) return bit16(1001);

  bool ls1FirstInLs2 = contains(ls1.first, ls2);
  bool ls1SecondInLs2 = contains(ls1.second, ls2);

  if (ls1FirstInLs2 && ls1SecondInLs2) {
    // this is not strictly correct, as we might overlap into the inside
    // of the polygon, but it doesn't matter, because this will be fetched
    // by the adjacent line segments of ls1
    return bit16(1001);
  }

  if (ls1.first == ls2.first && !ls1SecondInLs2 && !ls2SecondInLs1) {
    int16_t ang =
        (angBetween(
             ls1.second, ls1.first,
             {ls2.second.getX() - (ls1.first.getX() - ls1.second.getX()),
              ls2.second.getY() - (ls1.first.getY() - ls1.second.getY())}) /
         M_PI) *
        32766;
    if (ang > prevLs1Ang) return bit16(0101);
    if (ang == prevLs1Ang) return bit16(1001);

    return bit16(0011);
  }

  if (ls1.first == ls2.second && !ls1SecondInLs2 && !ls2FirstInLs1) {
    int16_t ang =
        (angBetween(
             ls1.second, ls1.first,
             {ls2.first.getX() - (ls1.first.getX() - ls1.second.getX()),
              ls2.first.getY() - (ls1.first.getY() - ls1.second.getY())}) /
         M_PI) *
        32766;
    if (ang > prevLs1Ang) return bit16(0101);
    if (ang == prevLs1Ang) return bit16(1001);
    return bit16(0011);
  }

  if (ls1.second == ls2.first && !ls1FirstInLs2 && !ls2SecondInLs1) {
    int16_t ang =
        (angBetween(
             ls1.first, ls1.second,
             {ls2.second.getX() - (ls1.second.getX() - ls1.first.getX()),
              ls2.second.getY() - (ls1.second.getY() - ls1.first.getY())}) /
         M_PI) *
        32766;

    if (ang < nextLs1Ang) return bit16(0101);
    if (ang == nextLs1Ang) return bit16(1001);
    return bit16(0011);
  }

  if (ls1.second == ls2.second && !ls1FirstInLs2 && !ls2FirstInLs1) {
    int16_t ang =
        (angBetween(
             ls1.first, ls1.second,
             {ls2.first.getX() - (ls1.second.getX() - ls1.first.getX()),
              ls2.first.getY() - (ls1.second.getY() - ls1.first.getY())}) /
         M_PI) *
        32766;

    if (ang < nextLs1Ang) return bit16(0101);
    if (ang == nextLs1Ang) return bit16(1001);
    return bit16(0011);
  }

  if (ls1FirstInLs2 && !ls1SecondInLs2 && !ls2FirstInLs1 && !ls2SecondInLs1) {
    // ls1.first is strictly (excluding end-points) on ls2
    int16_t ang1 =
        (angBetween(
             ls1.second, ls1.first,
             {ls2.first.getX() - (ls1.first.getX() - ls1.second.getX()),
              ls2.first.getY() - (ls1.first.getY() - ls1.second.getY())}) /
         M_PI) *
        32766;
    int16_t ang2 =
        (angBetween(
             ls1.second, ls1.first,
             {ls2.second.getX() - (ls1.first.getX() - ls1.second.getX()),
              ls2.second.getY() - (ls1.first.getY() - ls1.second.getY())}) /
         M_PI) *
        32766;

    if (ang1 > prevLs1Ang) return bit16(0111);
    if (ang2 > prevLs1Ang) return bit16(0111);
    return bit16(0011);
  }

  if (ls1SecondInLs2 && !ls1FirstInLs2 && !ls2FirstInLs1 && !ls2SecondInLs1) {
    // ls1.second is strictly (excluding end-points) on ls2
    int16_t ang1 =
        (angBetween(
             ls1.first, ls1.second,
             {ls2.first.getX() - (ls1.second.getX() - ls1.first.getX()),
              ls2.first.getY() - (ls1.second.getY() - ls1.first.getY())}) /
         M_PI) *
        32766;
    int16_t ang2 =
        (angBetween(
             ls1.first, ls1.second,
             {ls2.second.getX() - (ls1.second.getX() - ls1.first.getX()),
              ls2.second.getY() - (ls1.second.getY() - ls1.first.getY())}) /
         M_PI) *
        32766;
    if (ang1 < nextLs1Ang) return bit16(0111);
    if (ang2 < nextLs1Ang) return bit16(0111);
    return bit16(0011);
  }

  if (ls2FirstInLs1 && !ls1SecondInLs2 && !ls1FirstInLs2 && !ls2SecondInLs1) {
    // ls2.first is strictly (excluding end-points) on ls1
    if (crossProd(ls2.second, ls1) > 0) return bit16(0101);
    return bit16(0011);
  }

  if (ls2SecondInLs1 && !ls2FirstInLs1 && !ls1FirstInLs2 && !ls1SecondInLs2) {
    // ls2.second is strictly (excluding end-points) on ls1
    if (crossProd(ls2.first, ls1) > 0) return bit16(0101);
    return bit16(0011);
  }

  // the line segments strictly intersect
  if (((crossProd(ls1.first, ls2) < 0) ^ (crossProd(ls1.second, ls2) < 0)) &&
      ((crossProd(ls2.first, ls1) < 0) ^ (crossProd(ls2.second, ls1) < 0))) {
    return bit16(0111);
  }

  return 0;
}

// _____________________________________________________________________________
template <typename T>
uint8_t IntersectorPoly<T>::check(const LineSegment<T>& ls1, int16_t prevLs1Ang,
                                  bool, int16_t nextLs1Ang, bool,
                                  const LineSegment<T>& ls2, int16_t prevLs2Ang,
                                  bool, int16_t nextLs2Ang, bool) {
  // returns {intersects, strict intersects, inside, overlap}
  // returns a bit array, where bit 1 means "intersects", bit 2 means "strict
  // intersects" (that is the the segments intersects in a single point "into
  // the exterior") and bit 3 means "inside", which means that ls2 intersects
  // with the interior of the polygon of which ls1 is part of (but only
  // intersecting line segments are reported this way)

  // Examples: 001 means "intersects", 101 means "intersects" and "inside", but
  // not "strong intersect"

  // ls1 is a polygons line segment. we assume a clockwise ordering

  // trivial case: no intersect
  if (!intersects(getBoundingBox(ls1), getBoundingBox(ls2))) return 0;

  int8_t ret = 0;

  ret |= IntersectorPoly<T>::checkOneDir(ls1, prevLs1Ang, nextLs1Ang, ls2);
  ret |= (IntersectorPoly<T>::checkOneDir(ls2, prevLs2Ang, nextLs2Ang, ls1))
         << 4;

  return ret;
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const LineSegment<T>& ls1, const LineSegment<T>& ls2) {
  // check if two linesegments intersect

  // two line segments intersect of there is a single, well-defined intersection
  // point between them. If more than 1 endpoint is colinear with any line,
  // the segments have infinite intersections. We handle this case as non-
  // intersecting
  return intersects(getBoundingBox(ls1), getBoundingBox(ls2)) &&
         ((dist(ls1.first, ls2.first) < EPSILON && !contains(ls1.second, ls2) &&
           !contains(ls2.second, ls1)) ||
          (dist(ls1.second, ls2.first) < EPSILON && !contains(ls1.first, ls2) &&
           !contains(ls2.second, ls1)) ||
          (dist(ls1.second, ls2.second) < EPSILON &&
           !contains(ls1.first, ls2) && !contains(ls2.first, ls1)) ||
          (dist(ls1.first, ls2.first) < EPSILON && !contains(ls1.second, ls2) &&
           !contains(ls2.second, ls1)) ||

          (contains(ls1.first, ls2) && !contains(ls1.second, ls2) &&
           !contains(ls2.first, ls1) && !contains(ls2.second, ls1)) ||
          (contains(ls1.second, ls2) && !contains(ls1.first, ls2) &&
           !contains(ls2.first, ls1) && !contains(ls2.second, ls1)) ||
          (contains(ls2.first, ls1) && !contains(ls1.second, ls2) &&
           !contains(ls1.first, ls2) && !contains(ls2.second, ls1)) ||
          (contains(ls2.second, ls1) && !contains(ls2.first, ls1) &&
           !contains(ls1.first, ls2) && !contains(ls1.second, ls2)) ||

          (((crossProd(ls1.first, ls2) < 0) ^
            (crossProd(ls1.second, ls2) < 0)) &&
           ((crossProd(ls2.first, ls1) < 0) ^
            (crossProd(ls2.second, ls1) < 0))));
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Point<T>& a, const Point<T>& b, const Point<T>& c,
                const Point<T>& d) {
  // legacy function
  return intersects(LineSegment<T>(a, b), LineSegment<T>(c, d));
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Line<T>& ls1, const Line<T>& ls2) {
  for (size_t i = 1; i < ls1.size(); i++) {
    for (size_t j = 1; j < ls2.size(); j++) {
      if (intersects(LineSegment<T>(ls1[i - 1], ls1[i]),
                     LineSegment<T>(ls2[j - 1], ls2[j])))
        return true;
    }
  }

  return false;
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Line<T>& l, const Point<T>& p) {
  return contains(l, p);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Point<T>& p, const Line<T>& l) {
  return intersects(l, p);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Polygon<T>& l, const Point<T>& p) {
  return contains(p, l);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Point<T>& p, const Polygon<T>& l) {
  return intersects(l, p);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Box<T>& b1, const Box<T>& b2) {
  return b1.getLowerLeft().getX() <= b2.getUpperRight().getX() &&
         b1.getUpperRight().getX() >= b2.getLowerLeft().getX() &&
         b1.getLowerLeft().getY() <= b2.getUpperRight().getY() &&
         b1.getUpperRight().getY() >= b2.getLowerLeft().getY();
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Box<T>& b, const Polygon<T>& poly) {
  return intersects(poly, b);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Polygon<T>& poly, const Box<T>& b) {
  if (intersects(
          LineSegment<T>(b.getLowerLeft(), Point<T>(b.getUpperRight().getX(),
                                                    b.getLowerLeft().getY())),
          poly))
    return true;
  if (intersects(
          LineSegment<T>(b.getLowerLeft(), Point<T>(b.getLowerLeft().getX(),
                                                    b.getUpperRight().getY())),
          poly))
    return true;
  if (intersects(
          LineSegment<T>(b.getUpperRight(), Point<T>(b.getLowerLeft().getX(),
                                                     b.getUpperRight().getY())),
          poly))
    return true;
  if (intersects(
          LineSegment<T>(b.getUpperRight(), Point<T>(b.getUpperRight().getX(),
                                                     b.getLowerLeft().getY())),
          poly))
    return true;

  return contains(poly, b) || contains(b, poly);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const LineSegment<T>& ls, const Box<T>& b) {
  if (intersects(ls, LineSegment<T>(b.getLowerLeft(),
                                    Point<T>(b.getUpperRight().getX(),
                                             b.getLowerLeft().getY()))))
    return true;
  if (intersects(ls, LineSegment<T>(b.getLowerLeft(),
                                    Point<T>(b.getLowerLeft().getX(),
                                             b.getUpperRight().getY()))))
    return true;
  if (intersects(ls, LineSegment<T>(b.getUpperRight(),
                                    Point<T>(b.getLowerLeft().getX(),
                                             b.getUpperRight().getY()))))
    return true;
  if (intersects(ls, LineSegment<T>(b.getUpperRight(),
                                    Point<T>(b.getUpperRight().getX(),
                                             b.getLowerLeft().getY()))))
    return true;

  return contains(ls, b);
}

// _____________________________________________________________________________
template <typename T>
bool ringIntersects(const LineSegment<T>& ls, const Ring<T>& ph) {
  for (size_t i = 1; i < ph.size(); i++) {
    if (intersects(LineSegment<T>(ph[i - 1], ph[i]), ls)) return true;
  }

  // also check the last hop
  if (intersects(LineSegment<T>(ph.back(), ph.front()), ls)) return true;

  return ringContains(ls, ph);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const LineSegment<T>& ls, const Polygon<T>& p) {
  return ringIntersects(ls, p.getOuter());
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Polygon<T>& p, const LineSegment<T>& ls) {
  return intersects(ls, p);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Box<T>& b, const LineSegment<T>& ls) {
  return intersects(ls, b);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Line<T>& l, const Box<T>& b) {
  for (size_t i = 1; i < l.size(); i++) {
    if (intersects(LineSegment<T>(l[i - 1], l[i]), b)) return true;
  }
  return false;
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Box<T>& b, const Line<T>& l) {
  return intersects(l, b);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Point<T>& p, const Box<T>& b) {
  return contains(p, b);
}

// _____________________________________________________________________________
template <typename T>
bool intersects(const Box<T>& b, const Point<T>& p) {
  return intersects(p, b);
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T>
bool intersects(const std::vector<GeometryA<T>>& multigeom,
                const GeometryB<T>& b) {
  for (const auto& geom : multigeom)
    if (intersects(geom, b)) return true;
  return false;
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T>
bool intersects(const GeometryB<T>& b,
                const std::vector<GeometryA<T>>& multigeom) {
  return intersects(multigeom, b);
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T>
bool intersects(const std::vector<GeometryA<T>>& multigeomA,
                const std::vector<GeometryA<T>>& multigeomB) {
  for (const auto& geom : multigeomA)
    if (intersects(geom, multigeomB)) return true;
  return false;
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool, bool> intersectsContainsInner(
    const util::geo::XSortedLine<T>& a, const util::geo::XSortedRing<T>& b) {
  // returns {intersects, contains, covered, border intersects, border
  // intersects strict, border overlaps}

  size_t firstRel1 = 0;
  size_t firstRel2 = 0;

  if (a.rawLine().size() == 0 || b.rawRing().size() == 0)
    return {0, 0, 0, 0, 0, 0};

  auto borderInt = util::geo::intersectsPoly(
      a.rawLine(), b.rawRing(), a.getMaxSegLen(), b.getMaxSegLen(),
      util::geo::Box<T>(
          {std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()},
          {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()}),
      util::geo::Box<T>(
          {std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()},
          {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()}),
      &firstRel1, &firstRel2);

  if (std::get<1>(borderInt)) {
    // intersects
    return {true, false, !std::get<2>(borderInt),
            true, true,  std::get<3>(borderInt)};
  }

  if (util::geo::ringContains(a.rawLine().front().seg().second, b, firstRel2)
          .second) {
    return {true,
            !std::get<0>(borderInt),
            !std::get<2>(borderInt),
            std::get<0>(borderInt),
            std::get<1>(borderInt),
            std::get<3>(borderInt)};
  }

  // disjoint
  return {std::get<0>(borderInt),
          false,
          false,
          std::get<0>(borderInt),
          std::get<1>(borderInt),
          std::get<3>(borderInt)};
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool, bool, bool> intersectsContainsInner(
    const util::geo::XSortedRing<T>& a, const util::geo::XSortedRing<T>& b) {
  // returns {intersects, contains, covered, border isect, border isect strict,
  // border overlaps, border isect strict reverse}

  if (a.rawRing().size() == 0 || b.rawRing().size() == 0)
    return {0, 0, 0, 0, 0, 0, 0};

  size_t firstRel1 = 0;
  size_t firstRel2 = 0;
  auto borderInt = util::geo::intersectsPoly(
      a, b, a.getMaxSegLen(), b.getMaxSegLen(), &firstRel1, &firstRel2);

  if (std::get<1>(borderInt)) {
    // strictly intersects
    return {true,
            false,
            !std::get<2>(borderInt),
            true,
            true,
            std::get<3>(borderInt),
            std::get<5>(borderInt)};
  }

  if (util::geo::ringContains(a.rawRing().front().seg().second, b, firstRel2)
          .second) {
    return {true,
            !std::get<0>(borderInt),
            !std::get<2>(borderInt),
            std::get<0>(borderInt),
            false,
            std::get<3>(borderInt),
            std::get<5>(borderInt)};
  }

  if (util::geo::ringContains(b.rawRing().front().seg().second, a, firstRel1)
          .second) {
    // intersects (b is in a)
    return {true,
            false,
            false,
            std::get<0>(borderInt),
            false,
            std::get<3>(borderInt),
            std::get<5>(borderInt)};
  }

  // disjoint
  return {std::get<0>(borderInt),
          false,
          false,
          std::get<0>(borderInt),
          std::get<1>(borderInt),
          false,
          std::get<5>(borderInt)};
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::XSortedLine<T>& a,
                 const util::geo::XSortedPolygon<T>& b, size_t* firstRel1,
                 size_t* firstRel2) {
  if (a.rawLine().size() == 0) return MFF1FF0102;
  if (b.getOuter().rawRing().size() < 2) return MFF1FF0102;

  const auto boxA = a.boundingBox();
  const auto boxB = b.boundingBox();

  auto res = util::geo::intersectsLinePoly(
      a.rawLine(), b.getOuter().rawRing(), a.getMaxSegLen(),
      b.getOuter().getMaxSegLen(), boxA, boxB, firstRel1, firstRel2);

  uint16_t ib = std::get<3>(res) ? D1 : (std::get<2>(res) ? D0 : F);

  if (std::get<1>(res)) {
    uint16_t ii = std::get<2>(res) ? D1 : F;
    uint16_t ie = D1;

    uint16_t bi = F;
    uint16_t bb = F;

    if (std::get<4>(res)) {
      bi = F;
      bb = D0;
    } else {
      auto cc1 = containsCovers(a.firstPoint(), b, *firstRel2);
      auto cc2 = containsCovers(a.lastPoint(), b, *firstRel2);
      bi = std::get<0>(cc1) || std::get<0>(cc2) ? D0 : F;
      bb = (!std::get<0>(cc1) && std::get<1>(cc1)) ||
                   (!std::get<0>(cc2) && std::get<1>(cc2))
               ? D0
               : F;
    }

    uint16_t be = D0;
    uint16_t ei = D2;
    uint16_t eb = D1;

    return (ii << 0) | (ib << 2) | (ie << 4) | (bi << 6) | (bb << 8) |
           (be << 10) | (ei << 12) | (eb << 14);
  }

  if (util::geo::contains(boxA, boxB) &&
      (std::get<0>(res) ||
       util::geo::ringContains(a.rawLine().front().seg().second, b.getOuter(),
                               *firstRel2)
           .second)) {
    // a is contained or covered in the outer of B
    uint16_t ii = !std::get<0>(res) || std::get<2>(res) ? D1 : F;
    uint16_t ie = F;

    auto cc1 = containsCovers(a.firstPoint(), b, *firstRel2);
    auto cc2 = containsCovers(a.lastPoint(), b, *firstRel2);
    uint16_t bi = std::get<0>(cc1) || std::get<0>(cc2) ? D0 : F;
    uint16_t bb = (!std::get<0>(cc1) && std::get<1>(cc1)) ||
                          (!std::get<0>(cc2) && std::get<1>(cc2))
                      ? D0
                      : F;
    uint16_t be = (std::get<1>(cc1) && std::get<1>(cc2)) ? F : D0;
    uint16_t ei = D2;
    uint16_t eb = (!std::get<0>(res) || std::get<2>(res)) ? D1 : F;

    // check inner polygons
    for (size_t i = b.getFirstInner(a); i < b.getInners().size(); i++) {
      if (b.getInners()[i].rawRing().size() < 2) continue;
      const auto& innerBBox = b.getInnerBoxes()[i];
      if (b.getInnerBoxes()[i].getLowerLeft().getX() >
          a.rawLine().back().seg().second.getX())
        break;
      if (!util::geo::intersects(boxA, innerBBox)) continue;

      auto res = intersectsContainsInner(a, b.getInners()[i]);

      if (std::get<1>(res)) {
        return MFF1FF0102;  // a is completely contained by innerB
      }

      if (std::get<2>(res)) {
        // a is completely covered by inner b
        ii = F;
        ib = std::get<5>(res) ? D1 : F;
        ie = !std::get<5>(res) ? D1 : F;

        return (ii << 0) | (ib << 2) | (ie << 4) | (bi << 6) | (bb << 8) |
               (be << 10) | (ei << 12) | (eb << 14);
      }

      if (std::get<5>(res)) {
        ib = D1;
      }

      if (std::get<4>(res)) {
        // strictly intersects
        ii = D1;
        ie = D1;

        if (ib != D1) ib = std::get<5>(res) ? D1 : D0;
      }
    }

    return (ii << 0) | (ib << 2) | (ie << 4) | (bi << 6) | (bb << 8) |
           (be << 10) | (ei << 12) | (eb << 14);
  }

  return MFF1FF0102;
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::XSortedPolygon<T>& a,
                 const util::geo::XSortedLine<T>& b, size_t* firstRel1,
                 size_t* firstRel2) {
  DE9IMatrix de9im = DE9IM(b, a, firstRel2, firstRel1);
  return de9im.transpose();
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::XSortedLine<T>& a,
                 const util::geo::XSortedPolygon<T>& b) {
  size_t firstRel1 = 0;
  size_t firstRel2 = 0;
  return DE9IM(a, b, &firstRel1, &firstRel2);
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::XSortedPolygon<T>& a,
                 const util::geo::XSortedLine<T>& b) {
  size_t firstRel1 = 0;
  size_t firstRel2 = 0;
  return DE9IM(a, b, &firstRel1, &firstRel2);
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool> intersectsContainsCovers(
    const util::geo::XSortedPolygon<T>& a,
    const util::geo::XSortedPolygon<T>& b, size_t* firstRel1,
    size_t* firstRel2) {
  // returns {intersects, contains, covers, touches, overlaps}
  auto de9im = DE9IM(a, b, firstRel1, firstRel2);

  // see https://en.wikipedia.org/wiki/DE-9IM#Spatial_predicates
  return {de9im.intersects(), de9im.within(), de9im.covered(), de9im.touches(),
          de9im.overlaps02()};
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool> intersectsContainsCovers(
    const util::geo::XSortedPolygon<T>& a,
    const util::geo::XSortedPolygon<T>& b) {
  size_t firstA = 0;
  size_t firstB = 0;
  return intersectsContainsCovers(a, b, &firstA, &firstB);
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::XSortedPolygon<T>& a,
                 const util::geo::XSortedPolygon<T>& b, size_t* firstRel1,
                 size_t* firstRel2) {
  if (a.getOuter().rawRing().size() < 2) return MFF2FF1212;
  if (b.getOuter().rawRing().size() < 2) return MFF2FF1212;

  auto borderInt = util::geo::intersectsPoly(
      b.getOuter(), a.getOuter(), b.getOuter().getMaxSegLen(),
      a.getOuter().getMaxSegLen(), firstRel2, firstRel1);

  const Box<T>& boxA = a.boundingBox();
  const Box<T>& boxB = b.boundingBox();
  double outerAreaA = a.getOuter().area();
  double outerAreaB = b.getOuter().area();

  uint16_t bb = std::get<3>(borderInt) ? D1 : std::get<0>(borderInt) ? D0 : F;

  if (std::get<1>(borderInt) && std::get<5>(borderInt)) {
    uint16_t ii = (std::get<2>(borderInt) || std::get<6>(borderInt)) ? D2 : F;
    uint16_t ib = std::get<2>(borderInt) ? D1 : F;
    uint16_t ie = (!std::get<2>(borderInt) && std::get<6>(borderInt)) ? F : D2;
    uint16_t bi = (std::get<2>(borderInt) || std::get<6>(borderInt)) ? D1 : F;
    bb = std::get<3>(borderInt) ? D1 : D0;
    uint16_t be = (!std::get<2>(borderInt) && std::get<6>(borderInt)) ? F : D1;
    uint16_t ei = D2;
    uint16_t eb = D1;

    if (bb != D1) {
      // check inner polygons
      for (size_t i = b.getFirstInner(a); i < b.getInners().size(); i++) {
        const auto& innerBBox = b.getInnerBoxes()[i];
        if (innerBBox.getLowerLeft().getX() >
            a.getOuter().rawRing().back().seg().second.getX())
          break;

        if (!util::geo::intersects(innerBBox, boxA)) continue;
        const auto& innerB = b.getInners()[i];
        auto res = intersectsContainsInner(a.getOuter(), innerB);
        if (std::get<5>(res)) {
          bb = D1;
          break;
        }
      }
    }

    if (bb != D1) {
      for (size_t i = a.getFirstInner(b); i < a.getInners().size(); i++) {
        const auto& innerABox = a.getInnerBoxes()[i];
        if (innerABox.getLowerLeft().getX() >
            b.getOuter().rawRing().back().seg().second.getX())
          break;

        if (!util::geo::intersects(innerABox, boxB)) continue;

        const auto& innerA = a.getInners()[i];

        auto res = intersectsContainsInner(b.getOuter(), innerA);
        if (std::get<5>(res)) {
          bb = D1;
          break;
        }
      }
    }

    return (ii << 0) | (ib << 2) | (ie << 4) | (bi << 6) | (bb << 8) |
           (be << 10) | (ei << 12) | (eb << 14);
  }

  if (a.getOuter().rawRing().size() > 1 && outerAreaA <= outerAreaB &&
      util::geo::contains(boxA, boxB) &&
      (std::get<0>(borderInt) ||
       util::geo::ringContains(a.getOuter().rawRing().front().seg().second,
                               b.getOuter(), *firstRel2)
           .second)) {
    // the outer hull of A is inside the outer hull of B!
    bool equal = std::get<0>(borderInt) && !std::get<2>(borderInt) &&
                 !std::get<6>(borderInt);

    uint16_t ii = D2;
    uint16_t ib = F;
    uint16_t ie = F;
    uint16_t bi = equal ? F : D1;
    uint16_t be = F;
    uint16_t ei = equal ? F : D2;
    uint16_t eb = equal ? F : D1;

    // check inner polygons
    for (size_t i = b.getFirstInner(a); i < b.getInners().size(); i++) {
      const auto& innerBBox = b.getInnerBoxes()[i];
      if (innerBBox.getLowerLeft().getX() >
          a.getOuter().rawRing().back().seg().second.getX())
        break;

      if (!util::geo::intersects(innerBBox, boxA)) continue;

      const auto& innerB = b.getInners()[i];

      auto res = intersectsContainsInner(a.getOuter(), innerB);

      if (bb != D1) bb = std::get<5>(res) ? D1 : (std::get<3>(res) ? D0 : bb);

      // a ist completely contained in inner
      if (std::get<1>(res)) return MFF2FF1212;

      // a ist completely covered by inner
      if (std::get<2>(res)) {
        uint16_t be = outerAreaA < b.getInnerAreas()[i] ? D1 : F;
        return (F << 0) | (F << 2) | (D2 << 4) | (F << 6) | (bb << 8) |
               (be << 10) | (D2 << 12) | (D1 << 14);
      }

      // inner is completely covered by A
      if (std::get<0>(res) && !std::get<6>(res)) {
        // but inner may be covered/contained by an inner of A!
        bool contained = false;
        bool covered = false;
        bool borderOverlap = false;
        for (size_t i = a.getFirstInner(innerB); i < a.getInners().size();
             i++) {
          const auto& innerABox = a.getInnerBoxes()[i];
          if (innerABox.getLowerLeft().getX() >
              innerB.rawRing().back().seg().second.getX())
            break;

          if (!util::geo::intersects(innerABox, innerBBox)) continue;

          const auto& innerA = a.getInners()[i];

          auto res2 = intersectsContainsInner(innerB, innerA);

          if (std::get<2>(res2)) covered = true;
          if (std::get<1>(res2)) contained = true;
          if (std::get<5>(res2)) borderOverlap = true;
        }

        if (!covered) {
          ie = D2;
          ib = D1;
        }

        if (borderOverlap) bb = D1;
        if (bb == F && !contained && covered) bb = D0;
      }

      // inner and a truly overlap, but we must still check others for BB if
      // not yet maxed out
      if (!std::get<2>(res) && std::get<4>(res) && std::get<6>(res)) {
        ii = D2;
        ib = D1;
        ie = D2;
        bi = D1;
        be = D1;
        ei = D2;
        eb = D1;

        if (bb == D1) break;
      }
    }

    return (ii << 0) | (ib << 2) | (ie << 4) | (bi << 6) | (bb << 8) |
           (be << 10) | (ei << 12) | (eb << 14);
  }

  if (b.getOuter().rawRing().size() > 1 && outerAreaB <= outerAreaA &&
      util::geo::contains(boxB, boxA) &&
      (std::get<0>(borderInt) ||
       ringContains(b.getOuter().rawRing().front().seg().second, a.getOuter(),
                    *firstRel1)
           .second)) {
    // the outer hull of B is inside the outer hull of A!
    bool equal = std::get<0>(borderInt) && !std::get<2>(borderInt) &&
                 !std::get<6>(borderInt);

    uint16_t ii = D2;
    uint16_t ib = D1;
    uint16_t ie = equal ? F : D2;
    uint16_t bi = F;
    uint16_t be = equal ? F : D1;
    uint16_t ei = F;
    uint16_t eb = F;

    // check inner polygons
    for (size_t i = a.getFirstInner(b); i < a.getInners().size(); i++) {
      const auto& innerABox = a.getInnerBoxes()[i];
      if (innerABox.getLowerLeft().getX() >
          b.getOuter().rawRing().back().seg().second.getX())
        break;

      if (!util::geo::intersects(innerABox, boxB)) continue;

      const auto& innerA = a.getInners()[i];

      auto res = intersectsContainsInner(b.getOuter(), innerA);

      if (bb != D1) bb = std::get<5>(res) ? D1 : (std::get<3>(res) ? D0 : bb);

      // b ist completely contained in inner
      if (std::get<1>(res)) return MFF2FF1212;

      // b ist completely covered by inner
      if (std::get<2>(res)) {
        uint16_t eb = outerAreaB < a.getInnerAreas()[i] ? D1 : F;
        return (F << 0) | (F << 2) | (D2 << 4) | (F << 6) | (bb << 8) |
               (D1 << 10) | (D2 << 12) | (eb << 14);
      }

      // inner is completely covered by b
      if (std::get<0>(res) && !std::get<6>(res)) {
        // but inner may be covered/contained by an inner of B!
        bool contained = false;
        bool covered = false;
        bool borderOverlap = false;
        for (size_t i = b.getFirstInner(innerA); i < b.getInners().size();
             i++) {
          const auto& innerBBox = b.getInnerBoxes()[i];
          if (innerBBox.getLowerLeft().getX() >
              innerA.rawRing().back().seg().second.getX())
            break;

          if (!util::geo::intersects(innerBBox, innerABox)) continue;

          const auto& innerB = b.getInners()[i];

          auto res2 = intersectsContainsInner(innerA, innerB);

          if (std::get<2>(res2)) covered = true;
          if (std::get<1>(res2)) contained = true;
          if (std::get<5>(res2)) borderOverlap = true;
        }

        if (!covered) {
          ei = D2;
          bi = D1;
        }

        if (borderOverlap) bb = D1;
        if (bb == F && !contained && covered) bb = D0;
      }

      // inner and b truly overlap, but we must still check others for BB if
      // not yet maxed out
      if (!std::get<2>(res) && std::get<4>(res) && std::get<6>(res)) {
        ii = D2;
        ib = D1;
        ie = D2;
        bi = D1;
        be = D1;
        ei = D2;
        eb = D1;

        if (bb == D1) break;
      }
    }

    return (ii << 0) | (ib << 2) | (ie << 4) | (bi << 6) | (bb << 8) |
           (be << 10) | (ei << 12) | (eb << 14);
  }

  return MFF2FF1212;
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::Point<T>& a, const util::geo::Point<T>& b) {
  if (a == b) return M0FFFFFFF2;
  return MFF0FFF0F2;
}

// _____________________________________________________________________________
template <typename T>
DE9IMatrix DE9IM(const util::geo::XSortedPolygon<T>& a,
                 const util::geo::XSortedPolygon<T>& b) {
  size_t firstA = 0;
  size_t firstB = 0;
  return DE9IM(a, b, &firstA, &firstB);
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool> intersectsCovers(
    const util::geo::XSortedLine<T>& a, const util::geo::XSortedLine<T>& b,
    size_t* firstRelIn1, size_t* firstRelIn2) {
  auto de9im = DE9IM(a, b, firstRelIn1, firstRelIn2);

  // see https://en.wikipedia.org/wiki/DE-9IM#Spatial_predicates
  return {de9im.intersects(), de9im.covered(), de9im.touches(),
          de9im.overlaps1(), de9im.II() == D0};
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool> intersectsCovers(
    const util::geo::XSortedLine<T>& a, const util::geo::XSortedLine<T>& b) {
  return util::geo::intersectsCovers(a, b, 0, 0);
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool> intersectsContainsCovers(
    const util::geo::XSortedLine<T>& a, const util::geo::XSortedPolygon<T>& b,
    size_t* firstRel1, size_t* firstRel2) {
  // returns {intersects, contains, covers, touches, crosses}
  auto de9im = DE9IM(a, b, firstRel1, firstRel2);

  // see https://en.wikipedia.org/wiki/DE-9IM#Spatial_predicates
  return {de9im.intersects(), de9im.within(), de9im.covered(), de9im.touches(),
          de9im.II() && de9im.IE()};
}

// _____________________________________________________________________________
template <typename T>
std::tuple<bool, bool, bool, bool, bool> intersectsContainsCovers(
    const util::geo::XSortedLine<T>& a, const util::geo::XSortedPolygon<T>& b) {
  size_t firstRel1 = 0;
  size_t firstRel2 = 0;
  return intersectsContainsCovers(a, b, &firstRel1, &firstRel2);
}

// _____________________________________________________________________________
template <typename T>
Point<T> intersection(T p1x, T p1y, T q1x, T q1y, T p2x, T p2y, T q2x, T q2y) {
  /*
   * calculates the intersection between two line segments
   */
  if (doubleEq(p1x, q1x) && doubleEq(p1y, q1y))
    return Point<T>(p1x, p1y);  // TODO: <-- intersecting with a point??
  if (doubleEq(p2x, q1x) && doubleEq(p2y, q1y)) return Point<T>(p2x, p2y);
  if (doubleEq(p2x, q2x) && doubleEq(p2y, q2y))
    return Point<T>(p2x, p2y);  // TODO: <-- intersecting with a point??
  if (doubleEq(p1x, q2x) && doubleEq(p1y, q2y)) return Point<T>(p1x, p1y);

  double a = ((q2y - p2y) * (q1x - p1x)) - ((q2x - p2x) * (q1y - p1y));
  double u = (((q2x - p2x) * (p1y - p2y)) - ((q2y - p2y) * (p1x - p2x))) / a;

  return Point<T>(p1x + (q1x - p1x) * u, p1y + (q1y - p1y) * u);
}

// _____________________________________________________________________________
template <typename T>
Point<T> intersection(const Point<T>& p1, const Point<T>& q1,
                      const Point<T>& p2, const Point<T>& q2) {
  /*
   * calculates the intersection between two line segments
   */
  return intersection(p1.getX(), p1.getY(), q1.getX(), q1.getY(), p2.getX(),
                      p2.getY(), q2.getX(), q2.getY());
}

// _____________________________________________________________________________
template <typename T>
Point<T> intersection(const LineSegment<T>& s1, const LineSegment<T>& s2) {
  return intersection(s1.first, s1.second, s2.first, s2.second);
}

// _____________________________________________________________________________
template <typename T>
std::vector<Point<T>> intersection(const Line<T>& l1, const Line<T>& l2) {
  std::vector<Point<T>> ret;

  // TODO: better implementation than this naive baseline
  for (size_t i = 1; i < l1.size(); i++) {
    for (size_t j = 1; j < l2.size(); j++) {
      LineSegment<T> a = {l1[i - 1], l1[i]};
      LineSegment<T> b = {l2[j - 1], l2[j]};
      if (intersects(a, b)) ret.push_back(intersection(a, b));
    }
  }

  return ret;
}

// _____________________________________________________________________________
template <typename T>
Box<T> intersection(const Box<T>& b1, const Box<T>& b2) {
  if (!intersects(b1, b2)) return Box<T>();

  T llx, lly, urx, ury;

  if (b1.getLowerLeft().getX() > b2.getLowerLeft().getX())
    llx = b1.getLowerLeft().getX();
  else
    llx = b2.getLowerLeft().getX();

  if (b1.getLowerLeft().getY() > b2.getLowerLeft().getY())
    lly = b1.getLowerLeft().getY();
  else
    lly = b2.getLowerLeft().getY();

  if (b1.getUpperRight().getX() < b2.getUpperRight().getX())
    urx = b1.getUpperRight().getX();
  else
    urx = b2.getUpperRight().getX();

  if (b1.getUpperRight().getY() < b2.getUpperRight().getY())
    ury = b1.getUpperRight().getY();
  else
    ury = b2.getUpperRight().getY();

  return Box<T>{{llx, lly}, {urx, ury}};
}

// _____________________________________________________________________________
template <typename T>
bool lineIntersects(T p1x, T p1y, T q1x, T q1y, T p2x, T p2y, T q2x, T q2y) {
  /*
   * checks whether two lines intersect
   */
  double a = ((q2y - p2y) * (q1x - p1x)) - ((q2x - p2x) * (q1y - p1y));

  return a > EPSILON || a < -EPSILON;
}

// _____________________________________________________________________________
template <typename T>
bool lineIntersects(const Point<T>& p1, const Point<T>& q1, const Point<T>& p2,
                    const Point<T>& q2) {
  /*
   * checks whether two lines intersect
   */
  return lineIntersects(p1.getX(), p1.getY(), q1.getX(), q1.getY(), p2.getX(),
                        p2.getY(), q2.getX(), q2.getY());
}

// _____________________________________________________________________________
template <typename T>
double angBetween(const Point<T>& p1) {
  return atan2(p1.getX(), p1.getY());
}

// _____________________________________________________________________________
template <typename T>
double angBetween(const Point<T>& p1, const MultiPoint<T>& points) {
  double sinSum = 0;
  double cosSum = 0;
  for (auto q1 : points) {
    double a = angBetween(p1.getX(), p1.getY(), q1.getX(), q1.getY());
    sinSum += sin(a);
    cosSum += cos(a);
  }
  return atan2(sinSum / points.size(), cosSum / points.size());
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const AnyGeometry<T>& any1, const AnyGeometry<T>& any2) {
  if (any1.getType() == 0) return dist(any1.getPoint(), any2);
  if (any1.getType() == 1) return dist(any1.getLine(), any2);
  if (any1.getType() == 2) return dist(any1.getPolygon(), any2);
  if (any1.getType() == 3) return dist(any1.getMultiLine(), any2);
  if (any1.getType() == 4) return dist(any1.getMultiPolygon(), any2);
  if (any1.getType() == 5) return dist(any1.getCollection(), any2);
  if (any1.getType() == 6) return dist(any1.getMultiPoint(), any2);

  return std::numeric_limits<double>::infinity();
}

// _____________________________________________________________________________
template <template <typename> class GeometryB, typename T, typename DF>
double dist(const AnyGeometry<T>& any, const GeometryB<T>& geomB) {
  if (any.getType() == 0) return dist(any.getPoint(), geomB);
  if (any.getType() == 1) return dist(any.getLine(), geomB);
  if (any.getType() == 2) return dist(any.getPolygon(), geomB);
  if (any.getType() == 3) return dist(any.getMultiLine(), geomB);
  if (any.getType() == 4) return dist(any.getMultiPolygon(), geomB);
  if (any.getType() == 5) return dist(any.getCollection(), geomB);
  if (any.getType() == 6) return dist(any.getMultiPoint(), geomB);

  return std::numeric_limits<double>::infinity();
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const AnyGeometry<T>& any1, const AnyGeometry<T>& any2,
            DF&& distFunc) {
  if (any1.getType() == 0) return dist(any1.getPoint(), any2, distFunc);
  if (any1.getType() == 1) return dist(any1.getLine(), any2, distFunc);
  if (any1.getType() == 2) return dist(any1.getPolygon(), any2, distFunc);
  if (any1.getType() == 3) return dist(any1.getMultiLine(), any2, distFunc);
  if (any1.getType() == 4) return dist(any1.getMultiPolygon(), any2, distFunc);
  if (any1.getType() == 5) return dist(any1.getCollection(), any2, distFunc);
  if (any1.getType() == 6) return dist(any1.getMultiPoint(), any2, distFunc);

  return std::numeric_limits<double>::infinity();
}

// _____________________________________________________________________________
template <template <typename> class GeometryB, typename T, typename DF>
double dist(const AnyGeometry<T>& any, const GeometryB<T>& geomB,
            DF&& distFunc) {
  if (any.getType() == 0) return dist(any.getPoint(), geomB, distFunc);
  if (any.getType() == 1) return dist(any.getLine(), geomB, distFunc);
  if (any.getType() == 2) return dist(any.getPolygon(), geomB, distFunc);
  if (any.getType() == 3) return dist(any.getMultiLine(), geomB, distFunc);
  if (any.getType() == 4) return dist(any.getMultiPolygon(), geomB, distFunc);
  if (any.getType() == 5) return dist(any.getCollection(), geomB, distFunc);
  if (any.getType() == 6) return dist(any.getMultiPoint(), geomB, distFunc);

  return std::numeric_limits<double>::infinity();
}

// _____________________________________________________________________________
template <template <typename> class GeometryB, typename T, typename DF>
double dist(const GeometryB<T>& geomB, const AnyGeometry<T>& any,
            DF&& distFunc) {
  return dist(any, geomB, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const LineSegment<T>& ls, const Point<T>& p, DF&& distFunc) {
  return distToSegment(ls, p, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Point<T>& p, const LineSegment<T>& ls, DF&& distFunc) {
  return dist(ls, p, distFunc);
}

// _____________________________________________________________________________
template <typename T>
double distSquared(const LineSegment<T>& ls1, const LineSegment<T>& ls2) {
  if (intersects(ls1, ls2)) return 0;
  double d1 = distToSegmentSquared(ls2.first.getX(), ls2.first.getY(),
                                   ls2.second.getX(), ls2.second.getY(),
                                   ls1.first.getX(), ls1.first.getY());
  double d2 = distToSegmentSquared(ls2.first.getX(), ls2.first.getY(),
                                   ls2.second.getX(), ls2.second.getY(),
                                   ls1.second.getX(), ls1.second.getY());
  double d3 = distToSegmentSquared(ls1.first.getX(), ls1.first.getY(),
                                   ls1.second.getX(), ls1.second.getY(),
                                   ls2.first.getX(), ls2.first.getY());
  double d4 = distToSegmentSquared(ls1.first.getX(), ls1.first.getY(),
                                   ls1.second.getX(), ls1.second.getY(),
                                   ls2.second.getX(), ls2.second.getY());
  return (std::min(d1, std::min(d2, (std::min(d3, d4)))));
}

// _____________________________________________________________________________
template <typename T>
double dist(const LineSegment<T>& ls1, const LineSegment<T>& ls2) {
  if (intersects(ls1, ls2)) return 0;
  double d1 = distToSegmentSquared(ls2.first.getX(), ls2.first.getY(),
                                   ls2.second.getX(), ls2.second.getY(),
                                   ls1.first.getX(), ls1.first.getY());
  double d2 = distToSegmentSquared(ls2.first.getX(), ls2.first.getY(),
                                   ls2.second.getX(), ls2.second.getY(),
                                   ls1.second.getX(), ls1.second.getY());
  double d3 = distToSegmentSquared(ls1.first.getX(), ls1.first.getY(),
                                   ls1.second.getX(), ls1.second.getY(),
                                   ls2.first.getX(), ls2.first.getY());
  double d4 = distToSegmentSquared(ls1.first.getX(), ls1.first.getY(),
                                   ls1.second.getX(), ls1.second.getY(),
                                   ls2.second.getX(), ls2.second.getY());
  return sqrt(std::min(d1, std::min(d2, (std::min(d3, d4)))));
}

// _____________________________________________________________________________
template <typename T, typename DF>
double withinDist(const LineSegment<T>& ls1, const LineSegment<T>& ls2,
                  DF&& distFunc, double maxD) {
  if (intersects(ls1, ls2)) return 0;
  double d1 = distToSegmentSquared(ls2.first.getX(), ls2.first.getY(),
                                   ls2.second.getX(), ls2.second.getY(),
                                   ls1.first.getX(), ls1.first.getY());
  double d2 = distToSegmentSquared(ls2.first.getX(), ls2.first.getY(),
                                   ls2.second.getX(), ls2.second.getY(),
                                   ls1.second.getX(), ls1.second.getY());
  double d3 = distToSegmentSquared(ls1.first.getX(), ls1.first.getY(),
                                   ls1.second.getX(), ls1.second.getY(),
                                   ls2.first.getX(), ls2.first.getY());
  double d4 = distToSegmentSquared(ls1.first.getX(), ls1.first.getY(),
                                   ls1.second.getX(), ls1.second.getY(),
                                   ls2.second.getX(), ls2.second.getY());

  if (d1 <= d2 && d1 <= d3 && d1 <= d4) {
    auto p2 = projectOn(ls2.first, ls1.first, ls2.second);
    return distFunc(ls1.first, p2, maxD);
  }

  if (d2 <= d1 && d2 <= d3 && d2 <= d4) {
    auto p2 = projectOn(ls2.first, ls1.second, ls2.second);
    return distFunc(ls1.second, p2, maxD);
  }

  if (d3 <= d1 && d3 <= d2 && d3 <= d4) {
    auto p2 = projectOn(ls1.first, ls2.first, ls1.second);
    return distFunc(ls2.first, p2, maxD);
  }

  auto p2 = projectOn(ls1.first, ls2.second, ls1.second);
  return distFunc(ls2.second, p2, maxD);
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const LineSegment<T>& ls1, const LineSegment<T>& ls2,
            double padding, DF&& distFunc, double maxD, double& euclideanDist) {
  if (intersects(ls1, ls2)) {
    euclideanDist = 0;
    return 0;
  }

  double d1 = distToSegmentSquared(ls2.first.getX(), ls2.first.getY(),
                                   ls2.second.getX(), ls2.second.getY(),
                                   ls1.first.getX(), ls1.first.getY());
  double d2 = distToSegmentSquared(ls2.first.getX(), ls2.first.getY(),
                                   ls2.second.getX(), ls2.second.getY(),
                                   ls1.second.getX(), ls1.second.getY());
  double d3 = distToSegmentSquared(ls1.first.getX(), ls1.first.getY(),
                                   ls1.second.getX(), ls1.second.getY(),
                                   ls2.first.getX(), ls2.first.getY());
  double d4 = distToSegmentSquared(ls1.first.getX(), ls1.first.getY(),
                                   ls1.second.getX(), ls1.second.getY(),
                                   ls2.second.getX(), ls2.second.getY());

  double minSq = std::min(d1, std::min(d2, std::min(d3, d4)));
  euclideanDist = sqrt(minSq);

  // skip costly computation entirely
  if (euclideanDist > padding) return std::numeric_limits<double>::max();

  if (d1 <= d2 && d1 <= d3 && d1 <= d4) {
    auto p2 = projectOn(ls2.first, ls1.first, ls2.second);
    return distFunc(ls1.first, p2, maxD);
  }

  if (d2 <= d1 && d2 <= d3 && d2 <= d4) {
    auto p2 = projectOn(ls2.first, ls1.second, ls2.second);
    return distFunc(ls1.second, p2, maxD);
  }

  if (d3 <= d1 && d3 <= d2 && d3 <= d4) {
    auto p2 = projectOn(ls1.first, ls2.first, ls1.second);
    return distFunc(ls2.first, p2, maxD);
  }

  auto p2 = projectOn(ls1.first, ls2.second, ls1.second);
  return distFunc(ls2.second, p2, maxD);
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const LineSegment<T>& ls1, const LineSegment<T>& ls2,
            DF&& distFunc) {
  return withinDist(ls1, ls2, distFunc, std::numeric_limits<double>::max());
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Point<T>& p, const Line<T>& l, DF&& distFunc) {
  return withinDist(p, l, std::numeric_limits<double>::infinity(),
                    defaultPaddingFunc<T>(),
                    std::numeric_limits<double>::infinity(), distFunc);
}

// _____________________________________________________________________________
template <typename T>
double dist(const Point<T>& p, const Line<T>& l) {
  return withinDist(p, l, std::numeric_limits<double>::max());
}

// _____________________________________________________________________________
template <typename T>
double withinDist(const Point<T>& p, const Line<T>& l, double maxD) {
  if (l.size() > 1000) {
    return withinDist(p, XSortedLine<T>(l), maxD);
  }

  double d = maxD;
  for (size_t i = 1; i < l.size(); i++) {
    double dTmp = distToSegment(l[i - 1], l[i], p);
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const Point<T>& p, const Line<T>& l, double maxDist,
                  PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  if (l.size() > 1000) {
    return withinDist(p, XSortedLine<T>(l), maxDist, paddingFunc,
                      maxEuclideanDist, distFunc);
  }

  double d = maxDist;
  for (size_t i = 1; i < l.size(); i++) {
    double dTmp = distToSegment(l[i - 1], l[i], p, distFunc);
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Line<T>& l, const Point<T>& p, DF&& distFunc) {
  return dist(p, l, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const LineSegment<T>& ls, const Line<T>& l, DF&& distFunc) {
  double d = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < l.size(); i++) {
    double dTmp = dist(ls, LineSegment<T>(l[i - 1], l[i]), distFunc);
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
}

// _____________________________________________________________________________
template <typename T>
double dist(const LineSegment<T>& ls, const Line<T>& l) {
  double d = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < l.size(); i++) {
    double dTmp = dist(ls, LineSegment<T>(l[i - 1], l[i]));
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Line<T>& l, const LineSegment<T>& ls, DF&& distFunc) {
  return dist(ls, l, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Line<T>& la, const Line<T>& lb, DF&& distFunc) {
  double d = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < la.size(); i++) {
    double dTmp = dist(LineSegment<T>(la[i - 1], la[i]), lb, distFunc);
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
}

// _____________________________________________________________________________
template <typename T>
double dist(const Line<T>& la, const Line<T>& lb) {
  if (la.size() * lb.size() > EST_CHECKS_THRESHOLD_XSORTED) {
    return dist(XSortedLine<T>(la), XSortedLine<T>(lb));
  }
  double d = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < la.size(); i++) {
    double dTmp = dist(LineSegment<T>(la[i - 1], la[i]), lb);
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double withinDist(const LineSegment<T>& a, const Box<T>& b, DF&& distF,
                  double dMax) {
  if (util::geo::intersects(a, b)) return 0;

  double d = std::numeric_limits<double>::infinity();
  d = std::min(d, util::geo::withinDist(
                      a, LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()},
                      distF, dMax));
  d = std::min(d, util::geo::withinDist(
                      a, LineSegment<T>{b.getUpperLeft(), b.getUpperRight()},
                      distF, dMax));
  d = std::min(d, util::geo::withinDist(
                      a, LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                      distF, dMax));
  d = std::min(d, util::geo::withinDist(
                      a, LineSegment<T>{b.getLowerRight(), b.getLowerLeft()},
                      distF, dMax));
  return d;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double withinDist(const Box<T>& a, const Box<T>& b, DF&& distF, double dMax) {
  if (util::geo::intersects(a, b)) return 0;

  double d = std::numeric_limits<double>::infinity();
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getLowerLeft(), a.getUpperLeft()},
             LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()}, distF, dMax));
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getLowerLeft(), a.getUpperLeft()},
             LineSegment<T>{b.getUpperLeft(), b.getUpperRight()}, distF, dMax));
  d = std::min(d, util::geo::withinDist(
                      LineSegment<T>{a.getLowerLeft(), a.getUpperLeft()},
                      LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                      distF, dMax));
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getLowerLeft(), a.getUpperLeft()},
             LineSegment<T>{b.getLowerRight(), b.getLowerLeft()}, distF, dMax));

  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getUpperLeft(), a.getUpperRight()},
             LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()}, distF, dMax));
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getUpperLeft(), a.getUpperRight()},
             LineSegment<T>{b.getUpperLeft(), b.getUpperRight()}, distF, dMax));
  d = std::min(d, util::geo::withinDist(
                      LineSegment<T>{a.getUpperLeft(), a.getUpperRight()},
                      LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                      distF, dMax));
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getUpperLeft(), a.getUpperRight()},
             LineSegment<T>{b.getLowerRight(), b.getLowerLeft()}, distF, dMax));

  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getUpperRight(), a.getLowerRight()},
             LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()}, distF, dMax));
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getUpperRight(), a.getLowerRight()},
             LineSegment<T>{b.getUpperLeft(), b.getUpperRight()}, distF, dMax));
  d = std::min(d, util::geo::withinDist(
                      LineSegment<T>{a.getUpperRight(), a.getLowerRight()},
                      LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                      distF, dMax));
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getUpperRight(), a.getLowerRight()},
             LineSegment<T>{b.getLowerRight(), b.getLowerLeft()}, distF, dMax));

  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getLowerRight(), a.getLowerLeft()},
             LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()}, distF, dMax));
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getLowerRight(), a.getLowerLeft()},
             LineSegment<T>{b.getUpperLeft(), b.getUpperRight()}, distF, dMax));
  d = std::min(d, util::geo::withinDist(
                      LineSegment<T>{a.getLowerRight(), a.getLowerLeft()},
                      LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                      distF, dMax));
  d = std::min(
      d, util::geo::withinDist(
             LineSegment<T>{a.getLowerRight(), a.getLowerLeft()},
             LineSegment<T>{b.getLowerRight(), b.getLowerLeft()}, distF, dMax));

  return d;
}

// _____________________________________________________________________________
template <typename T>
double distSquared(const LineSegment<T>& a, const Box<T>& b) {
  double d = std::numeric_limits<double>::infinity();
  d = std::min(d, util::geo::distSquared(
                      a, LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()}));
  d = std::min(d, util::geo::distSquared(
                      a, LineSegment<T>{b.getUpperLeft(), b.getUpperRight()}));
  d = std::min(d, util::geo::distSquared(
                      a, LineSegment<T>{b.getUpperRight(), b.getLowerRight()}));
  d = std::min(d, util::geo::distSquared(
                      a, LineSegment<T>{b.getLowerRight(), b.getLowerLeft()}));

  return d;
}

// _____________________________________________________________________________
template <typename T>
double dist(const LineSegment<T>& a, const Box<T>& b) {
  double d = std::numeric_limits<double>::infinity();
  d = std::min(d, util::geo::dist(
                      a, LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()}));
  d = std::min(d, util::geo::dist(
                      a, LineSegment<T>{b.getUpperLeft(), b.getUpperRight()}));
  d = std::min(d, util::geo::dist(
                      a, LineSegment<T>{b.getUpperRight(), b.getLowerRight()}));
  d = std::min(d, util::geo::dist(
                      a, LineSegment<T>{b.getLowerRight(), b.getLowerLeft()}));

  return d;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Box<T>& a, const Box<T>& b, DF&& distFunc) {
  if (util::geo::intersects(a, b)) return 0;

  double d = std::numeric_limits<double>::infinity();
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getLowerLeft(), a.getUpperLeft()},
                         LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getLowerLeft(), a.getUpperLeft()},
                         LineSegment<T>{b.getUpperLeft(), b.getUpperRight()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getLowerLeft(), a.getUpperLeft()},
                         LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getLowerLeft(), a.getUpperLeft()},
                         LineSegment<T>{b.getLowerRight(), b.getLowerLeft()},
                         distFunc));

  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getUpperLeft(), a.getUpperRight()},
                         LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getUpperLeft(), a.getUpperRight()},
                         LineSegment<T>{b.getUpperLeft(), b.getUpperRight()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getUpperLeft(), a.getUpperRight()},
                         LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getUpperLeft(), a.getUpperRight()},
                         LineSegment<T>{b.getLowerRight(), b.getLowerLeft()},
                         distFunc));

  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getUpperRight(), a.getLowerRight()},
                         LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getUpperRight(), a.getLowerRight()},
                         LineSegment<T>{b.getUpperLeft(), b.getUpperRight()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getUpperRight(), a.getLowerRight()},
                         LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getUpperRight(), a.getLowerRight()},
                         LineSegment<T>{b.getLowerRight(), b.getLowerLeft()},
                         distFunc));

  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getLowerRight(), a.getLowerLeft()},
                         LineSegment<T>{b.getLowerLeft(), b.getUpperLeft()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getLowerRight(), a.getLowerLeft()},
                         LineSegment<T>{b.getUpperLeft(), b.getUpperRight()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getLowerRight(), a.getLowerLeft()},
                         LineSegment<T>{b.getUpperRight(), b.getLowerRight()},
                         distFunc));
  d = std::min(
      d, util::geo::dist(LineSegment<T>{a.getLowerRight(), a.getLowerLeft()},
                         LineSegment<T>{b.getLowerRight(), b.getLowerLeft()},
                         distFunc));

  return d;
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T, typename DF>
double dist(const std::vector<GeometryA<T>>& a, const GeometryB<T>& b,
            DF&& distFunc) {
  double d = std::numeric_limits<double>::infinity();
  for (const auto& geom : a) d = std::min(d, dist(geom, b, distFunc));
  return d;
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T, typename DF>
double dist(const GeometryB<T>& b, const std::vector<GeometryA<T>>& multigeom,
            DF&& distFunc) {
  return dist(multigeom, b, distFunc);
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T, typename DF>
double dist(const std::vector<GeometryA<T>>& a,
            const std::vector<GeometryB<T>>& b, DF&& distFunc) {
  double d = std::numeric_limits<double>::infinity();
  for (const auto& geom : b) d = std::min(d, dist(geom, a, distFunc));
  return d;
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T, typename DF,
          typename PF>
double dist(const std::vector<GeometryA<T>>& a,
            const std::vector<GeometryB<T>>& b, PF&& paddingFunc,
            DF&& distFunc) {
  if (a.size() * b.size() > EST_MULTI_CHECKS_THRESHOLD_XSORTED)
    return withinDist(XSortedCollection<T>(a), XSortedCollection<T>(b),
                      std::numeric_limits<double>::infinity(), paddingFunc,
                      std::numeric_limits<double>::infinity(), distFunc);
  double d = std::numeric_limits<double>::infinity();
  for (const auto& geom : b) d = std::min(d, dist(geom, a, distFunc));
  return d;
}

// _____________________________________________________________________________
template <typename T>
double innerProd(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
  return innerProd(a.getX(), a.getY(), b.getX(), b.getY(), c.getX(), c.getY());
}

// _____________________________________________________________________________
template <typename T>
double crossProd(const Point<T>& a, const Point<T>& b) {
  return crossProd(a.getX(), a.getY(), b.getX(), b.getY());
}

// _____________________________________________________________________________
template <typename T>
double crossProd(const Point<T>& p, const LineSegment<T>& ls) {
  return crossProd(
      Point<T>(ls.second.getX() - ls.first.getX(),
               ls.second.getY() - ls.first.getY()),
      Point<T>(p.getX() - ls.first.getX(), p.getY() - ls.first.getY()));
}

// _____________________________________________________________________________
template <typename T, typename DF, typename PF>
double withinDist(const Polygon<T>& poly1, const Polygon<T>& poly2,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  if (poly1.size() * poly2.size() > EST_CHECKS_THRESHOLD_XSORTED) {
    return withinDist(XSortedPolygon<T>(poly1), XSortedPolygon<T>(poly2),
                      maxDist, paddingFunc, maxEuclideanDist, distFunc);
  }

  if (intersects(poly1, poly2) || intersects(poly2, poly1)) return 0;

  double d = dist(poly1.getOuter(), poly2.getOuter(), distFunc);

  for (const auto& inner1 : poly1.getInners()) {
    d = std::min(d, dist(poly2.getOuter(), inner1, distFunc));
  }

  for (const auto& inner2 : poly2.getInners()) {
    d = std::min(d, dist(poly1.getOuter(), inner2, distFunc));
  }

  for (const auto& inner1 : poly1.getInners()) {
    for (const auto& inner2 : poly2.getInners()) {
      d = std::min(d, dist(inner1, inner2, distFunc));
    }
  }

  return d;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T, typename DF,
          typename PF>
double withinDist(const Geometry<T>& geom,
                  const std::vector<Geometry<T>>& multi, double maxDist,
                  PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  return withinDist(multi, geom, maxDist, paddingFunc, maxEuclideanDist,
                    distFunc);
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T, typename DF,
          typename PF>
double withinDist(const std::vector<GeometryA<T>>& multi,
                  const GeometryB<T>& geom, double maxDist, PF&& paddingFunc,
                  double maxEuclideanDist, DF&& distFunc) {
  // for larger multigeometries, fall back to pre-sorted implementation to avoid
  // unnecessary checks
  if (multi.size() > EST_MULTI_CHECKS_THRESHOLD_XSORTED)
    return withinDist(XSortedCollection<T>(multi), XSortedCollection<T>(geom),
                      maxDist, paddingFunc, maxEuclideanDist, distFunc);

  double minD = nextafter(maxDist, std::numeric_limits<double>::infinity());
  for (const auto& g1 : multi) {
    minD = std::min(minD, withinDist(g1, geom, minD, paddingFunc,
                                     maxEuclideanDist, distFunc));
  }

  return minD;
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T, typename DF,
          typename PF>
double withinDist(const std::vector<GeometryA<T>>& multi1,
                  const std::vector<GeometryB<T>>& multi2, double maxDist,
                  PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  // for larger multigeometries, fall back to pre-sorted implementation to avoid
  // quadratic pairwise checks
  if (multi1.size() * multi2.size() > EST_MULTI_CHECKS_THRESHOLD_XSORTED)
    return withinDist(XSortedCollection<T>(multi1),
                      XSortedCollection<T>(multi2), maxDist, paddingFunc,
                      maxEuclideanDist, distFunc);

  double minD = nextafter(maxDist, std::numeric_limits<double>::infinity());
  for (const auto& g1 : multi1) {
    for (const auto& g2 : multi2) {
      minD = std::min(minD, withinDist(g1, g2, minD, paddingFunc,
                                       maxEuclideanDist, distFunc));
    }
  }

  return minD;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Polygon<T>& poly1, const Polygon<T>& poly2, DF&& distFunc) {
  return withinDist(poly1, poly2, std::numeric_limits<double>::max(),
                    defaultPaddingFunc<T>(), std::numeric_limits<double>::max(),
                    distFunc);
}

// _____________________________________________________________________________
template <typename T>
double dist(const Polygon<T>& poly1, const Polygon<T>& poly2) {
  return withinDist(poly1, poly2, std::numeric_limits<double>::max());
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Line<T>& l, const Polygon<T>& poly, DF&& distFunc) {
  return withinDist(l, poly, std::numeric_limits<double>::max(),
                    defaultPaddingFunc<T>(), std::numeric_limits<double>::max(),
                    distFunc);
}

// _____________________________________________________________________________
template <typename T>
double dist(const Line<T>& l, const Polygon<T>& poly) {
  if (l.size() * poly.size() > EST_CHECKS_THRESHOLD_XSORTED) {
    return dist(XSortedLine<T>(l), XSortedPolygon<T>(poly));
  }
  return withinDist(l, poly, std::numeric_limits<double>::max());
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Polygon<T>& poly, const Line<T>& l, DF&& distFunc) {
  return dist(l, poly, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const Polygon<T>& poly, const Line<T>& l, double maxDist,
                  PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  if (l.size() * poly.size() > EST_CHECKS_THRESHOLD_XSORTED) {
    return withinDist(XSortedPolygon<T>(poly), XSortedLine<T>(l), maxDist,
                      paddingFunc, maxEuclideanDist, distFunc);
  }

  if (intersects(l, poly)) return 0;
  double d = dist(l, poly.getOuter());

  for (const auto& inner : poly.getInners()) {
    d = std::min(d, dist(l, inner, distFunc));
  }

  return d;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Point<T>& p, const Polygon<T>& poly, DF&& distFunc) {
  return withinDist(p, poly, std::numeric_limits<double>::max(),
                    defaultPaddingFunc<T>(), std::numeric_limits<double>::max(),
                    distFunc);
}

// _____________________________________________________________________________
template <typename T>
double dist(const Point<T>& p, const Polygon<T>& poly) {
  return withinDist(p, poly, std::numeric_limits<double>::max());
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Polygon<T>& poly, const Point<T>& p, DF&& distFunc) {
  return dist(p, poly, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const Point<T>& p, const Polygon<T>& poly, double maxDist,
                  PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  if (poly.size() > EST_CHECKS_THRESHOLD_XSORTED) {
    return withinDist(p, XSortedPolygon<T>(poly), maxDist, paddingFunc,
                      maxEuclideanDist, distFunc);
  }

  if (contains(p, poly)) return 0;
  double d = dist(p, poly.getOuter(), distFunc);

  for (const auto& inner : poly.getInners()) {
    d = std::min(d, dist(p, inner, distFunc));
  }

  return d;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double dist(const Point<T>& p1, const Point<T>& p2, DF&& distFunc) {
  return distFunc(p1, p2, std::numeric_limits<T>::max());
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const Point<T>& p1, const Point<T>& p2, double, PF&&, double,
                  DF&& distFunc) {
  return distFunc(p1, p2, std::numeric_limits<T>::max());
}

// _____________________________________________________________________________
template <typename T>
size_t numPoints(const Point<T>&) {
  return 1;
}

// _____________________________________________________________________________
template <typename T>
size_t numPoints(const Line<T>& l) {
  return l.size();
}

// _____________________________________________________________________________
template <typename T>
size_t numPoints(const Polygon<T>& p) {
  size_t ret = p.getOuter().size();
  for (const auto& i : p.getInners()) ret += i.size();
  return ret;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
size_t numPoints(const std::vector<Geometry<T>>& pol) {
  size_t ret = 0;
  for (const auto& g : pol) ret += numPoints(g);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
size_t numPoints(const Collection<T>& collection) {
  size_t ret = 0;
  for (const auto& g : collection) {
    if (g.getType() == 0) ret += numPoints(g.getPoint());
    if (g.getType() == 1) ret += numPoints(g.getLine());
    if (g.getType() == 2) ret += numPoints(g.getPolygon());
    if (g.getType() == 3) ret += numPoints(g.getMultiLine());
    if (g.getType() == 4) ret += numPoints(g.getMultiPolygon());
    if (g.getType() == 5) ret += numPoints(g.getCollection());
  }
  return ret;
}

// _____________________________________________________________________________
template <typename T>
bool empty(const Point<T>&) {
  return false;
}

// _____________________________________________________________________________
template <typename T>
bool empty(const Polygon<T>& g) {
  return g.getOuter().empty() && g.getInners().empty();
}

// _____________________________________________________________________________
template <typename T>
bool empty(const Line<T>& g) {
  return g.empty();
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
size_t empty(const std::vector<Geometry<T>>& pol) {
  for (const auto& g : pol) {
    if (!empty(g)) return false;
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
bool empty(const Collection<T>& g) {
  return g.empty();
}

// _____________________________________________________________________________
template <typename T, typename F>
Line<T> lineFromWKTProj(const char* c, const char** endr, F projFunc) {
  Line<T> line;

  c = strchr(c, '(');
  if (!c) {
    if (endr) (*endr) = 0;
    return line;  // parse error
  }
  c++;

  auto end = strchr(c, ')');
  if (endr) (*endr) = end;
  if (!end) {
    if (endr) (*endr) = 0;
    return line;  // parse error
  }

  line.reserve((end - c) / 20);

  while (true) {
    while (*c && *c != ')' &&
           (*c == ' ' || *c == '\n' || *c == '\t' || *c == '\r'))
      c++;

    double x = util::atof(c, 10);

    const char* next = strchr(c, ' ');

    if (!next || next >= end) return {};  // parse error

    while (*next && *next != ')' &&
           (*next == ' ' || *next == '\n' || *next == '\t' || *next == '\r'))
      next++;

    double y = util::atof(next, 10);

    line.push_back(projFunc(util::geo::DPoint(x, y)));

    auto n = strchr(next, ',');
    if (!n || n > end) break;
    c = n + 1;
  }

  return line;
}

// _____________________________________________________________________________
template <typename T>
Line<T> lineFromWKT(const char* c, const char** endr) {
  return lineFromWKTProj<T>(c, endr, [](const Point<double>& p) {
    return Point<T>{static_cast<T>(p.getX()), static_cast<T>(p.getY())};
  });
}

// _____________________________________________________________________________
template <typename T>
MultiLine<T> multiLineFromWKT(const char* c, const char** endr) {
  return multiLineFromWKTProj<T>(c, endr, [](const Point<double>& p) {
    return Point<T>{static_cast<T>(p.getX()), static_cast<T>(p.getY())};
  });
}

// _____________________________________________________________________________
template <typename T, typename F>
MultiPoint<T> multiPointFromWKTProj(const char* c, const char** endr,
                                    F projFunc) {
  // try MULTIPOINT((1 1), (2 2)) syntax
  const auto& mline = multiLineFromWKTProj<T, F>(c, endr, projFunc);

  MultiPoint<T> ret;
  for (const auto& l : mline) {
    if (l.size() == 1) ret.push_back(l[0]);
  }

  if (ret.size()) return ret;

  // try MULTIPOINT(1 1, 2 2) syntax
  const auto& line = lineFromWKTProj<T, F>(c, endr, projFunc);
  if (line.size() > 0) return MultiPoint<T>(std::move(line));

  return ret;
}

// _____________________________________________________________________________
template <typename T>
MultiPoint<T> multiPointFromWKT(const char* c, const char** endr) {
  return multiPointFromWKTProj<T>(c, endr, [](const Point<double>& p) {
    return Point<T>{static_cast<T>(p.getX()), static_cast<T>(p.getY())};
  });
}

// _____________________________________________________________________________
template <typename T>
MultiPoint<T> multiPointFromWKT(const std::string& wkt) {
  return multiPointFromWKT<T>(wkt.c_str(), 0);
}

// _____________________________________________________________________________
template <typename T, typename F>
MultiPoint<T> multiPointFromWKTProj(const std::string& wkt, F&& projFunc) {
  return multiPointFromWKTProj<T>(wkt.c_str(), 0, projFunc);
}

// _____________________________________________________________________________
template <typename T, typename F>
Point<T> pointFromWKTProj(const char* c, const char** endr, F projFunc) {
  c = strchr(c, '(');
  if (!c) {
    if (endr) (*endr) = 0;
    return {0, 0};
  }

  c += 1;
  while (*c && *c != ')' &&
         (*c == ' ' || *c == '\n' || *c == '\t' || *c == '\r'))
    c++;

  double x = util::atof(c, 10);
  const char* next = strchr(c, ' ');
  if (!next) return {0, 0};  // TODO!
  while (*next && *next != ')' &&
         (*next == ' ' || *next == '\n' || *next == '\t' || *next == '\r'))
    next++;
  double y = util::atof(next, 10);

  if (endr) (*endr) = strchr(next, ')');

  return projFunc(util::geo::DPoint(x, y));
}

// _____________________________________________________________________________
template <typename T>
Point<T> pointFromWKT(const char* c, const char** endr) {
  return pointFromWKTProj<T>(c, endr, [](const Point<double>& p) {
    return Point<T>{static_cast<T>(p.getX()), static_cast<T>(p.getY())};
  });
}

// _____________________________________________________________________________
template <typename T>
Point<T> pointFromWKT(std::string wkt) {
  return pointFromWKT<T>(wkt.c_str(), 0);
}

// _____________________________________________________________________________
template <typename T, typename F>
Point<T> pointFromWKTProj(std::string wkt, F projFunc) {
  return pointFromWKTProj<T>(wkt.c_str(), 0, projFunc);
}

// _____________________________________________________________________________
template <typename T, typename F>
Polygon<T> polygonFromWKTProj(const char* c, const char** endr, F projFunc) {
  c = strchr(c, '(');
  if (!c) {
    if (endr) (*endr) = 0;
    return {};  // parse error
  }
  c += 1;

  size_t i = 0;
  Polygon<T> poly;
  while ((c = strchr(c, '('))) {
    const char* end = 0;
    const auto& line = lineFromWKTProj<T, F>(c, &end, projFunc);

    if (!end) {
      if (endr) (*endr) = 0;
      return {};  // parse error
    }

    c = end;

    if (i == 0)
      poly.getOuter() = line;
    else
      poly.getInners().push_back(std::move(line));
    i++;

    auto q = strchr(c + 1, ')');
    auto cc = strchr(c + 1, '(');

    if ((!cc && q) || (q && cc && q < cc)) {
      // polygon closes at q
      if (endr) (*endr) = q;
      poly.fix();
      return poly;
    }

    if (cc) c = cc;
  }

  poly.fix();
  return poly;
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> polygonFromWKT(const char* c, const char** endr) {
  return polygonFromWKTProj<T>(c, endr, [](const Point<double>& p) {
    return Point<T>{static_cast<T>(p.getX()), static_cast<T>(p.getY())};
  });
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> polygonFromWKT(std::string wkt) {
  return polygonFromWKT<T>(wkt.c_str(), 0);
}

// _____________________________________________________________________________
template <typename T, typename F>
Polygon<T> polygonFromWKTProj(std::string wkt, F projFunc) {
  return polygonFromWKTProj<T>(wkt.c_str(), 0, projFunc);
}

// _____________________________________________________________________________
template <typename T, typename F>
MultiLine<T> multiLineFromWKTProj(const char* c, const char** endr,
                                  F projFunc) {
  c = strchr(c, '(');
  if (!c) {
    if (endr) (*endr) = 0;
    return {};  // parse error
  }
  c += 1;

  MultiLine<T> ml;
  while ((c = strchr(c, '('))) {
    const char* end = 0;
    const auto& line = lineFromWKTProj<T, F>(c, &end, projFunc);
    if (!end) break;
    if (line.size() != 0) ml.push_back(std::move(line));

    auto nextComma = strchr(end + 1, ',');
    auto nextCloseBracket = strchr(end + 1, ')');

    if (!nextComma ||
        (nextComma && nextCloseBracket && nextComma > nextCloseBracket)) {
      if (endr) (*endr) = nextCloseBracket;
      return ml;
    }

    c = nextComma;
  }

  return ml;
}

// _____________________________________________________________________________
template <typename T, typename F>
MultiPolygon<T> multiPolygonFromWKTProj(const char* c, const char** endr,
                                        F projFunc) {
  c = strchr(c, '(');
  if (!c) {
    if (endr) (*endr) = 0;
    return {};  // parse error
  }
  c += 1;

  MultiPolygon<T> mp;
  do {
    c = strchr(c, '(');
    if (!c) break;
    const char* end = 0;
    const auto& poly = polygonFromWKTProj<T, F>(c, &end, projFunc);

    if (!end) break;

    if (poly.getOuter().size() > 1) mp.push_back(std::move(poly));

    auto nextComma = strchr(end + 1, ',');
    auto nextCloseBracket = strchr(end + 1, ')');

    if (!nextComma ||
        (nextComma && nextCloseBracket && nextComma > nextCloseBracket)) {
      if (endr) (*endr) = nextCloseBracket;
      return mp;
    }

    c = nextComma;
  } while (c);

  return mp;
}

// _____________________________________________________________________________
template <typename T>
MultiPolygon<T> multiPolygonFromWKT(const char* c, const char** endr) {
  return multiPolygonFromWKTProj<T>(c, endr, [](const Point<double>& p) {
    return Point<T>{static_cast<T>(p.getX()), static_cast<T>(p.getY())};
  });
}

// _____________________________________________________________________________
template <typename T, typename F>
Collection<T> collectionFromWKTProj(const char* c, const char** endr,
                                    F&& projFunc) {
  Collection<T> col;

  c = strchr(c, '(');
  if (!c) {
    if (endr) (*endr) = 0;
    return col;
  }
  do {
    c++;
    while ((*c == ' ' || *c == '\n' || *c == '\t' || *c == '\r'))
      c++;  // skip possible whitespace

    auto wktType = getWKTType(c, &c);

    if (wktType == NONE) {
      if (endr) (*endr) = 0;
      return {};
    }

    if (wktType == POINT) {
      const char* end = 0;
      const auto& point = pointFromWKTProj<T, F>(c, &end, projFunc);

      if (!end) {
        if (endr) (*endr) = 0;
        return {};
      }

      col.push_back(point);
      c = const_cast<char*>(strchr(end, ','));
    } else if (wktType == POLYGON) {
      const char* end = 0;
      const auto& poly = polygonFromWKTProj<T, F>(c, &end, projFunc);

      if (!end) {
        if (endr) (*endr) = 0;
        return {};
      }
      if (poly.getOuter().size() > 1) col.push_back(poly);
      c = const_cast<char*>(strchr(end, ','));
    } else if (wktType == LINESTRING) {
      const char* end = 0;
      const auto& line = lineFromWKTProj<T, F>(c, &end, projFunc);

      if (!end) {
        if (endr) (*endr) = 0;
        return {};
      }
      if (line.size() > 1) col.push_back(line);
      c = const_cast<char*>(strchr(end, ','));
    } else if (wktType == MULTIPOINT) {
      const char* end = 0;
      const auto& line = lineFromWKTProj<T, F>(c, &end, projFunc);

      if (!end) {
        if (endr) (*endr) = 0;
        return {};
      }
      if (line.size() > 0) col.push_back(MultiPoint<T>(std::move(line)));
      c = const_cast<char*>(strchr(end, ','));
    } else if (wktType == MULTIPOLYGON) {
      const char* end = 0;
      const auto& mp = multiPolygonFromWKTProj<T, F>(c, &end, projFunc);

      if (!end) {
        if (endr) (*endr) = 0;
        return {};
      }
      if (mp.size()) col.push_back(mp);
      c = const_cast<char*>(strchr(end, ','));
    } else if (wktType == MULTILINESTRING) {
      const char* end = 0;
      const auto& ml = multiLineFromWKTProj<T, F>(c, &end, projFunc);

      if (!end) {
        if (endr) (*endr) = 0;
        return {};
      }
      if (ml.size()) col.push_back(ml);
      c = const_cast<char*>(strchr(end, ','));
    }
  } while (c && *c);

  if (endr) (*endr) = strchr(c, ')');

  return col;
}

// _____________________________________________________________________________
template <typename T>
Collection<T> collectionFromWKT(const char* c, const char** endr) {
  return collectionFromWKTProj<T>(c, endr, [](const Point<double>& p) {
    return Point<T>{static_cast<T>(p.getX()), static_cast<T>(p.getY())};
  });
}

// _____________________________________________________________________________
template <typename T>
Line<T> lineFromWKT(const std::string& wkt) {
  return lineFromWKT<T>(wkt.c_str(), 0);
}

// _____________________________________________________________________________
template <typename T, typename F>
Line<T> lineFromWKTProj(const std::string& wkt, F&& projFunc) {
  return lineFromWKTProj<T>(wkt.c_str(), 0, projFunc);
}

// _____________________________________________________________________________
template <typename T>
MultiLine<T> multiLineFromWKT(const std::string& wkt) {
  return multiLineFromWKT<T>(wkt.c_str(), 0);
}

// _____________________________________________________________________________
template <typename T, typename F>
MultiLine<T> multiLineFromWKTProj(const std::string& wkt, F&& projFunc) {
  return multiLineFromWKTProj<T>(wkt.c_str(), 0, projFunc);
}

// _____________________________________________________________________________
template <typename T>
MultiPolygon<T> multiPolygonFromWKT(const std::string& wkt) {
  return multiPolygonFromWKT<T>(wkt.c_str(), 0);
}

// _____________________________________________________________________________
template <typename T, typename F>
MultiPolygon<T> multiPolygonFromWKTProj(const std::string& wkt, F&& projFunc) {
  return multiPolygonFromWKTProj<T>(wkt.c_str(), 0, projFunc);
}

// _____________________________________________________________________________
template <typename T>
Collection<T> collectionFromWKT(const std::string& wkt) {
  return collectionFromWKT<T>(wkt.c_str(), 0);
}

// _____________________________________________________________________________
template <typename T, typename F>
Collection<T> collectionFromWKTProj(const std::string& wkt, F&& projFunc) {
  return collectionFromWKTProj<T>(wkt.c_str(), 0, projFunc);
}

// _____________________________________________________________________________
template <typename T>
double len(const Point<T>&) {
  return 0;
}

// _____________________________________________________________________________
template <typename T>
double len(const Line<T>& g) {
  double ret = 0;
  for (size_t i = 1; i < g.size(); i++) ret += dist(g[i - 1], g[i]);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
double len(const Polygon<T>& g) {
  double ret = 0;
  for (size_t i = 1; i < g.getOuter().size(); i++)
    ret += dist(g.getOuter()[i - 1], g.getOuter()[i]);
  for (const auto& inner : g.getInners()) {
    for (size_t i = 1; i < inner.size(); i++)
      ret += dist(inner[i - 1], inner[i]);
  }
  return ret;
}

// _____________________________________________________________________________
template <typename T>
double len(const std::vector<XSortedTuple<T>>& g) {
  double ret = 0;
  for (size_t i = 0; i < g.size(); i++) ret += len(g[i].seg());
  return ret / 2;
}

// _____________________________________________________________________________
template <typename T>
double len(const XSortedLine<T>& g) {
  return g.length();
}

// _____________________________________________________________________________
template <typename T>
double len(const LineSegment<T>& g) {
  return dist(g.first, g.second);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
size_t len(const std::vector<Geometry<T>>& pol) {
  double ret = 0;
  for (const auto& g : pol) ret += len(g);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
double len(const Collection<T>& collection) {
  double ret = 0;
  for (const auto& g : collection) {
    if (g.getType() == 0) ret += len(g.getPoint());
    if (g.getType() == 1) ret += len(g.getLine());
    if (g.getType() == 2) ret += len(g.getPolygon());
    if (g.getType() == 3) ret += len(g.getMultiLine());
    if (g.getType() == 4) ret += len(g.getMultiPolygon());
    if (g.getType() == 5) ret += len(g.getCollection());
  }
  return ret;
}

// _____________________________________________________________________________
template <typename T>
bool shorterThan(const Line<T>& g, double d) {
  double ret = 0;
  for (size_t i = 1; i < g.size(); i++) {
    ret += dist(g[i - 1], g[i]);
    if (ret >= d) return false;
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
bool longerThan(const Line<T>& g, double d) {
  double ret = 0;
  for (size_t i = 1; i < g.size(); i++) {
    ret += dist(g[i - 1], g[i]);
    if (ret > d) return true;
  }
  return false;
}

// _____________________________________________________________________________
template <typename T>
bool longerThan(const Line<T>& a, const Line<T>& b, double d) {
  double ret = 0;
  for (size_t i = 1; i < a.size(); i++) {
    ret += dist(a[i - 1], a[i]);
    if (ret > d) return true;
  }
  for (size_t i = 1; i < b.size(); i++) {
    ret += dist(b[i - 1], b[i]);
    if (ret > d) return true;
  }
  return false;
}

// _____________________________________________________________________________
template <typename T>
Point<T> simplify(const Point<T>& g, double) {
  return g;
}

// _____________________________________________________________________________
template <typename T>
LineSegment<T> simplify(const LineSegment<T>& g, double) {
  return g;
}

// _____________________________________________________________________________
template <typename T>
Box<T> simplify(const Box<T>& g, double) {
  return g;
}

// _____________________________________________________________________________
template <typename T>
RotatedBox<T> simplify(const RotatedBox<T>& g, double) {
  return g;
}

// _____________________________________________________________________________
template <typename T>
void simplifyInPlaceHelper(Line<T>* g, size_t from, size_t to, double d) {
  // douglas peucker
  double maxd = 0;
  size_t maxi = 0;
  for (size_t i = from + 1; i < to - 1; i++) {
    double dt = distToSegmentSquared((*g)[from], (*g)[to - 1], (*g)[i]);
    if (dt > maxd) {
      maxi = i;
      maxd = dt;
    }
  }

  if (maxd > d * d && maxi != 0) {
    // maxi is kept, simplify lower and upper half
    simplifyInPlaceHelper(g, from, maxi + 1, d);
    simplifyInPlaceHelper(g, maxi, to, d);
  } else {
    // no distance exceeds d, delete all (by setting it to the first,
    // will be deleted later on)
    for (size_t i = from + 1; i + 1 < to; i++) (*g)[i] = (*g)[from];
  }
}

// _____________________________________________________________________________
template <typename T>
Line<T> simplify(Line<T> g, double d) {
  if (g.size() < 3) return g;

  simplifyInPlaceHelper(&g, 0, g.size(), d);

  // deleted points have been stored in-place by setting them to their
  // previous point, remove them in-place
  size_t w = 1;
  for (size_t i = 1; i < g.size(); i++) {
    if (g[i] == g[w - 1]) continue;
    if (w != i) g[w] = g[i];
    w++;
  }

  g.resize(w);

  return g;
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> simplify(const Polygon<T>& g, double d) {
  auto simple = simplify(g.getOuter(), d);
  std::rotate(simple.begin(), simple.begin() + simple.size() / 2, simple.end());
  simple = simplify(std::move(simple), d);
  Polygon<T> ret(simple);

  for (const auto& inner : g.getInners()) {
    auto simple = simplify(inner, d);
    std::rotate(simple.begin(), simple.begin() + simple.size() / 2,
                simple.end());
    simple = simplify(std::move(simple), d);
    ret.getInners().push_back(simple);
  }

  return ret;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
std::vector<Geometry<T>> simplify(const std::vector<Geometry<T>>& pol,
                                  double d) {
  std::vector<Geometry<T>> ret;
  ret.reserve(pol.size());
  for (const auto& g : pol) ret.push_back(std::move(simplify(g, d)));
  return ret;
}

// _____________________________________________________________________________
template <typename T>
Collection<T> simplify(const Collection<T>& collection, double d) {
  Collection<T> ret;
  for (const auto& g : collection) {
    if (g.getType() == 0) ret.push_back(std::move(simplify(g.getPoint(), d)));
    if (g.getType() == 1) ret.push_back(std::move(simplify(g.getLine(), d)));
    if (g.getType() == 2) ret.push_back(std::move(simplify(g.getPolygon(), d)));
    if (g.getType() == 3)
      ret.push_back(std::move(simplify(g.getMultiLine(), d)));
    if (g.getType() == 4)
      ret.push_back(std::move(simplify(g.getMultiPolygon(), d)));
    if (g.getType() == 5)
      ret.push_back(std::move(simplify(g.getCollection(), d)));
  }

  return ret;
}

// _____________________________________________________________________________
template <typename T, typename DF>
double distToSegment(T lax, T lay, T lbx, T lby, T px, T py, DF&& distFunc) {
  double d = distFunc(Point<T>{lax, lay}, Point<T>{lbx, lby},
                      std::numeric_limits<double>::max()) *
             distFunc(Point<T>{lax, lay}, Point<T>{lbx, lby},
                      std::numeric_limits<double>::max());
  if (d == 0)
    return distFunc(Point<T>{px, py}, Point<T>{lax, lay},
                    std::numeric_limits<double>::max());

  double t = ((px - lax) * (lbx - lax) + (py - lay) * (lby - lay)) / d;

  if (t < 0) {
    return distFunc(Point<T>{px, py}, Point<T>{lax, lay},
                    std::numeric_limits<double>::max());
  } else if (t > 1) {
    return distFunc(Point<T>{px, py}, Point<T>{lbx, lby},
                    std::numeric_limits<double>::max());
  }

  return distFunc(Point<T>{px, py},
                  Point<T>{static_cast<T>(lax + t * (lbx - lax)),
                           static_cast<T>(lay + t * (lby - lay))},
                  std::numeric_limits<double>::max());
}

// _____________________________________________________________________________
template <typename T, typename DF>
double distToSegment(const Point<T>& la, const Point<T>& lb, const Point<T>& p,
                     DF&& distFunc) {
  return distToSegment(la.getX(), la.getY(), lb.getX(), lb.getY(), p.getX(),
                       p.getY(), distFunc);
}

// _____________________________________________________________________________
template <typename T, typename DF>
double distToSegment(const LineSegment<T>& ls, const Point<T>& p,
                     DF&& distFunc) {
  return distToSegment(ls.first.getX(), ls.first.getY(), ls.second.getX(),
                       ls.second.getY(), p.getX(), p.getY(), distFunc);
}

// _____________________________________________________________________________
template <typename T>
double distToSegment(const Point<T>& la, const Point<T>& lb,
                     const Point<T>& p) {
  return distToSegment(la.getX(), la.getY(), lb.getX(), lb.getY(), p.getX(),
                       p.getY());
}

// _____________________________________________________________________________
template <typename T>
double distToSegment(const LineSegment<T>& ls, const Point<T>& p) {
  return distToSegment(ls.first.getX(), ls.first.getY(), ls.second.getX(),
                       ls.second.getY(), p.getX(), p.getY());
}

// _____________________________________________________________________________
template <typename T, typename DF>
double distToSegmentSquared(const LineSegment<T>& ls, const Point<T>& p,
                            DF&& distFunc) {
  return distToSegmentSquared(ls.first.getX(), ls.first.getY(),
                              ls.second.getX(), ls.second.getY(), p.getX(),
                              p.getY(), distFunc);
}

// _____________________________________________________________________________
template <typename T>
double distToSegmentSquared(const Point<T>& la, const Point<T>& lb,
                            const Point<T>& p) {
  return distToSegmentSquared(la.getX(), la.getY(), lb.getX(), lb.getY(),
                              p.getX(), p.getY());
}

// _____________________________________________________________________________
template <typename T>
double distToSegmentSquared(const LineSegment<T>& ls, const Point<T>& p) {
  return distToSegmentSquared(ls.first.getX(), ls.first.getY(),
                              ls.second.getX(), ls.second.getY(), p.getX(),
                              p.getY());
}

// _____________________________________________________________________________
template <typename T>
Point<T> projectOn(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
  if (doubleEq(a.getX(), b.getX()) && doubleEq(a.getY(), b.getY())) return a;
  if (doubleEq(a.getX(), c.getX()) && doubleEq(a.getY(), c.getY())) return a;
  if (doubleEq(b.getX(), c.getX()) && doubleEq(b.getY(), c.getY())) return b;

  double x, y;

  if (c.getX() == a.getX()) {
    // infinite slope
    x = a.getX();
    y = b.getY();
  } else {
    double m = (double)(1.0 * c.getY() - 1.0 * a.getY()) /
               (1.0 * c.getX() - 1.0 * a.getX());
    double bb = (double)a.getY() - (m * a.getX());

    x = (m * b.getY() + 1.0 * b.getX() - m * bb) / (m * m + 1.0);
    y = (m * m * b.getY() + m * b.getX() + bb) / (m * m + 1.0);
  }

  Point<T> ret = Point<T>(x, y);

  bool isBetween = dist(a, c) > dist(a, ret) && dist(a, c) > dist(c, ret);
  bool nearer = dist(a, ret) < dist(c, ret);

  if (!isBetween) return nearer ? a : c;

  return ret;
}

// _____________________________________________________________________________
template <typename T>
double parallelity(const Box<T>& box, const Line<T>& line) {
  double ret = M_PI;

  double a = angBetween(
      box.getLowerLeft(),
      Point<T>(box.getLowerLeft().getX(), box.getUpperRight().getY()));
  double b = angBetween(
      box.getLowerLeft(),
      Point<T>(box.getUpperRight().getX(), box.getLowerLeft().getY()));
  double c = angBetween(
      box.getUpperRight(),
      Point<T>(box.getLowerLeft().getX(), box.getUpperRight().getY()));
  double d = angBetween(
      box.getUpperRight(),
      Point<T>(box.getUpperRight().getX(), box.getLowerLeft().getY()));

  double e = angBetween(line.front(), line.back());

  double vals[] = {a, b, c, d};

  for (double ang : vals) {
    double v = fabs(ang - e);
    if (v > M_PI) v = 2 * M_PI - v;
    if (v > M_PI_2) v = M_PI - v;
    if (v < ret) ret = v;
  }

  return 1 - (ret / (M_PI_4));
}

// _____________________________________________________________________________
template <typename T>
double parallelity(const Box<T>& box, const MultiLine<T>& multiline) {
  double ret = 0;
  for (const Line<T>& l : multiline) {
    ret += parallelity(box, l);
  }

  return ret / static_cast<float>(multiline.size());
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
RotatedBox<T> getOrientedEnvelope(const std::vector<Geometry<T>>& pol,
                                  const std::vector<double>& angles) {
  const Point<T> center = centroid(pol);
  Box<T> tmpBox = getBoundingBox(pol);
  double tmpArea = area(tmpBox);
  double rotateAngle = 0;

  // check each segment
  for (size_t i = 0; i < angles.size(); i++) {
    // rotate segment such that it is parallel to the x axis
    const auto p = rotateRAD(pol, angles[i], center);
    const Box<T>& e = getBoundingBox(p);
    const double newArea = area(e);
    if (tmpArea > newArea) {
      tmpBox = e;
      tmpArea = newArea;
      rotateAngle = angles[i];
    }
  }

  return RotatedBox<T>(tmpBox, -rotateAngle * -IRAD, center);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
RotatedBox<T> getOrientedEnvelope(const std::vector<Geometry<T>>& pol,
                                  double step) {
  std::vector<double> angles;
  double i = 0;
  while (i < 360) {
    angles.push_back(i * RAD);
    i += step;
  }

  return getOrientedEnvelope(pol, angles);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
RotatedBox<T> getOrientedEnvelope(const std::vector<Geometry<T>>& pol) {
  Line<T> hull = convexHull(pol).getOuter();
  std::vector<double> angles;

  angles.reserve(hull.size());

  for (size_t i = 1; i < hull.size(); i++) {
    const auto s = hull[i] - hull[i - 1];
    double ang = -std::atan2(s.getY(), s.getX());
    if (ang < 0) ang = M_PI + ang;
    if (fabs(ang) > EPSILON && fabs(ang - M_PI_2) > EPSILON)
      angles.push_back(ang);
  }

  // Check segment between the ends of the hull line
  const auto s = hull[0] - hull[hull.size() - 1];
  double ang = -std::atan2(s.getY(), s.getX());
  if (ang < 0) ang = M_PI + ang;
  if (fabs(ang) > EPSILON && fabs(ang - M_PI_2) > EPSILON)
    angles.push_back(ang);

  auto end = angles.end();

  // only interested in unique angles!
  std::sort(angles.begin(), angles.end());
  end = std::unique(angles.begin(), angles.end());

  return getOrientedEnvelope(pol, angles);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
RotatedBox<T> getOrientedEnvelope(const Geometry<T>& pol) {
  return getOrientedEnvelope(std::vector<Geometry<T>>{pol});
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
RotatedBox<T> getOrientedEnvelope(const Geometry<T>& pol,
                                  const std::vector<double>& angles) {
  return getOrientedEnvelope(std::vector<Geometry<T>>{pol}, angles);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
RotatedBox<T> getOrientedEnvelope(const Geometry<T>& pol, double step) {
  return getOrientedEnvelope(std::vector<Geometry<T>>{pol}, step);
}

// _____________________________________________________________________________
template <typename T>
RotatedBox<T> getOrientedEnvelope(const Collection<T>& collection) {
  MultiPolygon<T> p;
  for (const auto& g : collection) {
    if (g.getType() == 0) p.push_back(convexHull(g.getPoint()));
    if (g.getType() == 1) p.push_back(convexHull(g.getLine()));
    if (g.getType() == 2) p.push_back(convexHull(g.getPolygon()));
    if (g.getType() == 3) p.push_back(convexHull(g.getMultiLine()));
    if (g.getType() == 4) p.push_back(convexHull(g.getMultiPolygon()));
    if (g.getType() == 5) p.push_back(convexHull(g.getCollection()));
  }
  return getOrientedEnvelope(p);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> buffer(const Line<T>& line, double d, size_t points) {
  // so far only works correctly if the input polygon is convex
  MultiPoint<T> pSet;
  for (const auto& p : line) {
    Point<T> anchor{p.getX() + d, p.getY()};
    double deg = 0;
    pSet.push_back(p);
    for (size_t i = 0; i < points; i++) {
      pSet.push_back(rotate(anchor, deg, p));
      deg += 360 / (1.0 * points);
    }
  }

  return convexHull(pSet);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> buffer(const Polygon<T>& pol, double d, size_t points) {
  return buffer(pol.getOuter(), d, points);
}

// _____________________________________________________________________________
template <typename T>
Box<T> getBoundingBox(const Point<T>& p) {
  return Box<T>(p, p);
}

// _____________________________________________________________________________
template <typename T>
Box<T> getBoundingBox(const Line<T>& l) {
  Box<T> ret;
  for (const auto& p : l) ret = extendBox(p, ret);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
Box<T> getBoundingBox(const Polygon<T>& pol) {
  Box<T> ret;
  for (const auto& p : pol.getOuter()) ret = extendBox(p, ret);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
Box<T> getBoundingBox(const LineSegment<T>& ls) {
  return Box<T>({std::min(ls.first.getX(), ls.second.getX()),
                 std::min(ls.first.getY(), ls.second.getY())},
                {std::max(ls.first.getX(), ls.second.getX()),
                 std::max(ls.first.getY(), ls.second.getY())});
}

// _____________________________________________________________________________
template <typename T>
Box<T> getBoundingBox(const Box<T>& b) {
  return b;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
Box<T> getBoundingBox(const std::vector<Geometry<T>>& multigeo) {
  Box<T> b;
  b = extendBox(multigeo, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
Box<T> getBoundingBox(const Collection<T>& coll) {
  Box<T> b;
  b = extendBox(coll, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
Box<T> getBoundingRect(const Box<T>& b) {
  auto box = Box<T>();
  auto centroid = util::geo::centroid(b);
  box = extendBox(b, box);
  box = extendBox(rotate(convexHull(b), 180, centroid), box);
  return box;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
Box<T> getBoundingRect(const Geometry<T> geom) {
  return getBoundingRect<T>(getBoundingBox<T>(geom));
}

// _____________________________________________________________________________
template <typename T>
double getEnclosingRadius(const Point<T>& p, const Point<T>& pp) {
  return dist(p, pp);
}

// _____________________________________________________________________________
template <typename T>
double getEnclosingRadius(const Point<T>& p, const Line<T>& l) {
  double ret = 0;
  for (const auto& pp : l)
    if (getEnclosingRadius(p, pp) > ret) ret = getEnclosingRadius(p, pp);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
double getEnclosingRadius(const Point<T>& p, const Polygon<T>& pg) {
  double ret = 0;
  for (const auto& pp : pg.getOuter())
    if (getEnclosingRadius(p, pp) > ret) ret = getEnclosingRadius(p, pp);
  return ret;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
double getEnclosingRadius(const Point<T>& p,
                          const std::vector<Geometry<T>>& multigeom) {
  double ret = 0;
  for (const auto& pp : multigeom)
    if (getEnclosingRadius(p, pp) > ret) ret = getEnclosingRadius(p, pp);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const Point<T>& p) {
  return Polygon<T>({p});
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const Box<T>& b) {
  return Polygon<T>(b);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const LineSegment<T>& b) {
  return Polygon<T>(Line<T>{b.first, b.second});
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const RotatedBox<T>& b) {
  auto p = convexHull(b.getBox());
  p = rotate(p, b.getDegree(), b.getCenter());
  return p;
}

// _____________________________________________________________________________
template <typename T>
size_t convexHullImpl(const MultiPoint<T>& a, size_t p1, size_t p2,
                      Line<T>* h) {
  // emergency stop
  if (h->size() >= a.size() + 1) return 0;

  // quickhull by Barber, Dobkin & Huhdanpaa
  Point<T> pa;
  bool found = false;
  double maxDist = 0;
  for (const auto& p : a) {
    double tmpDist = distToSegment((*h)[p1], (*h)[p2], p);
    double cp = crossProd(p, LineSegment<T>((*h)[p1], (*h)[p2]));
    if ((cp > 0 + EPSILON) && tmpDist > maxDist) {
      pa = p;
      found = true;
      maxDist = tmpDist;
    }
  }

  if (!found) return 0;

  h->insert(h->begin() + p2, pa);
  size_t in = 1 + convexHullImpl(a, p1, p2, h);
  return in + convexHullImpl(a, p2 + in - 1, p2 + in, h);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const MultiPoint<T>& l) {
  if (l.size() == 2) return convexHull(LineSegment<T>(l[0], l[1]));
  if (l.size() == 1) return convexHull(l[0]);

  Point<T> left(std::numeric_limits<T>::max(), std::numeric_limits<T>::max());
  Point<T> right(std::numeric_limits<T>::lowest(),
                 std::numeric_limits<T>::lowest());
  for (const auto& p : l) {
    if (p.getX() < left.getX() ||
        (p.getX() == left.getX() && p.getY() < left.getY()))
      left = p;
    if (p.getX() > right.getX() ||
        (p.getX() == right.getX() && p.getY() > right.getY()))
      right = p;
  }

  Line<T> hull{left, right};
  convexHullImpl(l, 0, 1, &hull);
  hull.push_back(hull.front());
  convexHullImpl(l, hull.size() - 2, hull.size() - 1, &hull);
  hull.pop_back();

  return Polygon<T>(hull);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const Polygon<T>& p) {
  return convexHull(p.getOuter());
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const MultiPolygon<T>& ps) {
  MultiPoint<T> mp;
  for (const auto& p : ps)
    mp.insert(mp.end(), p.getOuter().begin(), p.getOuter().end());
  return convexHull(mp);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const MultiLine<T>& ls) {
  MultiPoint<T> mp;
  for (const auto& l : ls) mp.insert(mp.end(), l.begin(), l.end());
  return convexHull(mp);
}

// _____________________________________________________________________________
template <typename T>
Polygon<T> convexHull(const Collection<T>& collection) {
  MultiPolygon<T> p;
  for (const auto& g : collection) {
    if (g.getType() == 0) p.push_back(convexHull(g.getPoint()));
    if (g.getType() == 1) p.push_back(convexHull(g.getLine()));
    if (g.getType() == 2) p.push_back(convexHull(g.getPolygon()));
    if (g.getType() == 3) p.push_back(convexHull(g.getMultiLine()));
    if (g.getType() == 4) p.push_back(convexHull(g.getMultiPolygon()));
    if (g.getType() == 5) p.push_back(convexHull(g.getCollection()));
  }

  return convexHull(p);
}

// _____________________________________________________________________________
template <typename T>
Box<T> extendBox(const Line<T>& l, Box<T> b) {
  for (const auto& p : l) b = extendBox(p, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
Box<T> extendBox(const LineSegment<T>& ls, Box<T> b) {
  b = extendBox(ls.first, b);
  b = extendBox(ls.second, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
Box<T> extendBox(const Polygon<T>& ls, Box<T> b) {
  return extendBox(ls.getOuter(), b);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
Box<T> extendBox(const std::vector<Geometry<T>>& multigeom, Box<T> b) {
  for (const auto& g : multigeom) b = extendBox(g, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
Box<T> extendBox(const Collection<T>& collection, Box<T> b) {
  for (const auto& g : collection) {
    if (g.getType() == 0) b = extendBox(g.getPoint(), b);
    if (g.getType() == 1) b = extendBox(g.getLine(), b);
    if (g.getType() == 2) b = extendBox(g.getPolygon(), b);
    if (g.getType() == 3) b = extendBox(g.getMultiLine(), b);
    if (g.getType() == 4) b = extendBox(g.getMultiPolygon(), b);
    if (g.getType() == 5) b = extendBox(g.getCollection(), b);
  }
  return b;
}

// _____________________________________________________________________________
template <typename T>
Point<T> pointAt(const Line<T> l, double at) {
  return pointAtDist(l, at * len(l));
}

// _____________________________________________________________________________
template <typename T>
Point<T> pointAt(const Line<T> l, double at, size_t* lastI, double* totPos) {
  return pointAtDist(l, at * len(l), lastI, totPos);
}

// _____________________________________________________________________________
template <typename T>
Point<T> pointAtDist(const Line<T> l, double atDist) {
  return pointAtDist(l, atDist, nullptr, nullptr);
}

// _____________________________________________________________________________
template <typename T>
Point<T> pointAtDist(const Line<T> l, double atDist, size_t* lastI,
                     double* totPos) {
  if (l.size() == 1) {
    if (lastI) *lastI = 0;
    if (totPos) *totPos = 0;
    return l[1];
  }

  if (atDist > geo::len(l)) atDist = geo::len(l);
  if (atDist < 0) atDist = 0;

  double dist = 0;

  const Point<T>* last = &l[0];

  for (size_t i = 1; i < l.size(); i++) {
    const Point<T>& cur = l[i];
    double d = geo::dist(*last, cur);
    dist += d;

    if (dist > atDist) {
      double p = (d - (dist - atDist));
      if (lastI) *lastI = i - 1;
      if (totPos) *totPos = atDist / util::geo::len(l);
      return interpolate(*last, cur, p / dist);
    }

    last = &l[i];
  }

  if (lastI) *lastI = l.size() - 1;
  if (totPos) *totPos = 1;
  return l.back();
}

// _____________________________________________________________________________
template <typename T>
Point<T> interpolate(const Point<T>& a, const Point<T>& b, double d) {
  double n1 = b.getX() - a.getX();
  double n2 = b.getY() - a.getY();
  return Point<T>(1.0 * a.getX() + (n1 * d), 1.0 * a.getY() + (n2 * d));
}

// _____________________________________________________________________________
template <typename T>
Line<T> orthoLineAtDist(const Line<T>& l, double d, double length) {
  Point<T> avgP = pointAtDist(l, d);

  double angle = angBetween(pointAtDist(l, d - 5), pointAtDist(l, d + 5));

  double angleX1 = 1.0 * avgP.getX() + cos(angle + M_PI_2) * length / 2;
  double angleY1 = 1.0 * avgP.getY() + sin(angle + M_PI_2) * length / 2;

  double angleX2 = 1.0 * avgP.getX() + cos(angle + M_PI_2) * -length / 2;
  double angleY2 = 1.0 * avgP.getY() + sin(angle + M_PI_2) * -length / 2;

  return Line<T>{Point<T>(angleX1, angleY1), Point<T>(angleX2, angleY2)};
}

// _____________________________________________________________________________
template <typename T>
Line<T> segment(const Line<T>& line, double a, double b) {
  if (a > b) {
    double c = a;
    a = b;
    b = c;
  }
  size_t startI, endI;
  auto start = pointAt(line, a, &startI, 0);
  auto end = pointAt(line, b, &endI, 0);

  return segment(line, start, startI, end, endI);
}

// _____________________________________________________________________________
template <typename T>
Line<T> segment(const Line<T>& line, const Point<T>& start, size_t startI,
                const Point<T>& end, size_t endI) {
  Line<T> ret;
  ret.push_back(start);

  if (startI + 1 <= endI) {
    ret.insert(ret.end(), line.begin() + startI + 1, line.begin() + endI + 1);
  }
  ret.push_back(end);

  // find a more performant way to clear the result of above
  ret = util::geo::simplify(ret, 0);

  assert(ret.size());

  return ret;
}

// _____________________________________________________________________________
template <typename T>
Line<T> average(const std::vector<const Line<T>*>& lines) {
  return average(lines, std::vector<double>());
}

// _____________________________________________________________________________
template <typename T>
Line<T> average(const std::vector<const Line<T>*>& lines,
                const std::vector<double>& weights) {
  bool weighted = lines.size() == weights.size();
  double stepSize;

  double longestLength =
      std::numeric_limits<double>::min();  // avoid recalc of length on each
                                           // comparision
  for (auto p : lines) {
    if (len(*p) > longestLength) {
      longestLength = len(*p);
    }
  }

  Line<T> ret;
  double total = 0;

  for (size_t i = 0; i < lines.size(); ++i) {
    if (weighted) {
      total += weights[i];
    } else {
      total += 1;
    }
  }

  stepSize = AVERAGING_STEP / longestLength;
  bool end = false;
  for (double a = 0; !end; a += stepSize) {
    if (a > 1) {
      a = 1;
      end = true;
    }
    double x = 0, y = 0;

    for (size_t i = 0; i < lines.size(); ++i) {
      auto pl = lines[i];
      Point<T> p = pointAt(*pl, a);
      if (weighted) {
        x += p.getX() * weights[i];
        y += p.getY() * weights[i];
      } else {
        x += p.getX();
        y += p.getY();
      }
    }
    ret.push_back(Point<T>(x / total, y / total));
  }

  simplify(ret, 0);

  return ret;
}

// _____________________________________________________________________________
template <typename T>
double area(const Point<T>&) {
  return 0;
}

// _____________________________________________________________________________
template <typename T>
double area(const LineSegment<T>&) {
  return 0;
}

// _____________________________________________________________________________
template <typename T>
double area(const Line<T>&) {
  return 0;
}

// _____________________________________________________________________________
template <typename T>
double area(const Collection<T>& col) {
  double ret = 0;
  for (const auto& g : col) {
    if (g.getType() == 0) ret += area(g.getPoint());
    if (g.getType() == 1) ret += area(g.getLine());
    if (g.getType() == 2) ret += area(g.getPolygon());
    if (g.getType() == 3) ret += area(g.getMultiLine());
    if (g.getType() == 4) ret += area(g.getMultiPolygon());
    if (g.getType() == 5) ret += area(g.getCollection());
  }

  return ret;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
double area(const std::vector<Geometry<T>>& gs) {
  double ret = 0;
  for (const auto& g : gs) {
    ret += area(g);
  }

  return ret;
}

// _____________________________________________________________________________
template <typename T>
double area(const XSortedPolygon<T>& p) {
  return p.area();
}

// _____________________________________________________________________________
template <typename T>
double commonArea(const Box<T>& ba, const Box<T>& bb) {
  double l = std::max(ba.getLowerLeft().getX(), bb.getLowerLeft().getX());
  double r = std::min(ba.getUpperRight().getX(), bb.getUpperRight().getX());
  double b = std::max(ba.getLowerLeft().getY(), bb.getLowerLeft().getY());
  double t = std::min(ba.getUpperRight().getY(), bb.getUpperRight().getY());

  if (l > r || b > t) return 0;
  return (r - l) * (t - b);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
RotatedBox<T> getFullEnvelope(std::vector<Geometry<T>> pol) {
  Point<T> center = centroid(pol);
  Box<T> tmpBox = getBoundingBox(pol);
  double rotateDeg = 0;

  std::vector<Polygon<T>> ml;

  // rotate in 5 deg steps
  for (int i = 1; i < 360; i += 1) {
    pol = rotate(pol, 1, center);
    Polygon<T> hull = convexHull(pol);
    ml.push_back(hull);
    Box<T> e = getBoundingBox(pol);
    if (area(tmpBox) > area(e)) {
      tmpBox = e;
      rotateDeg = i;
    }
  }

  tmpBox = getBoundingBox(ml);

  return RotatedBox<T>(tmpBox, rotateDeg, center);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
RotatedBox<T> getFullEnvelope(const Geometry<T> pol) {
  std::vector<Geometry<T>> mult;
  mult.push_back(pol);
  return getFullEnvelope(mult);
}

// _____________________________________________________________________________
template <typename T>
RotatedBox<T> getOrientedEnvelopeAvg(MultiLine<T> ml) {
  MultiLine<T> orig = ml;
  // get oriented envelope for hull
  RotatedBox<T> rbox = getFullEnvelope(ml);
  Point<T> center = centroid(rbox.getBox());

  ml = rotate(ml, -rbox.getDegree() - 45, center);

  double bestDeg = -45;
  double score = parallelity(rbox.getBox(), ml);

  double i = -45;
  while (i <= 45) {
    ml = rotate(ml, -.5, center);
    double p = parallelity(rbox.getBox(), ml);
    if (parallelity(rbox.getBox(), ml) > score) {
      bestDeg = i;
      score = p;
    }

    i += .5;
  }

  rbox.setDegree(rbox.getDegree() + bestDeg);

  // move the box along 45deg angles from its origin until it fits the ml
  // = until the intersection of its hull and the box is largest
  Polygon<T> p = convexHull(rbox);
  p = rotate(p, -rbox.getDegree(), rbox.getCenter());

  Polygon<T> hull = convexHull(orig);
  hull = rotate(hull, -rbox.getDegree(), rbox.getCenter());

  Box<T> box = getBoundingBox(hull);
  rbox = RotatedBox<T>(box, rbox.getDegree(), rbox.getCenter());

  return rbox;
}

// _____________________________________________________________________________
template <typename T>
double haversine(T lat1, T lon1, T lat2, T lon2) {
  lat1 *= RAD;
  lat2 *= RAD;

  const double dLat = lat2 - lat1;
  const double dLon = (lon2 - lon1) * RAD;

  const double sDLat = sin(dLat / 2);
  const double sDLon = sin(dLon / 2);

  const double a = (sDLat * sDLat) + (sDLon * sDLon) * cos(lat1) * cos(lat2);
  return EQUATORIAL_RAD * 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
}

// _____________________________________________________________________________
template <typename T>
double haversine(const Point<T>& a, const Point<T>& b) {
  return haversine(a.getY(), a.getX(), b.getY(), b.getX());
}

// _____________________________________________________________________________
template <typename T>
Line<T> sparseify(const Line<T>& l, double mind) {
  if (l.size() < 3) return l;

  Line<T> ret;
  ret.push_back(l.front());

  for (size_t i = 1; i < l.size() - 1; i++) {
    double segd = dist(l[i], ret.back());
    if (segd >= mind) ret.push_back(l[i]);
  }

  ret.push_back(l.back());

  return ret;
}

// _____________________________________________________________________________
template <typename T>
Line<T> densify(const Line<T>& l, double d) {
  if (!l.size()) return l;
  if (d <= 0) return l;

  // compute number of required points
  size_t exp = l.size();
  for (size_t i = 1; i < l.size(); i++) {
    exp += dist(l[i - 1], l[i]) / d;
  }

  // shortcut
  if (exp == l.size()) return l;

  Line<T> ret;
  ret.reserve(exp);
  ret.push_back(l.front());

  for (size_t i = 1; i < l.size(); i++) {
    double segd = dist(l[i - 1], l[i]);
    double dx = (l[i].getX() - l[i - 1].getX()) / segd;
    double dy = (l[i].getY() - l[i - 1].getY()) / segd;
    double curd = d;
    while (curd < segd) {
      ret.push_back(Point<T>(1.0 * l[i - 1].getX() + dx * curd,
                             1.0 * l[i - 1].getY() + dy * curd));
      curd += d;
    }

    ret.push_back(l[i]);
  }

  return ret;
}

// _____________________________________________________________________________
template <typename T>
double frechetDistC(size_t i, size_t j, const Line<T>& p, const Line<T>& q,
                    std::vector<float>& ca) {
  // based on Eiter / Mannila
  // http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf

  if (ca[i * q.size() + j] > -1)
    return ca[i * q.size() + j];
  else if (i == 0 && j == 0)
    ca[i * q.size() + j] = dist(p[0], q[0]);
  else if (i > 0 && j == 0)
    ca[i * q.size() + j] =
        std::max(frechetDistC(i - 1, 0, p, q, ca), dist(p[i], q[0]));
  else if (i == 0 && j > 0)
    ca[i * q.size() + j] =
        std::max(frechetDistC(0, j - 1, p, q, ca), dist(p[0], q[j]));
  else if (i > 0 && j > 0)
    ca[i * q.size() + j] =
        std::max(std::min(std::min(frechetDistC(i - 1, j, p, q, ca),
                                   frechetDistC(i - 1, j - 1, p, q, ca)),
                          frechetDistC(i, j - 1, p, q, ca)),
                 dist(p[i], q[j]));
  else
    ca[i * q.size() + j] = std::numeric_limits<float>::infinity();

  return ca[i * q.size() + j];
}

// _____________________________________________________________________________
template <typename T>
double frechetDist(const Line<T>& a, const Line<T>& b, double d) {
  // based on Eiter / Mannila
  // http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf

  const auto& p = densify(a, d);
  const auto& q = densify(b, d);

  std::vector<float> ca(p.size() * q.size(), -1.0);
  double fd = frechetDistC(p.size() - 1, q.size() - 1, p, q, ca);

  return fd;
}

// _____________________________________________________________________________
template <typename T>
double accFrechetDistC(const Line<T>& a, const Line<T>& b, double d) {
  const auto& p = densify(a, d);
  const auto& q = densify(b, d);

  assert(p.size());
  assert(q.size());

  std::vector<float> ca(p.size() * q.size(), 0);

  for (size_t i = 0; i < p.size(); i++)
    ca[i * q.size() + 0] = std::numeric_limits<float>::infinity();
  for (size_t j = 0; j < q.size(); j++)
    ca[j] = std::numeric_limits<float>::infinity();
  ca[0] = 0;

  for (size_t i = 1; i < p.size(); i++) {
    for (size_t j = 1; j < q.size(); j++) {
      float d = util::geo::dist(p[i], q[j]) * util::geo::dist(p[i], p[i - 1]);
      ca[i * q.size() + j] =
          d + std::min(ca[(i - 1) * q.size() + j],
                       std::min(ca[i * q.size() + (j - 1)],
                                ca[(i - 1) * q.size() + (j - 1)]));
    }
  }

  return ca[p.size() * q.size() - 1];
}

// _____________________________________________________________________________
template <typename T>
double frechetDistCHav(size_t i, size_t j, const Line<T>& p, const Line<T>& q,
                       std::vector<float>& ca) {
  // based on Eiter / Mannila
  // http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf

  if (ca[i * q.size() + j] > -1)
    return ca[i * q.size() + j];
  else if (i == 0 && j == 0)
    ca[i * q.size() + j] = haversine(p[0], q[0]);
  else if (i > 0 && j == 0)
    ca[i * q.size() + j] =
        std::max(frechetDistCHav(i - 1, 0, p, q, ca), haversine(p[i], q[0]));
  else if (i == 0 && j > 0)
    ca[i * q.size() + j] =
        std::max(frechetDistCHav(0, j - 1, p, q, ca), haversine(p[0], q[j]));
  else if (i > 0 && j > 0)
    ca[i * q.size() + j] =
        std::max(std::min(std::min(frechetDistCHav(i - 1, j, p, q, ca),
                                   frechetDistCHav(i - 1, j - 1, p, q, ca)),
                          frechetDistCHav(i, j - 1, p, q, ca)),
                 haversine(p[i], q[j]));
  else
    ca[i * q.size() + j] = std::numeric_limits<float>::infinity();

  return ca[i * q.size() + j];
}

// _____________________________________________________________________________
template <typename T>
double frechetDistHav(const Line<T>& a, const Line<T>& b, double d) {
  // based on Eiter / Mannila
  // http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf

  const auto& p = densify(a, d);
  const auto& q = densify(b, d);

  std::vector<float> ca(p.size() * q.size(), -1.0);
  double fd = frechetDistCHav(p.size() - 1, q.size() - 1, p, q, ca);

  return fd;
}

// _____________________________________________________________________________
template <typename T>
double accFrechetDistCHav(const Line<T>& a, const Line<T>& b, double d) {
  const auto& p = densify(a, d);
  const auto& q = densify(b, d);

  assert(p.size());
  assert(q.size());

  std::vector<float> ca(p.size() * q.size(), 0);

  for (size_t i = 0; i < p.size(); i++)
    ca[i * q.size() + 0] = std::numeric_limits<float>::infinity();
  for (size_t j = 0; j < q.size(); j++)
    ca[j] = std::numeric_limits<float>::infinity();
  ca[0] = 0;

  for (size_t i = 1; i < p.size(); i++) {
    for (size_t j = 1; j < q.size(); j++) {
      float d = util::geo::haversine(p[i], q[j]) *
                util::geo::haversine(p[i], p[i - 1]);
      ca[i * q.size() + j] =
          d + std::min(ca[(i - 1) * q.size() + j],
                       std::min(ca[i * q.size() + (j - 1)],
                                ca[(i - 1) * q.size() + (j - 1)]));
    }
  }

  return ca[p.size() * q.size() - 1];
}

// _____________________________________________________________________________
template <typename T>
Point<T> latLngToWebMerc(double lat, double lng) {
  // clamp latitudes to web mercator range
  if (lat > 85.05112878) lat = 85.05112878;
  if (lat < -85.05112878) lat = -85.05112878;

  // restrict longitude to [-180, 180], wrap around at boundaries
  lng = std::remainder(lng, 360.0);

  double x = EQUATORIAL_RAD * lng * 0.017453292519943295;
  double sina = sin(lat * 0.017453292519943295);

  double y = 3189068.5 * log((1.0 + sina) / (1.0 - sina));
  return Point<T>(x, y);
}

// _____________________________________________________________________________
template <typename T>
// TODO: rename to lngLat
Point<T> latLngToWebMerc(Point<T> lngLat) {
  return latLngToWebMerc<T>(lngLat.getY(), lngLat.getX());
}

// _____________________________________________________________________________
template <typename T>
Point<T> webMercToLatLng(double x, double y) {
  const double lat =
      (1.5707963267948966 - (2.0 * atan(exp(-y / EQUATORIAL_RAD)))) * IRAD;
  const double lon = x / 111319.4907932735677;
  return Point<T>(lon, lat);
}

// _____________________________________________________________________________
template <typename T>
double webMercMeterDist(const Point<T>& a, const Point<T>& b) {
  const auto llA = webMercToLatLng<T>(a.getX(), a.getY());
  const auto llB = webMercToLatLng<T>(b.getX(), b.getY());
  return haversine(llA.getY(), llA.getX(), llB.getY(), llB.getX());
}

// _____________________________________________________________________________
template <typename G1, typename G2>
double webMercMeterDist(const G1& a, const G2& b) {
  // euclidean distance on web mercator is in meters on equator,
  // and proportional to cos(lat) in both y directions

  // this is just an approximation!

  auto pa = centroid(a);
  auto pb = centroid(b);

  double fA = webMercDistFactor(pa);
  double fB = webMercDistFactor(pb);

  return util::geo::dist(a, b) * cos((fA + fB) / 2.0);
}

// _____________________________________________________________________________
template <typename T>
double webMercLen(const Line<T>& g) {
  double ret = 0;
  for (size_t i = 1; i < g.size(); i++) ret += webMercMeterDist(g[i - 1], g[i]);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
double latLngLen(const Line<T>& g) {
  double ret = 0;
  for (size_t i = 1; i < g.size(); i++) ret += haversine(g[i - 1], g[i]);
  return ret;
}

// _____________________________________________________________________________
template <typename G>
double webMercDistFactor(const G& a) {
  // euclidean distance on web mercator is in meters on equator,
  // and proportional to cos(lat) in both y directions
  double et = exp(a.getY() / EQUATORIAL_RAD);
  return 2 * et / (et * et + 1);
}

// _____________________________________________________________________________
template <typename G>
double latLngDistFactor(const G& a) {
  // euclidean distance on web mercator is in meters on equator,
  // and proportional to cos(lat) in both y directions

  return cos(a.getY() * RAD);
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const std::vector<XSortedTuple<T>>& ls1,
                  const std::vector<XSortedTuple<T>>& ls2, T maxSegLenA,
                  T maxSegLenB, const Box<T>& boxA, const Box<T>& boxB,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  if (util::geo::dist(boxA, boxB) > maxEuclideanDist) {
    return std::numeric_limits<double>::max();
  }

  // always ensure that ls2 is smaller, because ls2 is padded
  if (ls1.size() < ls2.size())
    return withinDist(ls2, ls1, maxSegLenB, maxSegLenA, boxB, boxA, maxDist,
                      paddingFunc, maxEuclideanDist, distFunc);

  if (ls1.size() == 0 || ls2.size() == 0)
    return std::numeric_limits<double>::max();

  auto probeDists =
      probeDistanceUpperBound(100, ls1, ls2, boxA, boxB, maxSegLenB, maxSegLenA,
                              maxEuclideanDist, distFunc);

  // if we already have the correct distance from probing, return
  if (std::get<2>(probeDists)) {
    return std::get<0>(probeDists);
  }

  double minDist = std::get<0>(probeDists);
  double euclideanDistUpperBound = std::get<1>(probeDists);

  if (util::geo::dist(boxA, boxB) > euclideanDistUpperBound) {
    return std::numeric_limits<double>::max();
  }

  // minimum euclidean X dist based on bounding boxes for better y padding
  double euclideanXDistLowerBound =
      dist(LineSegment<T>{Point<T>{boxA.getLowerLeft().getX(), 0},
                          Point<T>{boxA.getUpperRight().getX(), 0}},
           LineSegment<T>{Point<T>{boxB.getLowerLeft().getX(), 0},
                          Point<T>{boxB.getUpperRight().getX(), 0}});
  double euclideanYDistLowerBound =
      dist(LineSegment<T>{Point<T>{0, boxA.getLowerLeft().getY()},
                          Point<T>{0, boxA.getUpperRight().getY()}},
           LineSegment<T>{Point<T>{0, boxB.getLowerLeft().getY()},
                          Point<T>{0, boxB.getUpperRight().getY()}});

  size_t i = 0;  // position in ls1
  size_t j = 0;  // position in ls2
  size_t k = 0;  // position in OUT ls2
  size_t ls2OutSize = 0;

  double padding = paddingFunc(euclideanDistUpperBound,
                               std::min(minDist, maxDist), boxA, boxB) *
                   1.0;

  T xPadding = std::min(
      std::numeric_limits<T>::max() * 1.0,
      (sqrt(std::max(0.0, padding * padding - euclideanYDistLowerBound *
                                                  euclideanYDistLowerBound))));
  T yPadding = std::min(
      std::numeric_limits<T>::max() * 1.0,
      (sqrt(std::max(0.0, padding * padding - euclideanXDistLowerBound *
                                                  euclideanXDistLowerBound))));

  // skip irrelevant parts in ls1
  if (maxSegLenA < std::numeric_limits<T>::max() &&
      ls1[i].p.getX() <
          boundedSub(boundedSub(ls2[j].p.getX(), maxSegLenA), xPadding)) {
    i = std::lower_bound(
            ls1.begin() + i, ls1.end(),
            XSortedTuple<T>{
                {boundedSub(boundedSub(ls2[j].p.getX(), maxSegLenA), xPadding),
                 0},
                false}) -
        ls1.begin();
  }

  // skip irrelevant parts in ls2
  if (maxSegLenB < std::numeric_limits<T>::max() &&
      ls2[j].p.getX() <
          boundedSub(boundedSub(ls1[i].p.getX(), maxSegLenB), xPadding)) {
    j = std::lower_bound(
            ls2.begin() + j, ls2.end(),
            XSortedTuple<T>{
                {boundedSub(boundedSub(ls1[i].p.getX(), maxSegLenB), xPadding),
                 0},
                false}) -
        ls2.begin();
  }

  while (i < ls1.size() && ls1[i].seg().second.getX() <
                               boundedSub(boxB.getLowerLeft().getX(), xPadding))
    i++;
  while (j < ls2.size() && ls2[j].seg().second.getX() <
                               boundedSub(boxA.getLowerLeft().getX(), xPadding))
    j++;

  // segments active by their padded bounding box
  util::geo::IntervalIdx<T, LineSegment<T>> activesA, activesB;

  // return buffer for interval idx
  std::vector<util::geo::IntervalVal<T, LineSegment<T>>> segs;

  // ls2 is always padded!
  while (i < ls1.size() || j < ls2.size() || k < ls2OutSize) {
    T ls1X = i < ls1.size() ? ls1[i].p.getX() : std::numeric_limits<T>::max();
    bool ls1Out = i < ls1.size() ? ls1[i].out() : true;
    T ls2X = j < ls2.size() ? boundedSub(ls2[j].p.getX(), xPadding)
                            : std::numeric_limits<T>::max();
    bool ls2Out = j < ls2.size() ? ls2[j].out() : true;

    T ls2OutX = k < ls2OutSize ? boundedAdd(ls2[k].p.getX(), xPadding)
                               : std::numeric_limits<T>::max();

    if ((ls1X < ls2X || (ls1X == ls2X && (!ls1Out || ls2Out))) &&
        (ls1X < ls2OutX || (ls1X == ls2OutX))) {
      // advance ls1
      const auto& ls1seg = ls1[i].seg();

      // we are past ls2
      if (ls1X > boundedAdd(ls2.back().p.getX(), xPadding)) break;

      // ignore segments out of the X range
      if (ls1seg.second.getX() <
          boundedSub(boxB.getLowerLeft().getX(), xPadding)) {
        i++;
        continue;
      }

      // ignore segments out of the Y range
      if (ls1seg.first.getY() < boxB.getLowerLeft().getY() - yPadding &&
          ls1seg.second.getY() < boxB.getLowerLeft().getY() - yPadding) {
        i++;
        continue;
      }

      if (ls1seg.first.getY() > boxB.getUpperRight().getY() + yPadding &&
          ls1seg.second.getY() > boxB.getUpperRight().getY() + yPadding) {
        i++;
        continue;
      }

      if (withinDist(ls1seg, boxB, distFunc, minDist) > minDist) {
        i++;
        continue;
      }

      const auto& box = getBoundingBox(ls1seg);

      if (!ls1Out) {
        // IN event
        activesA.insert({box.getLowerLeft().getY(), box.getUpperRight().getY()},
                        ls1seg);
      } else {
        // OUT event
        activesA.erase({box.getLowerLeft().getY(), box.getUpperRight().getY()},
                       ls1seg);

        if (processActives(activesB, ls1seg, euclideanDistUpperBound, minDist,
                           maxDist, padding, xPadding, yPadding,
                           euclideanXDistLowerBound, euclideanYDistLowerBound,
                           box, boxB, boxA, segs, paddingFunc, distFunc))
          return minDist;
      }

      i++;
    } else if ((ls2X < ls1X || (ls2X == ls1X && (!ls2Out || ls1Out))) &&
               (ls2X <= ls2OutX)) {
      // advance ls2
      const auto& ls2seg = ls2[j].seg();

      // we are past ls1, so simply break
      if (ls2X > ls1.back().p.getX()) break;

      // ignore segments out of the X range
      if (boundedAdd(ls2seg.second.getX(), xPadding) <
          boxA.getLowerLeft().getX()) {
        j++;
        continue;
      }

      // ignore segments out of the Y range
      if (boundedAdd(ls2seg.first.getY(), yPadding) <
              boxA.getLowerLeft().getY() &&
          boundedAdd(ls2seg.second.getY(), yPadding) <
              boxA.getLowerLeft().getY()) {
        j++;
        continue;
      }

      if (boundedSub(ls2seg.first.getY(), yPadding) >
              boxA.getUpperRight().getY() &&
          boundedSub(ls2seg.second.getY(), yPadding) >
              boxA.getUpperRight().getY()) {
        j++;
        continue;
      }

      const auto& box = getBoundingBox(ls2seg);
      if (withinDist(ls2seg, boxA, distFunc, minDist) > minDist) {
        j++;
        continue;
      }

      if (!ls2Out) {
        // IN event
        activesB.insert({box.getLowerLeft().getY(), box.getUpperRight().getY()},
                        ls2seg);
      } else {
        // OUT event, deferred because of padding
        ls2OutSize = j;
      }

      j++;
    } else if ((ls2OutX < ls1X) && (ls2OutX < ls2X)) {
      // advance ls2 out
      const auto& ls2OutSeg = ls2[k].seg();

      // ignore segments out of the X range
      if (boundedAdd(ls2OutSeg.second.getX(), xPadding) <
          boxA.getLowerLeft().getX()) {
        k++;
        while (k < ls2OutSize && !ls2[k].out()) k++;
        continue;
      }

      // ignore segments out of the Y range
      if (boundedSub(ls2OutSeg.first.getY(), yPadding) <
              boxA.getLowerLeft().getY() &&
          boundedAdd(ls2OutSeg.second.getY(), yPadding) <
              boxA.getLowerLeft().getY()) {
        k++;
        while (k < ls2OutSize && !ls2[k].out()) k++;
        continue;
      }

      if (boundedSub(ls2OutSeg.first.getY(), yPadding) >
              boxA.getUpperRight().getY() &&
          boundedAdd(ls2OutSeg.second.getY(), yPadding) >
              boxA.getUpperRight().getY()) {
        k++;
        while (k < ls2OutSize && !ls2[k].out()) k++;
        continue;
      }

      const auto& box = getBoundingBox(ls2OutSeg);

      activesB.erase({box.getLowerLeft().getY(), box.getUpperRight().getY()},
                     ls2OutSeg);

      if (processActives(activesA, ls2OutSeg, euclideanDistUpperBound, minDist,
                         maxDist, padding, xPadding, yPadding,
                         euclideanXDistLowerBound, euclideanYDistLowerBound,
                         box, boxA, boxB, segs, paddingFunc, distFunc))
        return minDist;

      // advance to next OUT
      k++;
      while (k < ls2OutSize && !ls2[k].out()) k++;
    }
  }

  return minDist;
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
std::pair<double, std::pair<bool, bool>> withinDist(
    const XSortedRing<T>& p1, const XSortedRing<T>& p2, double maxDist,
    PF&& paddingFunc, double maxEuclideanDist, DF&& distFunc) {
  if (p1.rawRing().size() == 1) {
    auto a = withinDist(p1.rawRing().front().p, p2, maxDist, paddingFunc,
                        maxEuclideanDist, distFunc);
    return {a.first, {false, a.second}};
  }
  if (p2.rawRing().size() == 1) {
    auto a = withinDist(p2.rawRing().front().p, p1, maxDist, paddingFunc,
                        maxEuclideanDist, distFunc);
    return {a.first, {a.second, false}};
  }

  if (p1.rawRing().size() == 0)
    return {std::numeric_limits<double>::max(), {false, false}};
  if (p2.rawRing().size() == 0)
    return {std::numeric_limits<double>::max(), {false, false}};

  const Box<T> boxA = p1.boundingBox();
  const Box<T> boxB = p2.boundingBox();

  double ringDist = util::geo::withinDist(
      p1.rawRing(), p2.rawRing(), p1.getMaxSegLen(), p2.getMaxSegLen(), boxA,
      boxB, maxDist, paddingFunc, maxEuclideanDist, distFunc);

  if (ringDist == 0) return {0, {false, false}};
  if (util::geo::ringContains(p1.rawRing().front().seg().second, p2, 0)
          .second) {
    return {ringDist, {false, true}};
  }
  if (util::geo::ringContains(p2.rawRing().front().seg().second, p1, 0)
          .second) {
    return {ringDist, {true, false}};
  }

  return {ringDist, {false, false}};
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
std::pair<double, bool> withinDist(const XSortedLine<T>& ls1,
                                   const XSortedRing<T>& p2, double maxDist,
                                   PF&& paddingFunc, double maxEuclideanDist,
                                   DF&& distFunc) {
  if (ls1.rawLine().size() == 1) {
    return withinDist(ls1.rawLine().front().p, p2, maxDist, paddingFunc,
                      maxEuclideanDist, distFunc);
  }

  if (p2.rawRing().size() == 1) {
    return {withinDist(p2.rawRing().front().p, ls1, maxDist, paddingFunc,
                       maxEuclideanDist, distFunc),
            false};
  }

  if (ls1.rawLine().size() == 0)
    return {std::numeric_limits<double>::max(), false};
  if (p2.rawRing().size() == 0)
    return {std::numeric_limits<double>::max(), false};

  double ringDist = util::geo::withinDist(
      ls1.rawLine(), p2.rawRing(), ls1.getMaxSegLen(), p2.getMaxSegLen(),
      ls1.boundingBox(), p2.boundingBox(), maxDist, paddingFunc,
      maxEuclideanDist, distFunc);

  if (ringDist == 0) return {ringDist, false};
  if (util::geo::ringContains(ls1.rawLine().front().seg().second, p2, 0).second)
    return {ringDist, true};

  return {ringDist, false};
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const XSortedLine<T>& a, const XSortedPolygon<T>& b,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  if (b.getInners().size() == 0 &&
      util::geo::ringContains(a.rawLine().front().seg().second, b.getOuter(), 0)
          .second)
    return 0;
  auto r = withinDist(a, b.getOuter(), maxDist, paddingFunc, maxEuclideanDist,
                      distFunc);
  if (!r.second) return r.first;

  // also check inner rings
  if (r.second && b.getInners().size()) {
    size_t i = 0;

    if (b.getInnerMaxSegLen() < std::numeric_limits<T>::max()) {
      i = std::lower_bound(
              b.getInnerBoxIdx().begin(), b.getInnerBoxIdx().end(),
              std::pair<T, size_t>{
                  boundedSub(a.rawLine().front().seg().first.getX(),
                             b.getInnerMaxSegLen()),
                  0}) -
          b.getInnerBoxIdx().begin();
    }

    for (; i < b.getInners().size(); i++) {
      if (b.getInners()[i].rawRing().size() < 2) continue;
      if (b.getInnerBoxes()[i].getLowerLeft().getX() >
          a.rawLine().back().seg().second.getX())
        break;
      if (!util::geo::contains(a.boundingBox(), b.getInnerBoxes()[i])) continue;

      auto r2 = withinDist(a, b.getInners()[i], maxDist, paddingFunc,
                           maxEuclideanDist, distFunc);

      // if we are contained in the inner ring, directly return the distance to
      // it
      if (r2.second) return r2.first;

      if (r2.first == 0) return 0;
    }
  }

  return r.second ? 0 : r.first;
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const XSortedPolygon<T>& p1, const XSortedPolygon<T>& p2,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  if (p1.getOuter().rawRing().size() == 0)
    return nextafter(maxDist, std::numeric_limits<double>::infinity());
  if (p2.getOuter().rawRing().size() == 0)
    return nextafter(maxDist, std::numeric_limits<double>::infinity());

  if (withinDist(p1.boundingBox(), p2.boundingBox(), distFunc, maxDist) >
      maxDist)
    return nextafter(maxDist, std::numeric_limits<double>::infinity());

  if (p1.getInners().size() == 0 &&
      util::geo::ringContains(p2.getOuter().rawRing().front().seg().second,
                              p1.getOuter(), 0)
          .second)
    return 0;
  if (p2.getInners().size() == 0 &&
      util::geo::ringContains(p1.getOuter().rawRing().front().seg().second,
                              p2.getOuter(), 0)
          .second)
    return 0;

  auto outerR = withinDist(p1.getOuter(), p2.getOuter(), maxDist, paddingFunc,
                           maxEuclideanDist, distFunc);

  // if we are not contained, directly return the distance (0 if we intersect)
  if (!outerR.second.first && !outerR.second.second) {
    return outerR.first;
  }

  // rationale: if we are contained in the outer ring, the distance can only be
  // nonzero if we are inside an inner ring. But the distance to this inner ring
  // cannot be greater than the distance to the outer ring, because otherwise
  // the inner ring would lie outside the outer ring, which we do not support
  maxDist = std::min(maxDist, outerR.first);

  const Box<T>& boxA = p1.boundingBox();
  const Box<T>& boxB = p2.boundingBox();

  // if p1 is contained in p2, check if p1 is in an inner ring of p2
  if (outerR.second.second && p2.getInners().size()) {
    size_t i = 0;

    if (p2.getInnerMaxSegLen() < std::numeric_limits<T>::max()) {
      i = std::lower_bound(
              p2.getInnerBoxIdx().begin(), p2.getInnerBoxIdx().end(),
              std::pair<T, size_t>{
                  boundedSub(p1.getOuter().rawRing().front().seg().first.getX(),
                             p2.getInnerMaxSegLen()),
                  0}) -
          p2.getInnerBoxIdx().begin();
    }

    for (; i < p2.getInners().size(); i++) {
      if (p2.getInnerBoxes()[i].getLowerLeft().getX() >
          p1.getOuter().rawRing().back().seg().second.getX())
        break;
      if (p2.getInners()[i].rawRing().size() < 2) continue;
      if (!util::geo::contains(boxA, p2.getInnerBoxes()[i])) continue;

      double maxLocalEuclideanDist = std::min(
          maxEuclideanDist,
          1.0 * std::min(p2.getInnerBoxes()[i].getUpperRight().getX() -
                             p2.getInnerBoxes()[i].getLowerLeft().getX(),
                         p2.getInnerBoxes()[i].getUpperRight().getY() -
                             p2.getInnerBoxes()[i].getLowerLeft().getY()));
      auto r2 = withinDist(p1.getOuter(), p2.getInners()[i], maxDist,
                           paddingFunc, maxLocalEuclideanDist, distFunc);

      // if we are contained in the inner ring, directly return the distance to
      // it
      if (r2.second.second) return r2.first;

      // if the distance is 0, also (we intersect)
      if (r2.first == 0) return 0;
    }
  }

  // if p2 is contained in p1, check if p2 is in an inner ring of p1
  if (outerR.second.first && p1.getInners().size()) {
    size_t i = 0;

    if (p1.getInnerMaxSegLen() < std::numeric_limits<T>::max()) {
      i = std::lower_bound(
              p1.getInnerBoxIdx().begin(), p1.getInnerBoxIdx().end(),
              std::pair<T, size_t>{
                  boundedSub(p2.getOuter().rawRing().front().seg().first.getX(),
                             p1.getInnerMaxSegLen()),
                  0}) -
          p1.getInnerBoxIdx().begin();
    }

    for (; i < p1.getInners().size(); i++) {
      if (p1.getInnerBoxes()[i].getLowerLeft().getX() >
          p2.getOuter().rawRing().back().seg().second.getX())
        break;
      if (p1.getInners()[i].rawRing().size() < 2) continue;
      if (!util::geo::contains(boxB, p1.getInnerBoxes()[i])) continue;

      double maxLocalEuclideanDist = std::min(
          maxEuclideanDist,
          1.0 * std::min(p1.getInnerBoxes()[i].getUpperRight().getX() -
                             p1.getInnerBoxes()[i].getLowerLeft().getX(),
                         p1.getInnerBoxes()[i].getUpperRight().getY() -
                             p1.getInnerBoxes()[i].getLowerLeft().getY()));
      auto r2 = withinDist(p2.getOuter(), p1.getInners()[i], maxDist,
                           paddingFunc, maxLocalEuclideanDist, distFunc);

      // if we are contained in the inner ring, directly return the distance to
      // it
      if (r2.second.second) return r2.first;

      // if the distance is 0, also (we intersect)
      if (r2.first == 0) return 0;
    }
  }

  return (outerR.second.first || outerR.second.second) ? 0 : outerR.first;
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
double withinDist(const XSortedLine<T>& ls1, const XSortedLine<T>& ls2,
                  double maxDist, PF&& paddingFunc, double maxEuclideanDist,
                  DF&& distFunc) {
  return withinDist(ls1.rawLine(), ls2.rawLine(), ls1.getMaxSegLen(),
                    ls2.getMaxSegLen(), ls1.boundingBox(), ls2.boundingBox(),
                    maxDist, paddingFunc, maxEuclideanDist, distFunc);
}

// _____________________________________________________________________________
template <typename T, typename DF>
std::tuple<double, double, bool> probeDistanceUpperBound(
    size_t maxComps, const std::vector<XSortedTuple<T>>& ls1,
    const std::vector<XSortedTuple<T>>& ls2, const Box<T>& boxA,
    const Box<T>& boxB, T maxSegLenA, T maxSegLenB, double euclideanUpperBound,
    DF&& distFunc) {
  double upperBound = std::numeric_limits<double>::infinity();
  double eucSquared = euclideanUpperBound * euclideanUpperBound;

  size_t samplesA = std::max(
      static_cast<size_t>(1),
      std::min(ls1.size(), static_cast<size_t>(sqrt(
                               maxComps * 1.0 *
                               ((ls1.size() * 1.0) / (ls2.size() * 1.0))))));
  size_t samplesB = std::max(
      static_cast<size_t>(1),
      std::min(ls2.size(),
               static_cast<size_t>((maxComps * 1.0) / (samplesA * 1.0))));
  samplesA = std::max(
      static_cast<size_t>(1),
      std::min(ls1.size(),
               static_cast<size_t>((maxComps * 1.0) / (samplesB * 1.0))));

  size_t stepA = ls1.size() / samplesA;
  size_t stepB = ls2.size() / samplesB;

  for (size_t i = 0; i < ls1.size(); i += stepA) {
    if (ls1[i].out()) continue;
    if (ls1[i].p.getX() + maxSegLenA + euclideanUpperBound <
        boxB.getLowerLeft().getX())
      continue;
    if (ls1[i].p.getX() - euclideanUpperBound > boxB.getUpperRight().getX())
      break;
    const auto& ls1seg = ls1[i].seg();
    if (distSquared(ls1seg, boxB) > eucSquared) continue;

    for (size_t j = 0; j < ls2.size(); j += stepB) {
      if (ls2[j].out()) continue;
      if (ls2[j].p.getX() + maxSegLenB + euclideanUpperBound <
          boxA.getLowerLeft().getX())
        continue;
      if (ls2[j].p.getX() - euclideanUpperBound > boxA.getUpperRight().getX())
        break;
      const auto& ls2seg = ls2[j].seg();
      double euD = distSquared(ls1seg, ls2seg);

      // early abort
      if (euD == 0) return {0, 0, true};

      if (euD > eucSquared) continue;

      eucSquared = euD;
      euclideanUpperBound = sqrt(eucSquared);
      double d = dist(ls1seg, ls2seg, distFunc);
      if (d < upperBound) upperBound = d;
    }
  }

  return {upperBound, euclideanUpperBound, stepA == 1 && stepB == 1};
}

// _____________________________________________________________________________
template <typename T, typename PF, typename DF>
inline bool processActives(util::geo::IntervalIdx<T, LineSegment<T>>& actives,
                           const LineSegment<T>& curSeg,
                           double& euclideanDistUpperBound, double& minDist,
                           const double maxDist, double& padding, T& xPadding,
                           T& yPadding, const double euclideanXDistLowerBound,
                           const double euclideanYDistLowerBound,
                           const Box<T>& box, const util::geo::Box<T>& otherBox,
                           const Box<T>& thisBox,
                           std::vector<IntervalVal<T, LineSegment<T>>>& segs,
                           PF&& paddingFunc, DF&& distFunc) {
  const T aLo = std::min(curSeg.first.getX(), curSeg.second.getX());
  const T aHi = std::max(curSeg.first.getX(), curSeg.second.getX());
  const T bLo = otherBox.getLowerLeft().getX();
  const T bHi = otherBox.getUpperRight().getX();
  double locMinEuclideanXDist = aLo > bHi
                                    ? boundedSub(aLo, bHi)
                                    : (bLo > aHi ? boundedSub(bLo, aHi) : 0.0);

  T locYPadding = std::min(
      std::numeric_limits<T>::max() * 1.0,
      (sqrt(std::max(0.0, padding * padding -
                              locMinEuclideanXDist * locMinEuclideanXDist))));

  actives.overlap_find_all(
      {boundedSub(box.getLowerLeft().getY(), locYPadding),
       boundedAdd(box.getUpperRight().getY(), locYPadding)},
      segs);

  bool minDistUpdated = false;

  for (const auto& seg : segs) {
    if (minDistUpdated) {
      // clear activesB along the way
      if (withinDist(seg.v, thisBox, distFunc, minDist) > minDist) {
        actives.erase({seg.l, seg.r}, seg.v);
        continue;
      }
      minDistUpdated = false;
    }

    bool localUpdated = false;

    double euclideanDist;
    double segD = util::geo::dist(curSeg, seg.v, padding, distFunc,
                                  std::min(minDist, maxDist), euclideanDist);
    if (euclideanDist <= padding) {
      if (segD < minDist) {
        minDist = segD;
        minDistUpdated = true;
        localUpdated = true;
      }

      if (minDist == 0) return true;

      if (euclideanDist < euclideanDistUpperBound) {
        euclideanDistUpperBound = euclideanDist;
        localUpdated = true;
      }

      if (localUpdated) {
        padding = paddingFunc(euclideanDistUpperBound,
                              std::min(minDist, maxDist), otherBox, thisBox) *
                  1.0;
        xPadding =
            std::min(std::numeric_limits<T>::max() * 1.0,
                     (sqrt(std::max(0.0, padding * padding -
                                             euclideanYDistLowerBound *
                                                 euclideanYDistLowerBound))));
        yPadding =
            std::min(std::numeric_limits<T>::max() * 1.0,
                     (sqrt(std::max(0.0, padding * padding -
                                             euclideanXDistLowerBound *
                                                 euclideanXDistLowerBound))));
      }
    }
  }

  return false;
}

// _____________________________________________________________________________
template <typename T>
double meterDistLocalSearchPadding(double euclideanDistanceUpperBound,
                                   double distanceUpperBound,
                                   const util::geo::Box<T>& boxA,
                                   const util::geo::Box<T>& boxB) {
  auto minScaleALow =
      getMinMaxLocalScaleFactors(util::geo::getBoundingBox(boxA.getLowerLeft()),
                                 boxB, distanceUpperBound)
          .second;
  auto minScaleAUp = getMinMaxLocalScaleFactors(
                         util::geo::getBoundingBox(boxA.getUpperRight()), boxB,
                         distanceUpperBound)
                         .second;
  auto minScaleBLow =
      getMinMaxLocalScaleFactors(util::geo::getBoundingBox(boxB.getLowerLeft()),
                                 boxA, distanceUpperBound)
          .second;
  auto minScaleBUp = getMinMaxLocalScaleFactors(
                         util::geo::getBoundingBox(boxB.getUpperRight()), boxA,
                         distanceUpperBound)
                         .second;

  double min = std::min(std::min(minScaleALow, minScaleAUp),
                        std::min(minScaleBLow, minScaleBUp));
  auto a = getMinMaxLocalScaleFactors(boxA, boxB, distanceUpperBound);
  double max = a.second;

  double factorNew2 = max / min;

  double minEuclideanXDist =
      util::geo::dist(LineSegment<T>{Point<T>{boxA.getLowerLeft().getX(), 0},
                                     Point<T>{boxA.getUpperRight().getX(), 0}},
                      LineSegment<T>{Point<T>{boxB.getLowerLeft().getX(), 0},
                                     Point<T>{boxB.getUpperRight().getX(), 0}});
  double minEuclideanYDist =
      util::geo::dist(LineSegment<T>{Point<T>{0, boxA.getLowerLeft().getY()},
                                     Point<T>{0, boxA.getUpperRight().getY()}},
                      LineSegment<T>{Point<T>{0, boxB.getLowerLeft().getY()},
                                     Point<T>{0, boxB.getUpperRight().getY()}});

  double padding = factorNew2 * euclideanDistanceUpperBound;
  auto xPadding = (sqrt(std::max(
      0.0, padding * padding - minEuclideanYDist * minEuclideanYDist)));
  auto yPadding = (sqrt(std::max(
      0.0, padding * padding - minEuclideanXDist * minEuclideanXDist)));

  auto paddedA = util::geo::pad(boxA, xPadding, yPadding);

  auto boxBStar = util::geo::intersection(paddedA, boxB);

  double min2 = std::numeric_limits<double>::infinity();
  double max2 = 0;

  std::vector<Point<T>> cornerA = {boxA.getLowerLeft(), boxA.getLowerRight(),
                                   boxA.getUpperRight(), boxA.getUpperLeft()};
  std::vector<Point<T>> cornerB = {
      boxBStar.getLowerLeft(), boxBStar.getLowerRight(),
      boxBStar.getUpperRight(), boxBStar.getUpperLeft()};

  for (size_t i = 0; i < cornerA.size(); i++) {
    for (size_t j = 0; j < cornerB.size(); j++) {
      double eucD = util::geo::dist(cornerA[i], cornerB[j]);
      double mD = haversine(cornerA[i], cornerB[j]);
      if (mD / eucD < min2) min2 = mD / eucD;
      if (mD / eucD > max2) max2 = mD / eucD;
    }
  }

  double factorNew3 = max2 / min2;

  if (factorNew2 < factorNew3) return factorNew2 * euclideanDistanceUpperBound;

  return factorNew3 * euclideanDistanceUpperBound;
}

// _____________________________________________________________________________
template <typename T>
std::vector<Point<T>> fill(const Polygon<T>& p, double d,
                           const Box<T>& bounds) {
  std::vector<Point<T>> ret;

  // y bounds
  auto minY = bounds.getLowerLeft().getY();
  auto maxY = bounds.getUpperRight().getY();

  std::vector<double> xs;

  std::vector<const Ring<T>*> rings;
  rings.reserve(p.getInners().size() + 1);
  rings.push_back(&p.getOuter());
  for (const auto& inner : p.getInners()) rings.push_back(&inner);

  // scanline for each y, in res steps
  for (double y = minY; y <= maxY; y += d) {
    xs.clear();

    for (const auto* ring : rings) {
      size_t from = ring->size() - 1;
      size_t to = 0;

      // for each segment [from, to]...
      while (to < ring->size()) {
        int yFr = (*ring)[from].getY();
        int yTo = (*ring)[to].getY();

        // .. check if we intersect the y-scanline, this is asymmetric (we
        // intersect if yFr touches y (then yTo > y and yFr !> y) but not if
        // yTo touches y (then yTo !> y and yFr !> y) to avoid counting
        // points were line segments touch twice
        if ((yTo > y) != (yFr > y)) {
          double t = (y - yFr) / static_cast<double>(yTo - yFr);

          // store the x intersection
          double xIsect = (*ring)[from].getX() +
                          t * ((*ring)[to].getX() - (*ring)[from].getX());
          xs.push_back(xIsect);
        }

        from = to;
        to++;
      }
    }

    // if we have less than 2 intersections, the scanline didnt cut into the
    // itnerior
    if (xs.size() < 2) continue;

    // now sort the intersection Xs, they are now pairs of IN,OUT events
    std::sort(xs.begin(), xs.end());

    // step over these pairs and fill in between in steps of AREA_FILL_RES
    // but keep inside the visible bounds
    for (size_t k = 0; k < xs.size() - 1; k += 2) {
      T xFr = std::max(static_cast<T>(xs[k]), bounds.getLowerLeft().getX());
      T xTo =
          std::min(static_cast<T>(xs[k + 1]), bounds.getUpperRight().getX());

      for (T x = xFr; x <= xTo; x += d) ret.push_back({x, y});
    }
  }

  return ret;
}

// _____________________________________________________________________________
template <typename T>
std::pair<double, double> getMinMaxLocalScaleFactors(
    const Box<T>& boxA, const Box<T>& boxB, double distanceUpperBound) {
  auto withinBox = util::geo::extendBox(boxA, boxB);

  // convert distanceUpperBound (meters) to maximum latitude padding (degrees)
  // we have to "pad" each box by -dy and dy because the distance path could
  // be within that padded box, and thus the distortions have to be computed
  // based on that path
  double dLat =
      distanceUpperBound / util::geo::MIN_METERS_PER_LAT_RAD * util::geo::IRAD;

  // the lower extrema of the bounding boxes are padded by dLat
  auto withinLow = Point<double>(0.0, withinBox.getLowerLeft().getY() * 1.0);
  auto aLow = Point<double>(0.0, boxA.getLowerLeft().getY() * 1.0);
  aLow.setY(aLow.getY() - dLat);
  auto bLow = Point<double>(0.0, boxB.getLowerLeft().getY() * 1.0);
  bLow.setY(bLow.getY() - dLat);

  // the upper extrema of the bounding boxes are padded by dLat
  auto withinUp = Point<double>(0.0, withinBox.getUpperRight().getY() * 1.0);
  auto aUp = Point<double>(0.0, boxA.getUpperRight().getY() * 1.0);
  aUp.setY(aUp.getY() + dLat);
  auto bUp = Point<double>(0.0, boxB.getUpperRight().getY() * 1.0);
  bUp.setY(bUp.getY() + dLat);

  double yRangeMin = std::min(
      90.0 - util::geo::EPSILON,
      std::max(withinLow.getY() * 1.0, std::max(aLow.getY(), bLow.getY())));
  double yRangeMax = std::max(
      -90.0 - util::geo::EPSILON,
      std::min(withinUp.getY() * 1.0, std::min(aUp.getY(), bUp.getY())));

  double a = cos(yRangeMin * util::geo::RAD);
  double b = cos(yRangeMax * util::geo::RAD);

  // if we crossed the pole, we encountered a scale factor of 1!
  if (withinLow.getY() < 0 && withinUp.getY() > 0) {
    return {std::min(a, b), std::max(1.0, std::max(a, b))};
  }

  return {std::min(a, b), std::max(a, b)};
}

// _____________________________________________________________________________
template <template <typename> class GeomA, template <typename> class GeomB,
          typename T>
double withinMeterDist(const std::vector<GeomA<T>>& a, const GeomB<T>& b,
                       double maxD) {
  if (a.size() > EST_MULTI_CHECKS_THRESHOLD_XSORTED)
    return withinMeterDist(XSortedCollection<T>(a), XSortedCollection<T>(b),
                           maxD);

  auto scale = getMinMaxLocalScaleFactors(util::geo::getBoundingBox(a),
                                          util::geo::getBoundingBox(b), maxD);

  double maxEuclideanDist = maxD / scale.first;
  return withinDist(a, b, maxD, &meterDistLocalSearchPadding<T>,
                    maxEuclideanDist,
                    [](const Point<T> a, const Point<T> b, double) -> double {
                      return haversine(a, b);
                    });
}

// _____________________________________________________________________________
template <template <typename> class GeomA, template <typename> class GeomB,
          typename T>
double withinMeterDist(const std::vector<GeomA<T>>& a,
                       const std::vector<GeomB<T>>& b, double maxD) {
  if (a.size() * b.size() > EST_MULTI_CHECKS_THRESHOLD_XSORTED)
    return withinMeterDist(XSortedCollection<T>(a), XSortedCollection<T>(b),
                           maxD);

  auto scale = getMinMaxLocalScaleFactors(util::geo::getBoundingBox(a),
                                          util::geo::getBoundingBox(b), maxD);

  double maxEuclideanDist = maxD / scale.first;
  return withinDist(a, b, maxD, &meterDistLocalSearchPadding<T>,
                    maxEuclideanDist,
                    [](const Point<T> a, const Point<T> b, double) -> double {
                      return haversine(a, b);
                    });
}

}  // namespace geo
}  // namespace util
