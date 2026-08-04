// Copyright 2025, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "./Geo.h"

// _____________________________________________________________________________
uint8_t util::geo::boolArrToInt8(const std::array<bool, 8> arr) {
  uint8_t ret = 0;
  for (size_t i = 0; i < 8; i++) ret |= (uint8_t)arr[i] << i;
  return ret;
}

// _____________________________________________________________________________
bool util::geo::doubleEq(double a, double b) { return fabs(a - b) < EPSILON; }

// _____________________________________________________________________________
double util::geo::innerProd(double x1, double y1, double x2, double y2,
                            double x3, double y3) {
  double dx21 = x2 - x1;
  double dx31 = x3 - x1;
  double dy21 = y2 - y1;
  double dy31 = y3 - y1;
  double m12 = sqrt(dx21 * dx21 + dy21 * dy21);
  double m13 = sqrt(dx31 * dx31 + dy31 * dy31);
  double theta = acos(std::min((dx21 * dx31 + dy21 * dy31) / (m12 * m13), 1.0));

  return theta * IRAD;
}

// _____________________________________________________________________________
double util::geo::crossProd(double x1, double y1, double x2, double y2) {
  return x1 * y2 - x2 * y1;
}

// _____________________________________________________________________________
util::geo::WKTType util::geo::getWKTType(const char* c, const char** endr) {
  bool measurement = false;
  while(true) {

    // TODO: how do we handle another <, is there need for further checks?

    // Check for possible IRI and skip it.
    if (*c == '<') {
      if (measurement) break; // Measurement M cannot be in front of IRI.
      while (*c && *c != '>') c++;
      if (*c == '>') c++; // Also skip '>'.
      continue;
    }
    if ((*c == ' ' || *c == '\n' || *c == '\t' || *c == '\r') ||
         ((*c) == '"') || ((*c) == '\'')) {
      c++; // skip possible whitespace
      continue;
    }
    if ((tolower(*c) == 'm' && ((*(c + 1) == ' ' || *(c + 1) == '\n' || *(c + 1) == '\t' || 
          *(c + 1) == '\r') || tolower(*(c + 1)) != 'u'))) {
      c++; // skip possible measurement M
      measurement = true;
      continue;
    }
    break;
  }
  if (strncicmp("POINT", c, 5) == 0) {
    if (endr) (*endr) = c + 5;
    return POINT;
  }
  if (strncicmp("LINESTRING", c, 10) == 0) {
    if (endr) (*endr) = c + 10;
    return LINESTRING;
  }
  if (strncicmp("POLYGON", c, 7) == 0) {
    if (endr) (*endr) = c + 7;
    return POLYGON;
  }
  if (strncicmp("MULTIPOINT", c, 10) == 0) {
    if (endr) (*endr) = c + 10;
    return MULTIPOINT;
  }
  if (strncicmp("MULTILINESTRING", c, 15) == 0) {
    if (endr) (*endr) = c + 15;
    return MULTILINESTRING;
  }
  if (strncicmp("MULTIPOLYGON", c, 12) == 0) {
    if (endr) (*endr) = c + 12;
    return MULTIPOLYGON;
  }
  if (strncicmp("GEOMETRYCOLLECTION", c, 18) == 0) {
    if (endr) (*endr) = c + 18;
    return COLLECTION;
  }

  if (endr) (*endr) = 0;
  return NONE;
}


// _____________________________________________________________________________
util::geo::CRSType util::geo::getCRSType(const char* c, const char** endr) {
  while (((*c == ' ' || *c == '\n' || *c == '\t' || *c == '\r') ||
         ((*c) == '"') || ((*c) == '\'')) && ((*c) != '\0'))
    c++; // Skip possible whitespace.

  // If 'endr == nullptr' this function should still update it, for this
  // 'endr' is replaced accordingly.
  const char* replacement = nullptr;
  endr = (endr != nullptr) ? endr : &replacement;

  if (*c != '<') {
    if (endr) (*endr) = c;
    return CRS84;  // Default.
  }

  if (strncicmp("<http://www.opengis.net/def/crs/OGC/1.3/CRS84>", c, 46) == 0) {
    if (endr) (*endr) = c + 46;
    return CRS84;
  }
  if (strncicmp("<http://www.opengis.net/def/crs/EPSG/0/4326>", c, 44) == 0) {
    if (endr) (*endr) = c + 44;
    return WGS84;
  }
  if (strncicmp("<http://www.opengis.net/def/crs/EPSG/0/3857>", c, 44) == 0) {
    if (endr) (*endr) = c + 44;
    return WEB_MERCATOR;
  }

  if (endr) (*endr) = c;
  return UNSUPPORTED;
}

// _____________________________________________________________________________
double util::geo::distToSegment(double lax, double lay, double lbx, double lby,
                                double px, double py) {
  double dx = lbx - lax;
  double dy = lby - lay;
  double d = dx * dx + dy * dy;
  if (d == 0) return dist(px, py, lax, lay);

  double dot = (px - lax) * dx + (py - lay) * dy;
  if (dot <= 0) return dist(px, py, lax, lay);
  if (dot >= d) return dist(px, py, lbx, lby);

  double t = dot / d;
  return dist(px, py, lax + t * dx, lay + t * dy);
}

// _____________________________________________________________________________
double util::geo::distToSegmentSquared(double lax, double lay, double lbx,
                                       double lby, double px, double py) {
  double dx = lbx - lax;
  double dy = lby - lay;
  double d = dx * dx + dy * dy;
  if (d == 0) return distSquared(px, py, lax, lay);

  double dot = (px - lax) * dx + (py - lay) * dy;
  if (dot <= 0) return distSquared(px, py, lax, lay);
  if (dot >= d) return distSquared(px, py, lbx, lby);

  double t = dot / d;
  return distSquared(px, py, lax + t * dx, lay + t * dy);
}

// _____________________________________________________________________________
util::geo::WKTType util::geo::getWKTType(const char* c) {
  return util::geo::getWKTType(c, 0);
}

// _____________________________________________________________________________
util::geo::CRSType util::geo::getCRSType(const char* c) { return util::geo::getCRSType(c, 0); }

// _____________________________________________________________________________
util::geo::WKTType util::geo::getWKTType(const std::string& str) {
  return util::geo::getWKTType(str.c_str(), 0);
}

// _____________________________________________________________________________
util::geo::CRSType util::geo::getCRSType(const std::string& str) {
  return util::geo::getCRSType(str.c_str(), 0);
}

// _____________________________________________________________________________
double util::geo::distSquared(double x1, double y1, double x2, double y2) {
  return (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
}

// _____________________________________________________________________________
double util::geo::dist(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

// _____________________________________________________________________________
double util::geo::angBetween(double p1x, double p1y) { return atan2(p1x, p1y); }
