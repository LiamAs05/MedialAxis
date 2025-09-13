#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/WKT.h>
#include <list>
#include <iostream>
#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_2;
using CgalLine = K::Line_2;
using CgalSegment = K::Segment_2;
using Polygon_2 = CGAL::Polygon_2<K>;
using CgalLinePair = std::pair<CgalLine, CgalLine>;
using CgalSegmentPair = std::pair<CgalSegment, CgalSegment>;

constexpr double CLIP_EPS = 1e-9;

/// Compare two points with a small coordinate tolerance
/// @param a, b points to compare
/// @param eps tolerance
/// @return true if equal
bool pointsEqualEps(const Point& a, const Point& b, double eps = CLIP_EPS);

/// Clips a segment to a given polygon
/// @param seg segment to clip
/// @param poly polygon to clip to
/// @return Segment representing the newly-clipped segment, or `std::nullopt` if the polygon and the segment don't intersect.
std::optional<CgalSegment> clipCgalSegmentToPolygon(const CgalSegment& seg, const Polygon_2& poly);

/// Removes a point if it exists
/// @param vertices container to lookup point in
/// @param toRemove candidate to remove
void removePointIfExists(std::vector<Point>& vertices, const Point& toRemove);

/// Computes the intersection point of two CgalLines
/// @param l1 first CgalLine
/// @param l2 second CgalLine
/// @return The intersection point, if such exists; otherwise `std::nullopt`
std::optional<Point> getIntersectionPoint(const CgalLine& l1, const CgalLine& l2);

/// Computes the angle bisector of edges (prev, curr) and (curr, next)
/// @param prev first vertex
/// @param curr middle vertex
/// @param next end vertex
/// @return Line object representing the angle bisector of the angle between said edges
CgalLine computeAngleBisector(const Point& prev, const Point& curr, const Point& next);

/// Computes the incenter of a triangle https://en.wikipedia.org/wiki/Incenter
/// @param A, B, C endpoints of the triangle
/// @return point representing the intersection of all angle bisectors of the triangle
Point computeTriangleIncenter(const Point& A, const Point& B, const Point& C);