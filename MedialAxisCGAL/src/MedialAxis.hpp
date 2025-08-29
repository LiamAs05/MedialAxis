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
using Line = K::Line_2;
using Segment = K::Segment_2;
using Polygon_2 = CGAL::Polygon_2<K>;
using LinePair = std::pair<Line, Line>;
using SegmentPair = std::pair<Segment, Segment>;

class MedialAxis
{
public:
    explicit MedialAxis(const Polygon_2& pgn);
    std::list<Segment> get() const;

private:
    /// Finds the pair of adjacent bisectors that meet first.
    /// @param vertices A convex polygon in CCW order
    /// @param meeting_edges (out) The two polygon edges adjacent to the collapse
    /// @param center (out) Intersection point of bisectors
    /// @return Two medial-axis segments connecting center to the two collapsing vertices
    SegmentPair findNextSegmentPair(const std::vector<Point>& vertices, LinePair& meeting_edges, Point& center) const;

    void addMedialAxisSegments(const SegmentPair& earliest_meeting_pair);

    void updateVertices(std::vector<Point>& vertices, const SegmentPair& earliest_meeting_pair, const LinePair& meeting_edges) const;

    void triangleMedialAxis(const std::vector<Point>& vertices);

    void clipToPolygon(const Polygon_2& poly);

    std::list<Segment> m_medialAxisSegments;
};