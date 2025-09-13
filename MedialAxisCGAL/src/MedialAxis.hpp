#pragma once

#include "GeometryUtils.hpp"

class MedialAxis
{
public:
    explicit MedialAxis(const Polygon_2& pgn);
    std::list<CgalSegment> get() const;

private:
    /// Finds the pair of adjacent bisectors that meet first.
    /// @param vertices A convex polygon in CCW order
    /// @param meeting_edges (out) The two polygon edges adjacent to the collapse
    /// @param center (out) Intersection point of bisectors
    /// @return Two medial-axis CgalSegments connecting center to the two collapsing vertices
    CgalSegmentPair findNextCgalSegmentPair(const std::vector<Point>& vertices, CgalLinePair& meeting_edges, Point& center) const;

    void addMedialAxisCgalSegments(const CgalSegmentPair& earliest_meeting_pair);

    void updateClipperPolygon(const Point& src1, const Point& src2, const Point& target);

    void updateVertices(std::vector<Point>& vertices, const CgalSegmentPair& earliest_meeting_pair, const CgalLinePair& meeting_edges);

    void triangleMedialAxis(const std::vector<Point>& vertices);

    void clipToPolygon(const Polygon_2& poly);

    std::list<CgalSegment> m_medialAxisSegments;
    Polygon_2 m_clipper;
    Polygon_2 m_originalPolygon;
};
