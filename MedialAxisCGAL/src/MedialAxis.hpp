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
    /// @param polygonEdgesToCollapse (out) The two polygon edges adjacent to the collapse
    /// @param bisectorIntersection (out) Intersection point of bisectors
    /// @return Two medial-axis CgalSegments connecting inter to the two collapsing vertices
    CgalSegmentPair findNextCgalSegmentPair(const std::vector<Point>& vertices, CgalLinePair& polygonEdgesToCollapse, Point& bisectorIntersection) const;

    /// Update the medial axis tree with two new segments
    /// @param earliestIntersectingBisectors the pair of bisectors to intersect next
    void addMedialAxisCgalSegments(const CgalSegmentPair& earliestIntersectingBisectors);

    /// Collapse the reduandent edge and update the polygons vertices
    /// @param vertices The current n-gon, which will be changed to an (n-1)-gon
    /// @param earliestIntersecting The pair of bisectors to intersect next
    /// @param polygonEdgesToCollapse The neighbor edges to that which will be deleted, used to calculate the new vertex to add
    void updateVertices(std::vector<Point>& vertices, const CgalSegmentPair& earliestIntersectingBisectors, const CgalLinePair& polygonEdgesToCollapse);

    /// Compute the medial axis of a triangle - it's incenter
    /// @param vertices the vertices of the triangle
    void triangleMedialAxis(const std::vector<Point>& vertices);

    std::list<CgalSegment> m_medialAxisSegments;
    Polygon_2 m_originalPolygon;
    std::unordered_map<Point, Point> m_polygonCorrespondents;
};
