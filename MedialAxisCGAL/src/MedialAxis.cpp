#include "MedialAxis.hpp"

#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iterator>

constexpr std::size_t TRIANGLE_VERTICES = 3;

MedialAxis::MedialAxis(const Polygon_2& pgn) : m_originalPolygon(pgn)
{
    if (!pgn.is_convex() || !pgn.is_simple()) {
        throw std::runtime_error("This algorithm only supports simple convex polygons.");
    }

    for (const auto vertex : pgn.vertices())
    {
        m_polygonCorrespondents[vertex] = vertex;
    }

    std::vector vertices(pgn.vertices_begin(), pgn.vertices_end());

    if (vertices.size() < TRIANGLE_VERTICES)
    {
        throw std::runtime_error("Not enough vertices in chosen polygon.");
    }

    while (vertices.size() > TRIANGLE_VERTICES) {
        CgalLinePair meeting_edges;
        Point center;

        auto next_meeting_bisectors = findNextCgalSegmentPair(vertices, meeting_edges, center);
        updateVertices(vertices, next_meeting_bisectors, meeting_edges);
        addMedialAxisCgalSegments(next_meeting_bisectors);
    }

    triangleMedialAxis(vertices);
}

CgalSegmentPair MedialAxis::findNextCgalSegmentPair(const std::vector<Point>& vertices, CgalLinePair& meeting_edges, Point& center) const
{
    std::pair<CgalSegment, CgalSegment> earliest_meeting_pair;
    double min_radius = std::numeric_limits<double>::max();
    std::size_t n = vertices.size();

    for (std::uint64_t i = 0; i < n; ++i) {
        const std::uint64_t i0 = (i - 1 + n) % n;
        const std::uint64_t i1 = i;
        const std::uint64_t i2 = (i + 1) % n;
        const std::uint64_t i3 = (i + 2) % n;

        CgalLine bisector1 = computeAngleBisector(vertices[i0], vertices[i1], vertices[i2]);
        CgalLine bisector2 = computeAngleBisector(vertices[i1], vertices[i2], vertices[i3]);

        auto inter = getIntersectionPoint(bisector1, bisector2);
        if (!inter.has_value())
        {
            continue;
        }

        CgalLine L1(vertices[i0], vertices[i1]);
        CgalLine L2(vertices[i1], vertices[i2]);
        CgalLine L3(vertices[i2], vertices[i3]);

        const double d1 = std::sqrt(CGAL::to_double(CGAL::squared_distance(inter.value(), L1)));
        const double d2 = std::sqrt(CGAL::to_double(CGAL::squared_distance(inter.value(), L2)));
        const double d3 = std::sqrt(CGAL::to_double(CGAL::squared_distance(inter.value(), L3)));

        if (std::abs(d1 - d2) > FP_TOLERANCE || std::abs(d2 - d3) > FP_TOLERANCE)
        {
            throw std::logic_error("Computation error: Incorrect results in computing radius of circle");
        }

        // Radius of circle centered at intersection point and tangent to L1, L2, L3
        const double r = (d1 + d2 + d3) / 3.0;
        if (r < min_radius) {
            min_radius = r;
            meeting_edges = { L1, L3 };
            center = inter.value();
            earliest_meeting_pair = { {vertices[i1], center}, {vertices[i2], center} };
        }
    }

    return earliest_meeting_pair;
}

void MedialAxis::addMedialAxisCgalSegments(const CgalSegmentPair& earliest_meeting_pair)
{
    if (m_polygonCorrespondents.find(earliest_meeting_pair.first.source()) == m_polygonCorrespondents.end())
    {
        throw std::runtime_error("Sad liam");
    }
    CgalSegment seg1 = {m_polygonCorrespondents[earliest_meeting_pair.first.source()], earliest_meeting_pair.first.target()};
    CgalSegment seg2 = {m_polygonCorrespondents[earliest_meeting_pair.second.source()], earliest_meeting_pair.second.target()};

    m_medialAxisSegments.push_back(seg1);
    m_medialAxisSegments.push_back(seg2);
}

void MedialAxis::updateVertices(std::vector<Point>& vertices,
                                const CgalSegmentPair& earliest_meeting_pair,
                                const CgalLinePair& meeting_edges)
{
    auto inter = getIntersectionPoint(meeting_edges.first, meeting_edges.second);
    if (!inter) {
        throw std::runtime_error("No valid edge intersection found; degenerate case.");
    }
    Point newVertex = *inter;

    // Find the collapsing vertices
    auto it1 = std::find(vertices.begin(), vertices.end(),
                         earliest_meeting_pair.first.source());
    auto it2 = std::find(vertices.begin(), vertices.end(),
                         earliest_meeting_pair.second.source());

    if (it1 == vertices.end() || it2 == vertices.end()) {
        throw std::runtime_error("Could not find collapsing vertices in polygon.");
    }

    // Ensure order: it1 comes before it2
    if (std::distance(vertices.begin(), it2) < std::distance(vertices.begin(), it1))
        std::swap(it1, it2);

    // Replace the first collapsed vertex with the intersection
    *it1 = newVertex;

    // Erase the second collapsed vertex
    vertices.erase(it2);

    m_polygonCorrespondents[newVertex] = earliest_meeting_pair.first.target();
}


void MedialAxis::triangleMedialAxis(const std::vector<Point>& vertices)
{
    const Point A = vertices[0], B = vertices[1], C = vertices[2];
    Point center = computeTriangleIncenter(A, B, C); // use incenter instead of centroid

    for (const auto& point : vertices)
    {
        CgalSegment clipped = {m_polygonCorrespondents[point], center};
        m_medialAxisSegments.push_back(clipped);
    }
}

std::list<CgalSegment> MedialAxis::get() const
{
    return m_medialAxisSegments;
}
