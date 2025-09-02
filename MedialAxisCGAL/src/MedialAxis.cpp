#include "MedialAxis.hpp"

static std::optional<CgalSegment> clipCgalSegmentToPolygon(const CgalSegment& seg, const Polygon_2& poly) {
    // Collect intersection points with polygon edges
    std::vector<Point> intersections;

    // Iterate over polygon edges
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        const CgalSegment& edge = *eit;

        CGAL::Object result = CGAL::intersection(seg, edge);

        if (const Point* ipoint = CGAL::object_cast<Point>(&result)) {
            intersections.push_back(*ipoint);
        }
        else if (const CgalSegment* iseg = CGAL::object_cast<CgalSegment>(&result)) {
            // The CgalSegment lies on the polygon edge: keep endpoints
            intersections.push_back(iseg->source());
            intersections.push_back(iseg->target());
        }
    }

    // Also include the endpoints that are inside the polygon
    if (poly.bounded_side(seg.source()) == CGAL::ON_BOUNDED_SIDE ||
        poly.bounded_side(seg.source()) == CGAL::ON_BOUNDARY) {
        intersections.push_back(seg.source());
    }
    if (poly.bounded_side(seg.target()) == CGAL::ON_BOUNDED_SIDE ||
        poly.bounded_side(seg.target()) == CGAL::ON_BOUNDARY) {
        intersections.push_back(seg.target());
    }

    // If we found less than 2 points, nothing to keep
    if (intersections.size() < 2) {
        return std::nullopt;
    }

    // Pick the two extreme points along the original CgalSegment
    auto cmp = [&seg](const Point& a, const Point& b) {
        return CGAL::squared_distance(seg.source(), a) < CGAL::squared_distance(seg.source(), b);
        };
    auto [pmin, pmax] = std::minmax_element(intersections.begin(), intersections.end(), cmp);

    return CgalSegment(*pmin, *pmax);
}

static void removePointIfExists(std::vector<Point>& vertices, const Point& toRemove)
{
    auto it = std::find(vertices.begin(), vertices.end(), toRemove);
    if (it != vertices.end())
    {
        vertices.erase(it);
    }
}

/// Computes the intersection point of two CgalLines
/// @param l1 first CgalLine
/// @param l2 second CgalLine
/// @return an intersection point, if such exists; otherwise `std::nullopt`
static std::optional<Point> getIntersectionPoint(const CgalLine& l1, const CgalLine& l2) {
    const CGAL::Object obj = intersection(l1, l2);
    if (const Point* p = CGAL::object_cast<Point>(&obj)) {
        return *p;
    }
    return std::nullopt;
}

/// Computes the angle bisector of edges (prev, curr) and (curr, next)
/// @param prev first vertex
/// @param curr middle vertex
/// @param next end vertex
/// @return a CgalLine object representing the angle bisector of the angle between said edges
static CgalLine get_angle_bisector(const Point& prev, const Point& curr, const Point& next) {
    // Construct vectors from current point to neighbors
    K::Vector_2 u1 = prev - curr;
    K::Vector_2 u2 = next - curr;

    if (u1.squared_length() == 0 || u2.squared_length() == 0) {
        throw std::runtime_error("Degenerate edge when computing angle bisector");
    }

    // Normalize
    u1 = u1 / std::sqrt(u1.squared_length());
    u2 = u2 / std::sqrt(u2.squared_length());

    // Bisector direction
    K::Vector_2 bisector = u1 + u2;
    if (bisector.squared_length() == 0) {
        bisector = u2;
    }

    return { curr, curr + bisector };
}

MedialAxis::MedialAxis(const Polygon_2& pgn)
{
    if (!pgn.is_convex() || !pgn.is_simple()) {
        throw std::runtime_error("This algorithm only supports simple convex polygons.");
    }

    std::vector vertices(pgn.vertices_begin(), pgn.vertices_end());
    constexpr std::size_t TRIANGLE_VERTICES = 3;

    if (vertices.size() < TRIANGLE_VERTICES)
    {
        throw std::runtime_error("Not enough vertices in chosen polygon.");
    }

    while (vertices.size() > TRIANGLE_VERTICES) {
        CgalLinePair meeting_edges;
        Point center;

        auto next_meeting_bisectors = findNextCgalSegmentPair(vertices, meeting_edges, center);
        addMedialAxisCgalSegments(next_meeting_bisectors);
        updateVertices(vertices, next_meeting_bisectors, meeting_edges);
    }

    // Final triangle case
    triangleMedialAxis(vertices);
    clipToPolygon(pgn);
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

        CgalLine bisector1 = get_angle_bisector(vertices[i0], vertices[i1], vertices[i2]);
        CgalLine bisector2 = get_angle_bisector(vertices[i1], vertices[i2], vertices[i3]);

        auto inter = getIntersectionPoint(bisector1, bisector2);
        if (!inter)
        {
            continue;
        }

        const double r = std::sqrt(CGAL::squared_distance(*inter, vertices[i1]));
        if (r < min_radius) {
            min_radius = r;
            meeting_edges = { CgalLine(vertices[i0], vertices[i1]), CgalLine(vertices[i2], vertices[i3]) };
            center = *inter;
            earliest_meeting_pair = { {vertices[i1], center}, {vertices[i2], center} };
        }
    }

    return earliest_meeting_pair;
}

void MedialAxis::addMedialAxisCgalSegments(const CgalSegmentPair& earliest_meeting_pair)
{
    m_medialAxisCgalSegments.emplace_back(earliest_meeting_pair.first);
    m_medialAxisCgalSegments.emplace_back(earliest_meeting_pair.second);
}

void MedialAxis::updateVertices(std::vector<Point>& vertices, const CgalSegmentPair& earliest_meeting_pair, const CgalLinePair& meeting_edges) const
{
    auto point = getIntersectionPoint(meeting_edges.first, meeting_edges.second);
    if (point.has_value())
    {
        vertices.push_back(point.value());
    }
    else
    {
        throw std::runtime_error("No valid edge intersection found; The supplied polygon probably fits a degenerate case.");
    }
    removePointIfExists(vertices, earliest_meeting_pair.first.point(0));
    removePointIfExists(vertices, earliest_meeting_pair.second.point(0));
}

void MedialAxis::triangleMedialAxis(const std::vector<Point>& vertices)
{
    const Point A = vertices[0], B = vertices[1], C = vertices[2];
    Point center = CGAL::centroid(A, B, C);
    m_medialAxisCgalSegments.emplace_back(A, center);
    m_medialAxisCgalSegments.emplace_back(B, center);
    m_medialAxisCgalSegments.emplace_back(C, center);
}

void MedialAxis::clipToPolygon(const Polygon_2& poly) {
    std::list<CgalSegment> clipped;
    for (const auto& seg : m_medialAxisCgalSegments) {
        auto clippedSeg = clipCgalSegmentToPolygon(seg, poly);
        if (clippedSeg) {
            clipped.push_back(*clippedSeg);
        }
    }
    m_medialAxisCgalSegments = clipped;
}

std::list<CgalSegment> MedialAxis::get() const
{
    return m_medialAxisCgalSegments;
}
