#include "MedialAxis.hpp"

#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iterator>

// tolerance for point comparisons (tweak if needed)
static constexpr double CLIP_EPS = 1e-9;

// Compare two points with a small coordinate tolerance
static bool points_equal_eps(const Point& a, const Point& b, double eps = CLIP_EPS) {
    return (std::abs(CGAL::to_double(a.x()) - CGAL::to_double(b.x())) <= eps) &&
           (std::abs(CGAL::to_double(a.y()) - CGAL::to_double(b.y())) <= eps);
}


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

MedialAxis::MedialAxis(const Polygon_2& pgn) : m_originalPolygon(pgn), m_clipper(pgn) 
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
    for (auto& vertex : vertices)
    {
        std::cout << vertex << std::endl;  
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
    auto seg1 = clipCgalSegmentToPolygon(earliest_meeting_pair.first, m_clipper);
    auto seg2 = clipCgalSegmentToPolygon(earliest_meeting_pair.second, m_clipper);

    if (seg1.has_value())
        m_medialAxisSegments.emplace_back(seg1.value());
    if (seg2.has_value())
        m_medialAxisSegments.emplace_back(seg2.value());
}

void MedialAxis::update_clipper_incremental(const Point& src1, const Point& src2, const Point& target)
{
    // Extract vertices of current clipping polygon
    std::vector<Point> pts(m_clipper.vertices_begin(), m_clipper.vertices_end());
    const std::size_t n = pts.size();
    if (n < 3) {
        return;
    }

    // find indices of the two source vertices with tolerance
    int idx1 = -1, idx2 = -1;
    for (std::size_t i = 0; i < n; ++i) {
        if (idx1 < 0 && points_equal_eps(pts[i], src1)) idx1 = static_cast<int>(i);
        if (idx2 < 0 && points_equal_eps(pts[i], src2)) idx2 = static_cast<int>(i);
        if (idx1 >= 0 && idx2 >= 0) break;
    }

    if (idx1 < 0 || idx2 < 0) {
        // Couldn't find one or both sources in m_clipper.
        // This can happen if numeric drift removed them earlier.
        // Fallback: try to find approximate matches by nearest neighbor.
        auto find_nearest = [&](const Point &s)->int {
            double bestd = std::numeric_limits<double>::max();
            int besti = -1;
            for (std::size_t i = 0; i < n; ++i) {
                double d = CGAL::to_double(CGAL::squared_distance(pts[i], s));
                if (d < bestd) { bestd = d; besti = static_cast<int>(i); }
            }
            // only accept if reasonably close
            if (bestd <= (CLIP_EPS*CLIP_EPS * 100.0)) return besti;
            return -1;
        };
        if (idx1 < 0) idx1 = find_nearest(src1);
        if (idx2 < 0) idx2 = find_nearest(src2);
    }

    if (idx1 < 0 || idx2 < 0) {
        // Give up and repair m_clipper by convex hull of (current pts + target).
        std::vector<Point> fallback_pts = pts;
        fallback_pts.push_back(target);
        std::vector<Point> hull;
        CGAL::convex_hull_2(fallback_pts.begin(), fallback_pts.end(), std::back_inserter(hull));
        if (hull.size() >= 3) m_clipper = Polygon_2(hull.begin(), hull.end());
        else m_clipper = Polygon_2();
        return;
    }

    // rotate so idx1 < idx2 in linear order (handle wrap-around)
    if (idx2 < idx1) {
        std::rotate(pts.begin(), pts.begin() + idx1, pts.end());
        // recompute indices after rotation
        idx2 = (idx2 + static_cast<int>(n) - idx1);
        idx1 = 0;
    }

    // Now idx1 < idx2 (in the rotated pts array). Replace pts[idx1] = target, erase pts[idx2].
    pts[idx1] = target;
    pts.erase(pts.begin() + idx2);

    // Build a polygon and check validity (convex/simple). If it fails, fall back to hull.
    Polygon_2 candidate(pts.begin(), pts.end());
    if (candidate.is_simple() && candidate.is_convex() && candidate.size() >= 3) {
        m_clipper = std::move(candidate);
    } else {
        // Repair: compute convex hull from pts
        std::vector<Point> hull;
        CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull));
        if (hull.size() >= 3) m_clipper = Polygon_2(hull.begin(), hull.end());
        else m_clipper = Polygon_2();
    }
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

    // after you compute newVertex and update `vertices` (your existing code)
    update_clipper_incremental(earliest_meeting_pair.first.source(),
                           earliest_meeting_pair.second.source(),
                           earliest_meeting_pair.first.target());
}


static Point compute_incenter(const Point& A, const Point& B, const Point& C)
{
    double a = std::sqrt(CGAL::squared_distance(B, C)); // length opposite A
    double b = std::sqrt(CGAL::squared_distance(C, A)); // length opposite B
    double c = std::sqrt(CGAL::squared_distance(A, B)); // length opposite C

    double px = a * A.x() + b * B.x() + c * C.x();
    double py = a * A.y() + b * B.y() + c * C.y();
    double denom = a + b + c;

    return Point(px / denom, py / denom);
}

void MedialAxis::triangleMedialAxis(const std::vector<Point>& vertices)
{
    const Point A = vertices[0], B = vertices[1], C = vertices[2];
    Point center = compute_incenter(A, B, C); // use incenter instead of centroid

    CGAL::write_polygon_WKT(std::cout, m_clipper);
    auto a = clipCgalSegmentToPolygon({A, center}, m_clipper);
    auto b = clipCgalSegmentToPolygon({B, center}, m_clipper);
    auto c = clipCgalSegmentToPolygon({C, center}, m_clipper);
    
    if (a.has_value())
        m_medialAxisSegments.push_back(a.value());
    if (b.has_value())
        m_medialAxisSegments.push_back(b.value());
    if (c.has_value())
        m_medialAxisSegments.push_back(c.value());
}

void MedialAxis::clipToPolygon(const Polygon_2& poly) {
    std::list<CgalSegment> clipped;
    for (const auto& seg : m_medialAxisSegments) {
        auto clippedSeg = clipCgalSegmentToPolygon(seg, poly);
        if (clippedSeg.has_value()) {
            clipped.push_back(clippedSeg.value());
        }
    }
    m_medialAxisSegments = clipped;
}

std::list<CgalSegment> MedialAxis::get() const
{
    return m_medialAxisSegments;
}
