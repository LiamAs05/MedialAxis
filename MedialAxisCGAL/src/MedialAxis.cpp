#include "MedialAxis.hpp"

#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iterator>

constexpr std::size_t TRIANGLE_VERTICES = 3;

MedialAxis::MedialAxis(const Polygon_2& pgn) : m_originalPolygon(pgn), m_clipper(pgn) 
{
    if (!pgn.is_convex() || !pgn.is_simple()) {
        throw std::runtime_error("This algorithm only supports simple convex polygons.");
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
        addMedialAxisCgalSegments(next_meeting_bisectors);
        updateVertices(vertices, next_meeting_bisectors, meeting_edges);
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

        if (std::abs(d1 - d2) > CLIP_EPS || std::abs(d2 - d3) > CLIP_EPS)
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
    auto seg1 = clipCgalSegmentToPolygon(earliest_meeting_pair.first, m_clipper);
    auto seg2 = clipCgalSegmentToPolygon(earliest_meeting_pair.second, m_clipper);

    if (seg1.has_value())
        m_medialAxisSegments.emplace_back(seg1.value());
    if (seg2.has_value())
        m_medialAxisSegments.emplace_back(seg2.value());
}

void MedialAxis::updateClipperPolygon(const Point& src1, const Point& src2, const Point& target)
{
    // Extract vertices of current clipping polygon
    std::vector<Point> pts(m_clipper.vertices_begin(), m_clipper.vertices_end());
    const std::size_t n = pts.size();

    if (n < TRIANGLE_VERTICES) {
        return;
    }

    // find indices of the two source vertices with tolerance
    int idx1 = -1, idx2 = -1;
    for (std::size_t i = 0; i < n; ++i) {
        if (idx1 < 0 && pointsEqualEps(pts[i], src1)) idx1 = static_cast<int>(i);
        if (idx2 < 0 && pointsEqualEps(pts[i], src2)) idx2 = static_cast<int>(i);
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

    CGAL::write_polygon_WKT(std::cout, m_clipper);
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
    updateClipperPolygon(earliest_meeting_pair.first.source(),
                           earliest_meeting_pair.second.source(),
                           earliest_meeting_pair.first.target());
}


void MedialAxis::triangleMedialAxis(const std::vector<Point>& vertices)
{
    const Point A = vertices[0], B = vertices[1], C = vertices[2];
    Point center = computeTriangleIncenter(A, B, C); // use incenter instead of centroid

    for (const auto& point : vertices)
    {
        auto clippedSegment = clipCgalSegmentToPolygon({point, center}, m_clipper);
        if (clippedSegment.has_value())
            m_medialAxisSegments.push_back(clippedSegment.value());
    }
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
