#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/WKT.h>
#include <iostream>
#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_2;
using Polygon_2 = CGAL::Polygon_2<K> ;
using Line = K::Line_2;
using Segment = K::Segment_2;
using LinePair = std::pair<Line, Line>;
using SegmentPair = std::pair<Segment, Segment>;

std::optional<Segment> clipSegmentToPolygon(const Segment& seg, const Polygon_2& poly) {
    // Collect intersection points with polygon edges
    std::vector<Point> intersections;

    // Iterate over polygon edges
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        const Segment& edge = *eit;

        CGAL::Object result = CGAL::intersection(seg, edge);

        if (const Point* ipoint = CGAL::object_cast<Point>(&result)) {
            intersections.push_back(*ipoint);
        }
        else if (const Segment* iseg = CGAL::object_cast<Segment>(&result)) {
            // The segment lies on the polygon edge: keep endpoints
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

    // Pick the two extreme points along the original segment
    auto cmp = [&seg](const Point& a, const Point& b) {
        return CGAL::squared_distance(seg.source(), a) < CGAL::squared_distance(seg.source(), b);
        };
    auto [pmin, pmax] = std::minmax_element(intersections.begin(), intersections.end(), cmp);

    return Segment(*pmin, *pmax);
}

void removePointIfExists(std::vector<Point>& vertices, const Point& toRemove)
{
    auto it = std::find(vertices.begin(), vertices.end(), toRemove);
    if (it != vertices.end())
    {
        vertices.erase(it);
    }
}

class MedialAxis
{
public:
    explicit MedialAxis(const Polygon_2& pgn);
    std::list<Segment> get() const;

private:
    SegmentPair findNextSegmentPair(const std::vector<Point>& vertices, LinePair& meeting_edges, Point& center) const;

	void addMedialAxisSegments(const SegmentPair& earliest_meeting_pair);

    void updateVertices(std::vector<Point>& vertices, const SegmentPair& earliest_meeting_pair, const LinePair& meeting_edges) const;

	void triangleMedialAxis(const std::vector<Point>& vertices);

    void clipToPolygon(const Polygon_2& poly);

	std::list<Segment> m_medialAxisSegments;
};

/// Computes the intersection point of two lines
/// @param l1 
/// @param l2 
/// @return an intersection point, if such exists; otherwise `std::nullopt`
static std::optional<Point> getIntersectionPoint(const Line& l1, const Line& l2) {
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
/// @return a Line object representing the angle bisector of the angle between said edges
Line get_angle_bisector(const Point& prev, const Point& curr, const Point& next) {
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

/// Finds the pair of adjacent bisectors that meet first.
/// @param vertices A convex polygon in CCW order
/// @param meeting_edges (out) The two polygon edges adjacent to the collapse
/// @param center (out) Intersection point of bisectors
/// @return Two medial-axis segments connecting center to the two collapsing vertices
SegmentPair MedialAxis::findNextSegmentPair(const std::vector<Point>& vertices, LinePair& meeting_edges, Point& center) const
{
    std::pair<Segment, Segment> earliest_meeting_pair;
    double min_radius = std::numeric_limits<double>::max();
    std::size_t n = vertices.size();

	for (std::uint64_t i = 0; i < n; ++i) {
		const std::uint64_t i0 = (i - 1 + n) % n;
		const std::uint64_t i1 = i;
		const std::uint64_t i2 = (i + 1) % n;
		const std::uint64_t i3 = (i + 2) % n;

		Line bisector1 = get_angle_bisector(vertices[i0], vertices[i1], vertices[i2]);
		Line bisector2 = get_angle_bisector(vertices[i1], vertices[i2], vertices[i3]);

		auto inter = getIntersectionPoint(bisector1, bisector2);
		if (!inter)
		{
			continue;
		}

		const double r = std::sqrt(CGAL::squared_distance(*inter, vertices[i1]));
		if (r < min_radius) {
			min_radius = r;
			meeting_edges = { Line(vertices[i0], vertices[i1]), Line(vertices[i2], vertices[i3]) };
			center = *inter;
            earliest_meeting_pair = { {vertices[i1], center}, {vertices[i2], center} };
		}
	}

    return earliest_meeting_pair;
}

void MedialAxis::addMedialAxisSegments(const SegmentPair& earliest_meeting_pair)
{
    m_medialAxisSegments.emplace_back(earliest_meeting_pair.first);
    m_medialAxisSegments.emplace_back(earliest_meeting_pair.second);
}

void MedialAxis::updateVertices(std::vector<Point>& vertices, const SegmentPair& earliest_meeting_pair, const LinePair& meeting_edges) const
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
	m_medialAxisSegments.emplace_back(A, center);
	m_medialAxisSegments.emplace_back(B, center);
	m_medialAxisSegments.emplace_back(C, center);
}

void MedialAxis::clipToPolygon(const Polygon_2& poly) {
    std::list<Segment> clipped;
    for (const auto& seg : m_medialAxisSegments) {
        auto clippedSeg = clipSegmentToPolygon(seg, poly);
        if (clippedSeg) {
            clipped.push_back(*clippedSeg);
        }
    }
    m_medialAxisSegments = clipped;
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
        LinePair meeting_edges;
        Point center;

        auto next_meeting_bisectors = findNextSegmentPair(vertices, meeting_edges, center);
        addMedialAxisSegments(next_meeting_bisectors);
        updateVertices(vertices, next_meeting_bisectors, meeting_edges);
    }

    // Final triangle case
    triangleMedialAxis(vertices);
    clipToPolygon(pgn);
}

std::list<Segment> MedialAxis::get() const
{
    return m_medialAxisSegments;
}

Polygon_2 polygon_from_wkt()
{
    const std::string coordinates_file_path = R"(C:\Users\liamd\Desktop\ipe-7.2.29\bin\coordinates.txt)";
    std::ifstream wkt_file(coordinates_file_path);

    if (!wkt_file) {
        throw std::runtime_error("Could not open WKT file: " + coordinates_file_path);
    }

    Polygon_2 poly;
    if (!read_polygon_WKT(wkt_file, poly)) {
        throw std::runtime_error("Failed to parse WKT polygon.");
    }

    return poly;
}

int main()
{
    try {
        const Polygon_2 pgn = polygon_from_wkt();

        MedialAxis m(pgn);
        const auto segments = m.get();
        std::cout << "Medial Axis Edges: " << segments.size() << std::endl;

        for (const auto& seg : segments) {
            std::cout << seg.point(0) << " -> "
                << seg.point(1) << std::endl;
        }
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}


//namespace medial_axis_ipelet
//{
//    const std::string labels[] = { "Medial Axis", "Help" };
//
//    const std::string hmsg[] = {
//          "Compute the Medial Axis transformation of a simple convex polygon."
//    };
//
//    class MedialAxisIpelet : public CGAL::Ipelet_base<K, 2>
//    {
//    public:
//        MedialAxisIpelet() : Ipelet_base("Medial Axis", labels, hmsg) {}
//        void protected_run(int) override;
//    };
//
//    void MedialAxisIpelet::protected_run(int fn)
//    {
//        switch (fn) {
//        case 1:
//            show_help();//print an help message
//            return;
//        default:
//            std::list<Point_2> pt_lst;
//
//            // Recovering points using output iterator of typ Dispatch_or_drop_output_iterator
//            read_active_objects(
//                CGAL::dispatch_or_drop_output<Point_2>(std::back_inserter(pt_lst))
//            );
//
//            if (pt_lst.empty()) {
//                print_error_message("No mark selected");
//                return;
//            }
//
//            const Polygon_2 pgn = Polygon_2(pt_lst.begin(), pt_lst.end());
//            MedialAxis m(pgn);
//
//            for (const auto& seg : m.get())
//            {
//                draw_in_ipe(seg);
//            }
//        };
//    }
//
//}
//
//CGAL_IPELET(medial_axis_ipelet::MedialAxisIpelet)
