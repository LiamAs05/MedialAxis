#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/IO/WKT.h>
#include <iostream>
#include <algorithm>
//#include <CGAL/CGAL_Ipelet_base.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
typedef K::Line_2 Line;
typedef K::Segment_2 Segment;

using std::cout;
using std::endl;

class MedialAxis
{
public:
    explicit MedialAxis(const Polygon_2& pgn);
    std::list<Segment> getMedialAxis() const;

private:
    static std::pair<Segment, Segment> findNextSegmentPair(std::vector<CGAL::Point_2<CGAL::Epick>> vertices,
        std::pair<Line, Line>& meeting_edges, Point& center);

	void addMedialAxisSegments(Polygon_2 poly, const std::pair<Segment, Segment>& earliest_meeting_pair);

	static void updateVertices(std::vector<CGAL::Point_2<CGAL::Epick>>& vertices, const std::pair<Segment, Segment>& earliest_meeting_pair,
        const std::pair<Line, Line>& meeting_edges);

	void triangleMedialAxis(const std::vector<CGAL::Point_2<CGAL::Epick>>& vertices);

	std::list<Segment> m_medialAxisSegments;
};

std::optional<Point> get_intersection_point(const Line& l1, const Line& l2) {
    const CGAL::Object obj = intersection(l1, l2);
    if (const Point* p = CGAL::object_cast<Point>(&obj)) {
        return *p;
    }
    return std::nullopt;
}

Line get_angle_bisector(const Point& prev, const Point& curr, const Point& next) {
    // Construct vectors from current point to neighbors
    K::Vector_2 u1 = prev - curr;
    K::Vector_2 u2 = next - curr;

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

std::pair<Segment, Segment> MedialAxis::findNextSegmentPair(std::vector<CGAL::Point_2<CGAL::Epick>> vertices, std::pair<Line, Line>& meeting_edges, Point& center)
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

		auto inter = get_intersection_point(bisector1, bisector2);
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

void MedialAxis::addMedialAxisSegments(Polygon_2 poly, const std::pair<Segment, Segment>& earliest_meeting_pair)
{
    m_medialAxisSegments.emplace_back(earliest_meeting_pair.first);
    m_medialAxisSegments.emplace_back(earliest_meeting_pair.second);
}

void MedialAxis::updateVertices(std::vector<CGAL::Point_2<CGAL::Epick>>& vertices, const std::pair<Segment, Segment>& earliest_meeting_pair, const std::pair<Line, Line>& meeting_edges)
{
	vertices.push_back(get_intersection_point(meeting_edges.first, meeting_edges.second).value());
	vertices.erase(std::ranges::find(vertices, earliest_meeting_pair.first.point(0)));
	vertices.erase(std::ranges::find(vertices, earliest_meeting_pair.second.point(0)));
}

void MedialAxis::triangleMedialAxis(const std::vector<CGAL::Point_2<CGAL::Epick>>& vertices)
{
	const Point A = vertices[0], B = vertices[1], C = vertices[2];
	Point center = CGAL::centroid(A, B, C);
	m_medialAxisSegments.emplace_back(A, center);
	m_medialAxisSegments.emplace_back(B, center);
	m_medialAxisSegments.emplace_back(C, center);
}

MedialAxis::MedialAxis(const Polygon_2& pgn)
{
    if (!pgn.is_convex() || !pgn.is_simple()) {
        throw std::runtime_error("This algorithm only supports simple convex polygons.");
    }

    std::vector vertices(pgn.vertices_begin(), pgn.vertices_end());

    if (vertices.size() < 3)
    {
        throw std::runtime_error("Not enough vertices in chosen polygon.");
    }

    while (vertices.size() > 3) {
        std::pair<Line, Line> meeting_edges;
        Point center;

        auto next_meeting_bisectors= findNextSegmentPair(vertices, meeting_edges, center);
        addMedialAxisSegments(Polygon_2(vertices.begin(), vertices.end()), next_meeting_bisectors);
        updateVertices(vertices, next_meeting_bisectors, meeting_edges);
    }

    // Final triangle case
    triangleMedialAxis(vertices);
}

std::list<Segment> MedialAxis::getMedialAxis() const
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
        const auto segments = m.getMedialAxis();
        cout << "Medial Axis Edges: " << segments.size() << endl;

        for (const auto& seg : segments) {
            cout << seg.point(0) << " -> "
                << seg.point(1) << endl;
        }
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << endl;
        return 1;
    }

    return 0;
}

//
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
//            Recovering points using output iterator of type
//            Dispatch_or_drop_output_iterator
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
//            for (const auto& seg : m.getMedialAxis())
//            {
//                draw_in_ipe(seg);
//            }
//        };
//    }
//
//}
//
//CGAL_IPELET(medial_axis_ipelet::MedialAxisIpelet)
