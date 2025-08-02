#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/WKT.h>
#include <iostream>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Line_2 Line;

using std::cout;
using std::endl;

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

std::optional<Point> get_intersection_point(const Line& l1, const Line& l2) {
    const CGAL::Object obj = CGAL::intersection(l1, l2);
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

std::vector<std::pair<Point, Point>> compute_medial_axis(const Polygon_2& poly) {
    std::vector<std::pair<Point, Point>> result;
    std::vector vertices(poly.vertices_begin(), poly.vertices_end());
    std::size_t n = vertices.size();

    if (n < 3)
    {
        throw std::runtime_error("Not enough vertices in chosen polygon.");
    }

    while (n > 3) {
        double min_radius = std::numeric_limits<double>::max();
        std::pair<Point, Point> earliest_meeting_pair;
        std::pair<Line, Line> meeting_edges;
        Point center;

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
                earliest_meeting_pair = { vertices[i1], vertices[i2] };
                meeting_edges = { Line(vertices[i0], vertices[i1]), Line(vertices[i2], vertices[i3]) };
                center = *inter;
            }
        }

        result.emplace_back(earliest_meeting_pair.first, center);
        result.emplace_back(earliest_meeting_pair.second, center);
        vertices.push_back(get_intersection_point(meeting_edges.first, meeting_edges.second).value());

        vertices.erase(std::ranges::find(vertices, earliest_meeting_pair.first));
        auto second_idx = std::ranges::find(vertices, earliest_meeting_pair.second);
        vertices.erase(second_idx);

        n = vertices.size();
    }

    // Final triangle case
    const Point A = vertices[0], B = vertices[1], C = vertices[2];
    Point center = CGAL::centroid(A, B, C);
    result.emplace_back(A, center);
    result.emplace_back(B, center);
    result.emplace_back(C, center);

    return result;
}

int main()
{
    try {
        const Polygon_2 pgn = polygon_from_wkt();

        cout << "The polygon is " << (pgn.is_simple() ? "" : "not ") << "simple." << endl;
        cout << "The polygon is " << (pgn.is_convex() ? "" : "not ") << "convex." << endl;

        if (!pgn.is_convex() || !pgn.is_simple()) {
            throw std::runtime_error("This algorithm only supports simple convex polygons.");
        }

        const auto segments = compute_medial_axis(pgn);
        cout << "Medial Axis Edges: " << segments.size() << endl;

        for (const auto& [fst, snd] : segments) {
            cout << fst.x() << " " << fst.y() << " -> "
                << snd.x() << " " << snd.y() << endl;
        }
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << endl;
        return 1;
    }

    return 0;
}
