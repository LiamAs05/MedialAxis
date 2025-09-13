#include "GeometryUtils.hpp"

bool pointsEqualEps(const Point& a, const Point& b, double eps) {
    return (std::abs(CGAL::to_double(a.x()) - CGAL::to_double(b.x())) <= eps) &&
           (std::abs(CGAL::to_double(a.y()) - CGAL::to_double(b.y())) <= eps);
}

std::optional<CgalSegment> clipCgalSegmentToPolygon(const CgalSegment& seg, const Polygon_2& poly) {
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

void removePointIfExists(std::vector<Point>& vertices, const Point& toRemove)
{
    auto it = std::find(vertices.begin(), vertices.end(), toRemove);
    if (it != vertices.end())
    {
        vertices.erase(it);
    }
}

std::optional<Point> getIntersectionPoint(const CgalLine& l1, const CgalLine& l2) {
    const CGAL::Object obj = intersection(l1, l2);
    if (const Point* p = CGAL::object_cast<Point>(&obj)) {
        return *p;
    }
    return std::nullopt;
}

CgalLine computeAngleBisector(const Point& prev, const Point& curr, const Point& next) {
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

Point computeTriangleIncenter(const Point& A, const Point& B, const Point& C)
{
    double a = std::sqrt(CGAL::squared_distance(B, C)); // length opposite A
    double b = std::sqrt(CGAL::squared_distance(C, A)); // length opposite B
    double c = std::sqrt(CGAL::squared_distance(A, B)); // length opposite C

    double px = a * A.x() + b * B.x() + c * C.x();
    double py = a * A.y() + b * B.y() + c * C.y();
    double denom = a + b + c;

    return Point(px / denom, py / denom);
}