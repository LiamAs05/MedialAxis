#include "GeometryUtils.hpp"

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