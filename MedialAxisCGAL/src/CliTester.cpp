#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/WKT.h>

#include "MedialAxis.hpp"

// Aliases matching your header
using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point     = K::Point_2;
using Segment   = K::Segment_2;
using Polygon_2 = CGAL::Polygon_2<K>;

// For WKT output (CGAL treats a LineString as a std::vector<Point>)
using LineString       = std::vector<Point>;
using MultiLineString  = std::vector<LineString>;

int main(int argc, char** argv) {
  try {
    Polygon_2 poly;

    // Read POLYGON WKT from file or stdin
    if (argc > 1) {
      std::ifstream in(argv[1]);
      if (!in) {
        std::cerr << "Error: cannot open file: " << argv[1] << "\n";
        return 1;
      }
      if (!CGAL::IO::read_polygon_WKT(in, poly)) {
        std::cerr << "Error: failed to read POLYGON WKT from file.\n";
        return 1;
      }
    } else {
      if (!CGAL::IO::read_polygon_WKT(std::cin, poly)) {
        std::cerr << "Error: failed to read POLYGON WKT from stdin.\n";
        return 1;
      }
    }

    // Compute medial axis
    MedialAxis m(poly);
    std::list<Segment> segs = m.get();

    // Collect into MultiLineString
    MultiLineString mls;
    mls.reserve(segs.size());
    for (const auto& s : segs) {
      LineString ls;
      ls.push_back(s.point(0));
      ls.push_back(s.point(1));
      mls.push_back(std::move(ls));
    }

    // Output WKT
    CGAL::IO::write_multi_linestring_WKT(std::cout, mls);
    std::cout << "\n";
    return 0;

  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}

