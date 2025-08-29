#include "MedialAxis.hpp"

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

int main() { 
	try { 
		const Polygon_2 pgn = polygon_from_wkt();
		MedialAxis m(pgn);
		const auto segments = m.get();
		std::cout << "Medial Axis Edges: " << segments.size() << std::endl;
        for (const auto& seg : segments)
        {
            std::cout << seg.point(0) << " -> " << seg.point(1) << std::endl;
        }
    } catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}