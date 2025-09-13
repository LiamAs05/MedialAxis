// --------------------------------------------------------------------
// Ipelet for computing the Medial Axis of a convex polygon
// --------------------------------------------------------------------

#include "ipelib.h"
#include "MedialAxis.hpp"
#include <vector>

using namespace ipe;
using vertices_t = std::vector<Point>;
using cgal_segments_t = std::list<CgalSegment>;

class MedialAxisIpelet : public Ipelet
{
public:
    MedialAxisIpelet() = default;
    virtual int ipelibVersion() const { return IPELIB_VERSION; }
    virtual bool run(int function, IpeletData *data, IpeletHelper *helper);

private:
    static Path *findSelectedPath(Page *page);
    static std::optional<vertices_t> extractVerticesFromShape(Shape path, IpeletHelper *helper);
    static Group *createSegmentGroup(cgal_segments_t segments, IpeletData *data);
};

Path *MedialAxisIpelet::findSelectedPath(Page *page)
{
  Path *polyPath = nullptr;
  for (int i = 0; i < page->count(); ++i)
  {
    if (page->select(i))
    {
      Object *obj = page->object(i);
      if (obj->type() == Object::EPath)
      {
        polyPath = obj->asPath();
        break;
      }
    }
  }
  return polyPath;
}

std::optional<vertices_t> MedialAxisIpelet::extractVerticesFromShape(Shape shape, IpeletHelper *helper)
{
    vertices_t vertices;
    const SubPath *sp = shape.subPath(0);
    const Curve *curve = sp->asCurve();
    
    if (!curve)
    {
        helper->messageBox("Unsupported path type", nullptr, 0);
        return std::nullopt;
    }

    for (int j = 0; j < curve->countSegments(); ++j)
    {
        Vector p = curve->segment(j).cp(0);
        vertices.emplace_back(p.x, p.y);
    }
    Vector last = curve->segment(curve->countSegments() - 1).last();
    vertices.emplace_back(last.x, last.y);

    if (vertices.size() < 3)
    {
        helper->messageBox("Polygon must have at least 3 vertices", nullptr, 0);
        return std::nullopt;
    }

    return vertices;
}

Group *MedialAxisIpelet::createSegmentGroup(cgal_segments_t segments, IpeletData *data)
{
    Group *group = new Group;
    for (const auto &s : segments)
    {
        Vector a(s.source().x(), s.source().y());
        Vector b(s.target().x(), s.target().y());
        group->push_back(new Path(data->iAttributes, Shape(Segment(a, b))));
    }
    return group;
}

bool MedialAxisIpelet::run(int function, IpeletData *data, IpeletHelper *helper)
{
    Page *page = data->iPage;
    if (!page)
    {
        helper->messageBox("No page open", nullptr, 0);
        return false;
    }

    auto polyPath = MedialAxisIpelet::findSelectedPath(page);
    if (!polyPath)
    {
        helper->messageBox("Please select a polygon path", nullptr, 0);
        return false;
    }

    Shape shape = polyPath->shape();
    if (shape.countSubPaths() == 0 || !shape.subPath(0)->closed())
    {
        helper->messageBox("Selected path is not a closed polygon", nullptr, 0);
        return false;
    }

    auto vertices = extractVerticesFromShape(shape, helper);
    if (!vertices.has_value())
    {
        return false;
    }

    try
    {
        Polygon_2 poly(vertices.value().begin(), vertices.value().end());
        MedialAxis m(poly);
        auto segments = m.get();
        auto segmentGroup = createSegmentGroup(segments, data);
        page->append(ESecondarySelected, data->iLayer, segmentGroup);
        // for (const auto &s : segments)
        // {
        //     Vector a(s.source().x(), s.source().y());
        //     Vector b(s.target().x(), s.target().y());
        //     page->append(ESecondarySelected, data->iLayer, new Path(data->iAttributes, Shape(Segment(a, b))));
        //     helper->messageBox("Added a segment", nullptr, 0);
        // }
    }
    catch (const std::exception &e)
    {
        helper->messageBox(e.what(), nullptr, 0);
        return false;
    }
    return true;
}

IPELET_DECLARE Ipelet *newIpelet()
{
  return new MedialAxisIpelet;
}
