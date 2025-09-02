// --------------------------------------------------------------------
// Ipelet for computing the Medial Axis of a convex polygon
// --------------------------------------------------------------------

#include "ipelib.h"
#include "MedialAxis.hpp"   // your C++ medial axis implementation
#include <vector>
#include <stdexcept>

using namespace ipe;

// --------------------------------------------------------------------

class MedialAxisIpelet : public Ipelet {
public:
  MedialAxisIpelet();
  virtual int ipelibVersion() const { return IPELIB_VERSION; }
  virtual bool run(int function, IpeletData *data, IpeletHelper *helper);
};

MedialAxisIpelet::MedialAxisIpelet()
{
  // nothing needed
}

bool MedialAxisIpelet::run(int function, IpeletData *data, IpeletHelper *helper)
{
  Page *page = data->iPage;
  if (!page) {
    helper->messageBox("No page open", nullptr, 0);
    return false;
  }

  // Collect a selected polygon path
  Path *polyPath = nullptr;
  for (int i = 0; i < page->count(); ++i) {
    if (page->select(i)) {
      Object *obj = page->object(i);
      if (obj->type() == Object::EPath) {
        polyPath = obj->asPath();
        break;
      }
    }
  }

  if (!polyPath) {
    helper->messageBox("Please select a polygon path", nullptr, 0);
    return false;
  }

  Shape shape = polyPath->shape();
  if (shape.countSubPaths() == 0 || !shape.subPath(0)->closed()) {
    helper->messageBox("Selected path is not a closed polygon", nullptr, 0);
    return false;
  }

  // Extract polygon vertices
  std::vector<Point> vertices;
  const SubPath *sp = shape.subPath(0);
  const Curve *curve = sp->asCurve();
  if (!curve) {
    helper->messageBox("Unsupported path type", nullptr, 0);
    return false;
  }
  for (int j = 0; j < curve->countSegments(); ++j) {
    Vector p = curve->segment(j).cp(0);
    vertices.emplace_back(p.x, p.y);
  }
  Vector last = curve->segment(curve->countSegments() - 1).last();
  vertices.emplace_back(last.x, last.y);

  if (vertices.size() < 3) {
    helper->messageBox("Polygon must have at least 3 vertices", nullptr, 0);
    return false;
  }

  try {
    Polygon_2 poly(vertices.begin(), vertices.end());
    MedialAxis m(poly);
    auto segs = m.get();

    // Build group of medial axis segments
    Group *group = new Group;
    for (const auto &s : segs) {
      Vector a(s.source().x(), s.source().y());
      Vector b(s.target().x(), s.target().y());
      group->push_back(new Path(data->iAttributes, Shape(Segment(a, b))));
    }

    page->append(ESecondarySelected, data->iLayer, group);
  }
  catch (const std::exception &e) {
    helper->messageBox(e.what(), nullptr, 0);
    return false;
  }

  return true;
}

// --------------------------------------------------------------------

IPELET_DECLARE Ipelet *newIpelet()
{
  return new MedialAxisIpelet;
}

