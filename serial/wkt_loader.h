#include <fstream>
#include <boost/geometry.hpp>

#include "File1.cu"


Points* LoadPoints(const std::string& path) {
  std::ifstream ifs(path);
  std::string line;
  using boost_point_t =boost::geometry::model::point<COORD_T, N_DIMS, boost::geometry::cs::cartesian>;
  using boost_polygon_t = boost::geometry::model::polygon<boost_point_t>;
  using point_t = typename cuda_vec<COORD_T, N_DIMS>::type;
  std::vector<point_t> points;

  while (std::getline(ifs, line)) {
    if (!line.empty()) {
      if (line.rfind("MULTIPOLYGON", 0) == 0) {
        boost::geometry::model::multi_polygon<boost_polygon_t> multi_poly;
        boost::geometry::read_wkt(line, multi_poly);

        for (auto& poly : multi_poly) {
          for (auto& p : poly.outer()) {
            points.push_back(*reinterpret_cast<point_t*>(&p));
          }
        }
      } else if (line.rfind("POLYGON", 0) == 0) {
        boost_polygon_t poly;
        boost::geometry::read_wkt(line, poly);

        for (auto& p : poly.outer()) {
          points.push_back(*reinterpret_cast<point_t*>(&p));
        }
      } else if (line.rfind("POINT", 0) == 0) {
        boost_point_t p;
        boost::geometry::read_wkt(line, p);

        points.push_back(*reinterpret_cast<point_t*>(&p));
      } else {
        std::cerr << "Bad Geometry " << line << "\n";
        abort();
      }
      if (points.size() % 1000 == 0) {
        std::cout << "Loaded geometries " << points.size() / 1000 << std::endl;
      }
      if (points.size() >= limit) {
        break;
      }
    }
  }
  ifs.close();
  return points;
}
