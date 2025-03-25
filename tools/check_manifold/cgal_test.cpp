#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/basic.h>

// Updated includes for CGAL 5.3+
// We have tested this so far for 5.3.2
// TODO:
// - Try GCAL 6.x
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/STL.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <cstddef>

namespace CGAL {

template <class HDS>
class Polyhedron_builder_from_STL : public CGAL::Modifier_base<HDS> {
  typedef typename HDS::Vertex::Point Point_3;
  typedef std::vector<std::array<double, 3>> Points_3;
  typedef std::array<int, 3> Facet;
  typedef std::vector<Facet> Surface;

  std::istream& is;
  Points_3 meshPoints;
  Surface mesh;

public:
  Polyhedron_builder_from_STL(std::istream& is_)
    : is(is_)
  {}

  void operator()(HDS& hds) {
    if (!CGAL::IO::read_STL(is, meshPoints, mesh)) return;

    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds);
    B.begin_surface(meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;

    for (size_type i = 0; i < meshPoints.size(); i++) {
      B.add_vertex(
        Point_3(meshPoints[i][0], meshPoints[i][1], meshPoints[i][2])
      );
    }

    for (size_type i = 0; i < mesh.size(); i++) {
      B.begin_facet();
      B.add_vertex_to_facet(mesh[i][0]);
      B.add_vertex_to_facet(mesh[i][1]);
      B.add_vertex_to_facet(mesh[i][2]);
      B.end_facet();
    }
    
    if (B.error()) {
      std::cerr << "An error occurred while creating a Polyhedron" << std::endl;
      B.rollback();
      is.clear(std::ios::badbit);
      return;
    }

    B.end_surface();
    is.clear(std::ios::goodbit);
  }
};

} // namespace CGAL

// Choose a geometry kernel
typedef CGAL::Simple_cartesian<double> Kernel;

// Make a short-hand for the geometry Kernel type
typedef Kernel::FT FT;
typedef Kernel::Point_3 cgalPoint;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Vector_3 Vec;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;

// Get the modifier for STL files
typedef CGAL::Polyhedron_builder_from_STL<HalfedgeDS> PolyhedronBuilder;

std::string getFileExtension(const std::string &fileName) {
  std::size_t i = fileName.rfind('.', fileName.length());

  std::string extension;
  if (i != std::string::npos) {
    extension = fileName.substr(i + 1, fileName.length() - 1);
  }
  else {
    std::cout << "Could not determine file type of: " << fileName << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "File extension: " << extension << std::endl;
  return extension;
}

void checkOffFile(const std::string &fileName) {
  // Open the file in a stream
  std::ifstream in(fileName);
  
  if (!in) {
    std::cerr << "Unable to open file" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  Polyhedron polyhedron;
  try {
    // Using the new IO functions in CGAL 5.3+
    if (CGAL::IO::read_OFF(in, polyhedron)) {
      if (polyhedron.is_valid(true, 0)) {
        std::cout << "CGAL mesh check: OK" << std::endl;
      } else {
        std::cout << "CGAL mesh check: FAIL (non-valid polyhedron)" << std::endl;
      }
    } else {
      std::cout << "CGAL mesh check: FAIL (could not read file)" << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "Exception during OFF file processing: " << e.what() << std::endl;
    std::cout << "CGAL mesh check: FAIL" << std::endl;
  }
  in.close();
}

void checkStlFile(const std::string &fileName) {
  // Open the file in a stream
  std::ifstream in(fileName);
  
  if (!in) {
    std::cerr << "Unable to open file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Polyhedron polyhedron;
  try {
    PolyhedronBuilder builder(in);
    polyhedron.delegate(builder);
    
    if (polyhedron.is_valid(true, 0) && in.good()) {
      std::cout << "CGAL mesh check: OK" << std::endl;
    } else {
      std::cout << "CGAL mesh check: FAIL" << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "Exception during STL file processing: " << e.what() << std::endl;
    std::cout << "CGAL mesh check: FAIL" << std::endl;
  }
  in.close();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "You need to supply the mesh file as a command line argument." << std::endl;
    std::cout << "Please call the program in the following way:" << std::endl;
    std::cout << "check_manifold myfilename.off or myfilename.stl" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "Name of the file: " << argv[1] << std::endl;

  std::string fileName(argv[1]);
  std::string extension = getFileExtension(fileName);

  if (extension == "off") {
    checkOffFile(fileName);
  }
  else if (extension == "stl") {
    checkStlFile(fileName);
  }
  else {
    std::cout << "Cannot handle file extension: " << extension << std::endl;
    std::exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}