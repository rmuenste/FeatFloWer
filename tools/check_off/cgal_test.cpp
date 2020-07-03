#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/config.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <fstream>

// Choose a geometry kernel
typedef CGAL::Simple_cartesian<double> Kernel;

// Make a short-hand for the geometry Kernel type
typedef Kernel::FT FT;

typedef Kernel::Point_3 cgalPoint;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Vector_3 Vec;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char * argv[]) {
  
  if (argc < 2) {
    std::cout << "You need to supply the mesh file as a command line argument." << std::endl;
    std::cout << "Please call the program in the following way:" << std::endl;
    std::cout << "check_off myfilename.off" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "Name of off file: " << argv[1] << std::endl;

  // Load a mesh from file in the CGAL off format
  std::ifstream in(argv[1]);

  if (!in)
  {
    std::cerr << "unable to open file" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  Polyhedron polyhedron;
  // Read the polyhedron from the stream
  in >> polyhedron;

  if (!in)
  {
    std::cerr << "invalid OFF file" << std::endl;
    in.close();
    std::exit(EXIT_FAILURE);
  }

  in.close();
  std::cout << "OFF file loaded successfully" << std::endl;

  return EXIT_SUCCESS;
}
