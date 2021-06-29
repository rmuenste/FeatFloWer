#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/config.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/IO/Polyhedron_scan_OFF.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <fstream>
#include <string>


#include <CGAL/license/Polyhedron.h>


#include <CGAL/basic.h>
#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/STL_reader.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <iostream>
#include <cstddef>

namespace CGAL {

template <class HDS>
class Polyhedron_builder_from_STL : public CGAL::Modifier_base<HDS> {
  typedef typename HDS::Vertex::Point Point_3;
  typedef std::vector<cpp11::array<double, 3> > Points_3;
  typedef cpp11::array<int,3> Facet;
  typedef std::vector<Facet> Surface;

  std::istream& is;
  Points_3 meshPoints;
  Surface mesh;

public:

  Polyhedron_builder_from_STL(std::istream& is_)
    : is(is_)
  {}

  void operator()( HDS& hds) {
    if(!read_STL(is, meshPoints, mesh)) return;

    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds);
    B.begin_surface( meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;

    for(size_type i=0; i < meshPoints.size(); i++){
      B.add_vertex(
        Point_3(meshPoints[i][0], meshPoints[i][1], meshPoints[i][2])
      );
    }

    for(size_type i=0; i < mesh.size(); i++){
      B.begin_facet();
      B.add_vertex_to_facet( mesh[i][0]);
      B.add_vertex_to_facet( mesh[i][1]);
      B.add_vertex_to_facet( mesh[i][2]);
      B.end_facet();
    }
    if(B.error())
      {
        std::cerr << "An error occured while creating a Polyhedron" << std::endl;
        B.rollback();
        is.clear( std::ios::badbit);
        return;
      }

    B.end_surface();
    is.clear( std::ios::goodbit);
  }
};

template < class HDS>
class Polyhedron_builder_from_OFF :  public Modifier_base<HDS> {
protected:
    std::istream&    m_in;
    File_header_OFF  m_file_header;
public:

    typedef HDS Halfedge_data_structure;

// DEFINITION
//
// Polyhedron_builder_from_OFF<Traits> is a polyhedral surface builder.
// It scans a polyhedron given in OFF from a stream and appends it
// incrementally using the incremental builder.

    Polyhedron_builder_from_OFF( std::istream& in, bool verbose = false)
        : m_in(in), m_file_header( verbose) {}

    // Activation
    void operator()( HDS& hds);

    const File_header_OFF&  header() const { return m_file_header; }
};

template < class HDS >
void
Polyhedron_builder_from_OFF<HDS>:: operator()( HDS& target) {
    File_scanner_OFF scanner( m_in, m_file_header.verbose());
    if ( ! m_in) {
        if ( scanner.verbose()) {
            std::cerr << " " << std::endl;
            std::cerr << "Polyhedron_builder_from_OFF<HDS>::" << std::endl;
            std::cerr << "operator(): input error: file format is not in "
                         "OFF." << std::endl;
        }
        return;
    }
    m_file_header = scanner;  // Remember file header after return.

    Polyhedron_incremental_builder_3<HDS> B( target, scanner.verbose());
    B.begin_surface( scanner.size_of_vertices(),
                     scanner.size_of_facets(),
                     scanner.size_of_halfedges());

    typedef typename HDS::Traits     Traits;
    typedef typename Traits::Point_3 Point;

    // read in all vertices
    std::size_t  i;
    for ( i = 0; i < scanner.size_of_vertices(); i++) {
        Point p;
        file_scan_vertex( scanner, p);
        B.add_vertex( p);
        if(scanner.has_colors())
        {
         Color c;
         file_scan_color(scanner, c);
        }
        else
         scanner.skip_to_next_vertex( i);
    }
    if ( ! m_in  || B.error()) {
        B.rollback();
        m_in.clear( std::ios::badbit);
        return;
    }

    // read in all facets
    for ( i = 0; i < scanner.size_of_facets(); i++) {
        B.begin_facet();
        std::size_t no;
        scanner.scan_facet( no, i);
        if( ! m_in || B.error() || no < 3) {
            if ( scanner.verbose()) {
                std::cerr << " " << std::endl;
                std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
                std::cerr << "operator()(): input error: facet " << i
                     << " has less than 3 vertices." << std::endl;
            }
            B.rollback();
            m_in.clear( std::ios::badbit);
            return;
        }
        for ( std::size_t j = 0; j < no; j++) {
            std::size_t index;
            scanner.scan_facet_vertex_index( index, i);
            B.add_vertex_to_facet( index);
        }
        //TO DO : Insert read color
        B.end_facet();
        scanner.skip_to_next_facet( i);
    }
    if ( ! m_in  || B.error()) {
        B.rollback();
        m_in.clear( std::ios::badbit);
    }
    if ( B.check_unconnected_vertices()) {
        if ( ! B.remove_unconnected_vertices()) {
            if ( scanner.verbose()) {
                std::cerr << " " << std::endl;
                std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
                std::cerr << "operator()(): input error: cannot "
                             "successfully remove isolated vertices."
                          << std::endl;
            }
            B.rollback();
            m_in.clear( std::ios::badbit);
        }
    }
    B.end_surface();
}

} //namespace CGAL




// Choose a geometry kernel
typedef CGAL::Simple_cartesian<double> Kernel;

// Make a short-hand for the geometry Kernel type
typedef Kernel::FT FT;

typedef Kernel::Point_3 cgalPoint;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Vector_3 Vec;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS     HalfedgeDS;

// Get the modifier for OFF files
typedef CGAL::Polyhedron_builder_from_OFF<HalfedgeDS> PolyhedronBuilderOff;

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
  
  if (!in)
  {
    std::cerr << "unable to open file" << std::endl;
    std::exit(EXIT_FAILURE);
  }
	
  Polyhedron polyhedron;
  PolyhedronBuilderOff scanner(in, true);
  polyhedron.delegate(scanner);
  CGAL_assertion(polyhedron.is_valid(true, 0));
  if (in.good())
   std::cout << "CGAL mesh check: OK" << std::endl;
  else
   std::cout << "CGAL mesh check: FAIL" << std::endl;
  in.close();

}

void checkStlFile(const std::string &fileName) {
  // Open the file in a stream
  std::ifstream in(fileName);
  
  if (!in)
  {
    std::cerr << "unable to open file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Polyhedron polyhedron;
  PolyhedronBuilder builder(in);
  polyhedron.delegate(builder);
  CGAL_assertion(polyhedron.is_valid(true, 0));
  if (in.good())
   std::cout << "CGAL mesh check: OK" << std::endl;
  else
   std::cout << "CGAL mesh check: FAIL" << std::endl;
  in.close();
}

int main(int argc, char * argv[]) {
  
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
