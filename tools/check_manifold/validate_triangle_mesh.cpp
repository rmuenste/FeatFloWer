#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/version.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;

struct TriangleQualityMetrics {
    double area;
    double min_angle;
    double max_angle;
    double edge_ratio;
    double circum_inradius_ratio;
    bool is_degenerate;
    
    TriangleQualityMetrics() : area(0), min_angle(0), max_angle(0), 
                              edge_ratio(0), circum_inradius_ratio(0), 
                              is_degenerate(false) {}
};

class MeshValidator {
private:
    double area_threshold;
    double min_angle_threshold;
    double max_angle_threshold;
    double aspect_ratio_threshold;
    double circum_inradius_threshold;

public:
    MeshValidator(double area_tol, double min_ang_deg, double max_ang_deg, 
                  double edge_ratio, double rr_ratio) 
        : area_threshold(area_tol),
          min_angle_threshold(min_ang_deg * M_PI / 180.0),
          max_angle_threshold(max_ang_deg * M_PI / 180.0),
          aspect_ratio_threshold(edge_ratio),
          circum_inradius_threshold(rr_ratio) {}
    
    MeshValidator() 
        : area_threshold(1e-12),
          min_angle_threshold(1.0 * M_PI / 180.0),
          max_angle_threshold(179.0 * M_PI / 180.0),
          aspect_ratio_threshold(30.0),
          circum_inradius_threshold(10.0) {}
    bool load_mesh(const std::string& filename, Mesh& mesh) {
        if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh)) {
            std::cerr << "Error: Cannot read mesh from " << filename << std::endl;
            return false;
        }
        
        std::cout << "Loaded mesh with " << mesh.number_of_vertices() 
                  << " vertices and " << mesh.number_of_faces() << " faces" << std::endl;
        return true;
    }
    
    bool is_manifold_raw_check(const std::string& filename) {
        // Use direct CGAL::IO::read_OFF with Polyhedron_3 for raw manifold detection
        // This avoids any automatic repair that PMP might do
        
        std::ifstream in(filename);
        if (!in) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return false;
        }
        
        Polyhedron polyhedron;
        try {
            // Using direct OFF reader - no automatic repair
            if (CGAL::IO::read_OFF(in, polyhedron)) {
                in.close();
                
                // Check validity using Polyhedron's built-in validation
                if (polyhedron.is_valid()) {
                    return true;
                } else {
                    std::cout << "Non-manifold: Polyhedron validity check failed" << std::endl;
                    return false;
                }
            } else {
                std::cout << "Non-manifold: Could not parse OFF file (may contain non-manifold elements)" << std::endl;
                in.close();
                return false;
            }
        } catch (const std::exception& e) {
            std::cerr << "Exception during OFF file processing: " << e.what() << std::endl;
            std::cout << "Non-manifold: Exception during loading (likely non-manifold)" << std::endl;
            in.close();
            return false;
        }
    }
    
    bool is_manifold(const Mesh& mesh) {
        // Use CGAL's built-in functions for manifold checking on already loaded mesh
        // Note: This operates on the mesh that may have been auto-repaired during PMP loading
        
        // 1. Check overall mesh validity (includes edge manifoldness)
        if (!CGAL::is_valid_polygon_mesh(mesh)) {
            std::cout << "Non-manifold: Mesh fails CGAL validity check (edges/structure)" << std::endl;
            return false;
        }
        
        // 2. Check for non-manifold vertices using CGAL's specialized function
        namespace PMP = CGAL::Polygon_mesh_processing;
        for (auto v : mesh.vertices()) {
            if (PMP::is_non_manifold_vertex(v, mesh)) {
                std::cout << "Non-manifold: Vertex with non-manifold neighborhood found" << std::endl;
                return false;
            }
        }
        
        // 3. Check for isolated vertices
        for (auto v : mesh.vertices()) {
            if (mesh.is_isolated(v)) {
                std::cout << "Non-manifold: Isolated vertex found" << std::endl;
                return false;
            }
        }
        
        return true;
    }
    
    TriangleQualityMetrics analyze_triangle(const Point_3& p1, const Point_3& p2, const Point_3& p3) {
        TriangleQualityMetrics metrics;
        
        // Calculate edge vectors and lengths
        Vector_3 v21 = p2 - p1;
        Vector_3 v31 = p3 - p1;
        Vector_3 v32 = p3 - p2;
        
        double a = std::sqrt(CGAL::to_double(v32.squared_length())); // |p3 - p2|
        double b = std::sqrt(CGAL::to_double(v31.squared_length())); // |p3 - p1|
        double c = std::sqrt(CGAL::to_double(v21.squared_length())); // |p2 - p1|
        
        // Check for degenerate triangle (zero-length edges)
        if (a < 1e-12 || b < 1e-12 || c < 1e-12) {
            metrics.is_degenerate = true;
            return metrics;
        }
        
        // 1. Calculate triangle area using cross product
        Vector_3 cross = CGAL::cross_product(v21, v31);
        metrics.area = 0.5 * std::sqrt(CGAL::to_double(cross.squared_length()));
        
        // 2. Calculate angles using law of cosines
        double cos_alpha = (b*b + c*c - a*a) / (2*b*c);
        double cos_beta = (a*a + c*c - b*b) / (2*a*c);
        double cos_gamma = (a*a + b*b - c*c) / (2*a*b);
        
        // Clamp cosine values to [-1, 1] to avoid numerical errors
        cos_alpha = std::max(-1.0, std::min(1.0, cos_alpha));
        cos_beta = std::max(-1.0, std::min(1.0, cos_beta));
        cos_gamma = std::max(-1.0, std::min(1.0, cos_gamma));
        
        double alpha = std::acos(cos_alpha);
        double beta = std::acos(cos_beta);
        double gamma = std::acos(cos_gamma);
        
        metrics.min_angle = std::min({alpha, beta, gamma});
        metrics.max_angle = std::max({alpha, beta, gamma});
        
        // 3. Calculate aspect ratio (longest/shortest edge)
        double longest_edge = std::max({a, b, c});
        double shortest_edge = std::min({a, b, c});
        metrics.edge_ratio = longest_edge / shortest_edge;
        
        // 4. Calculate circumradius to inradius ratio
        if (metrics.area > 1e-12) {
            double semi_perimeter = (a + b + c) / 2.0;
            double inradius = metrics.area / semi_perimeter;
            double circumradius = (a * b * c) / (4.0 * metrics.area);
            metrics.circum_inradius_ratio = circumradius / inradius;
        } else {
            metrics.circum_inradius_ratio = std::numeric_limits<double>::infinity();
        }
        
        return metrics;
    }
    
    void validate_triangle_quality(const Mesh& mesh) {
        int total_triangles = 0;
        int degenerate_triangles = 0;
        int small_area_triangles = 0;
        int extreme_angle_triangles = 0;
        int high_aspect_ratio_triangles = 0;
        int high_circum_inradius_triangles = 0;
        
        std::cout << "\n=== Triangle Quality Analysis ===" << std::endl;
        std::cout << "Area threshold: " << area_threshold << std::endl;
        std::cout << "Min angle threshold: " << min_angle_threshold * 180.0 / M_PI << "°" << std::endl;
        std::cout << "Max angle threshold: " << max_angle_threshold * 180.0 / M_PI << "°" << std::endl;
        std::cout << "Aspect ratio threshold: " << aspect_ratio_threshold << std::endl;
        std::cout << "Circumradius/inradius threshold: " << circum_inradius_threshold << std::endl;
        
        for (auto f : mesh.faces()) {
            auto h = mesh.halfedge(f);
            auto v1 = mesh.target(h);
            auto v2 = mesh.target(mesh.next(h));
            auto v3 = mesh.target(mesh.next(mesh.next(h)));
            
            Point_3 p1 = mesh.point(v1);
            Point_3 p2 = mesh.point(v2);
            Point_3 p3 = mesh.point(v3);
            
            TriangleQualityMetrics metrics = analyze_triangle(p1, p2, p3);
            total_triangles++;
            
            if (metrics.is_degenerate) {
                degenerate_triangles++;
                std::cout << "Degenerate triangle found (face " << f << ")" << std::endl;
            }
            
            if (metrics.area < area_threshold) {
                small_area_triangles++;
                std::cout << "Small area triangle (face " << f << "): area = " << metrics.area << std::endl;
            }
            
            if (metrics.min_angle < min_angle_threshold || metrics.max_angle > max_angle_threshold) {
                extreme_angle_triangles++;
                std::cout << "Extreme angle triangle (face " << f << "): min = " 
                         << metrics.min_angle * 180.0 / M_PI << "°, max = " 
                         << metrics.max_angle * 180.0 / M_PI << "°" << std::endl;
            }
            
            if (metrics.edge_ratio > aspect_ratio_threshold) {
                high_aspect_ratio_triangles++;
                std::cout << "High aspect ratio triangle (face " << f << "): ratio = " 
                         << metrics.edge_ratio << std::endl;
            }
            
            if (metrics.circum_inradius_ratio > circum_inradius_threshold) {
                high_circum_inradius_triangles++;
                std::cout << "High circumradius/inradius triangle (face " << f << "): ratio = " 
                         << metrics.circum_inradius_ratio << std::endl;
            }
        }
        
        std::cout << "\n=== Summary ===" << std::endl;
        std::cout << "Total triangles: " << total_triangles << std::endl;
        std::cout << "Degenerate triangles: " << degenerate_triangles << std::endl;
        std::cout << "Small area triangles: " << small_area_triangles << std::endl;
        std::cout << "Extreme angle triangles: " << extreme_angle_triangles << std::endl;
        std::cout << "High aspect ratio triangles: " << high_aspect_ratio_triangles << std::endl;
        std::cout << "High circumradius/inradius triangles: " << high_circum_inradius_triangles << std::endl;
    }
};

int main(int argc, char* argv[]) {
    std::cout << "CGAL Triangle Mesh Validator" << std::endl;
    std::cout << "CGAL version: " << CGAL_VERSION_STR << std::endl;
    
    if (argc < 2) {
        std::cerr << "Usage  : " << argv[0]
                  << " input.off  [areaTol] [minAng] [maxAng] [edgeRatio] [RrRatio]\n"
                  << "Defaults: areaTol=1e-12  minAng=1°  maxAng=179°  edgeRatio=30  RrRatio=10\n";
        return EXIT_FAILURE;
    }

    // ── thresholds (user-overridable from command line) ──────────────────────
    const double AREA_EPS   = (argc > 2) ? std::stod(argv[2]) : 1e-12;
    const double MIN_ANG    = (argc > 3) ? std::stod(argv[3]) : 1.0;   // degrees
    const double MAX_ANG    = (argc > 4) ? std::stod(argv[4]) : 179.0; // degrees
    const double EDGE_RATIO = (argc > 5) ? std::stod(argv[5]) : 30.0;
    const double RR_RATIO   = (argc > 6) ? std::stod(argv[6]) : 10.0;
    
    std::string filename = argv[1];
    Mesh mesh;
    MeshValidator validator(AREA_EPS, MIN_ANG, MAX_ANG, EDGE_RATIO, RR_RATIO);
    
    // Load mesh
    if (!validator.load_mesh(filename, mesh)) {
        return 1;
    }
    
    // Check if mesh is manifold using raw check first (no auto-repair)
    std::cout << "\n=== Manifold Check ===" << std::endl;
    std::cout << "Raw manifold check (no auto-repair): ";
    bool is_manifold_raw = validator.is_manifold_raw_check(filename);
    std::cout << (is_manifold_raw ? "MANIFOLD" : "NON-MANIFOLD") << std::endl;
    
    // Also check the loaded mesh (may have been auto-repaired)
    std::cout << "Loaded mesh manifold check: ";
    bool is_manifold_loaded = validator.is_manifold(mesh);
    std::cout << (is_manifold_loaded ? "MANIFOLD" : "NON-MANIFOLD") << std::endl;
    
    // Final verdict - use the raw check as authoritative
    bool is_manifold = is_manifold_raw;
    std::cout << "Final verdict: Mesh is " << (is_manifold ? "MANIFOLD" : "NON-MANIFOLD") << std::endl;
    
    // Analyze triangle quality
    validator.validate_triangle_quality(mesh);
    
    return 0;
}