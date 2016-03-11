#include "loader.h"

#include <fstream>
#include <CGAL/Gmpq.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
// #include <CGAL/Cartesian_converter.h>

typedef CGAL::Extended_cartesian< CGAL::Lazy_exact_nt<CGAL::Gmpq> > Kernel;
// typedef CGAL::Simple_cartesian<double>                           DoubleKernel;
// typedef CGAL::Cartesian_converter<Kernel, DoubleKernel>                         ConverterToDouble;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron_3;
typedef Nef_polyhedron_3::Plane_3  Plane_3;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Vertex_handle Vertex_handle;

// http://jamesgregson.blogspot.com/2012/05/example-code-for-building.html
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
    std::vector<GLfloat> &coords;
    std::vector<GLuint>    &tris;
    polyhedron_builder( std::vector<GLfloat> &_coords, std::vector<GLuint> &_tris ) : coords(_coords), tris(_tris) {}
    void operator()( HDS& hds) {
        typedef typename HDS::Vertex Vertex;
        typedef typename Vertex::Point Point;

    // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( coords.size()/3, tris.size()/3 );

        // add the polyhedron vertices
        for( int i=0; i<(int)coords.size(); i+=3 ){
            B.add_vertex( Point( coords[i+0], coords[i+1], coords[i+2] ) );
        }

        // add the polyhedron triangles
        for( int i=0; i<(int)tris.size(); i+=3 ){
            B.begin_facet();
            B.add_vertex_to_facet( tris[i+0] );
            B.add_vertex_to_facet( tris[i+1] );
            B.add_vertex_to_facet( tris[i+2] );
            B.end_facet();
        }

        // finish up the surface
        B.end_surface();
    }
};

Loader::Loader(QObject* parent, const QString& filename)
    : QThread(parent), filename(filename)
{
    // Nothing to do here
}

void Loader::run()
{
    Mesh* mesh = load_stl();
    if (mesh)
    {
        emit got_mesh(mesh);
        emit loaded_file(filename);
    }
}


////////////////////////////////////////////////////////////////////////////////

struct Vec3
{
    GLfloat x, y, z;
    bool operator!=(const Vec3& rhs) const
    {
        return x != rhs.x || y != rhs.y || z != rhs.z;
    }
    bool operator<(const Vec3& rhs) const
    {
        if      (x != rhs.x)    return x < rhs.x;
        else if (y != rhs.y)    return y < rhs.y;
        else if (z != rhs.z)    return z < rhs.z;
        else                    return false;
    }
};

typedef std::pair<Vec3, GLuint> Vec3i;

////////////////////////////////////////////////////////////////////////////////

Mesh* Loader::load_stl()
{
    QFile file(filename);
    file.open(QIODevice::ReadOnly);
    if (file.read(5) == "solid")
    {
        emit error_ascii_stl();
        return NULL;
    }
    // Skip the rest of the header material
    file.read(75);

    QDataStream data(&file);
    data.setByteOrder(QDataStream::LittleEndian);
    data.setFloatingPointPrecision(QDataStream::SinglePrecision);

    // Load the triangle count from the .stl file
    uint32_t tri_count;
    data >> tri_count;

    // Verify that the file is the right size
    if (file.size() != 84 + tri_count*50)
    {
        emit error_bad_stl();
        return NULL;
    }

    // Extract vertices into an array of xyz, unsigned pairs
    QVector<Vec3i> verts(tri_count*3);

    // Dummy array, because readRawData is faster than skipRawData
    char buffer[sizeof(float)*3];

    // Store vertices in the array, processing one triangle at a time.
    for (auto v=verts.begin(); v != verts.end(); v += 3)
    {
        // Skip face's normal vector
        data.readRawData(buffer, 3*sizeof(float));

        // Load vertex data from .stl file into vertices
        data >> v[0].first.x >> v[0].first.y >> v[0].first.z;
        data >> v[1].first.x >> v[1].first.y >> v[1].first.z;
        data >> v[2].first.x >> v[2].first.y >> v[2].first.z;

        // Skip face attribute
        data.readRawData(buffer, sizeof(uint16_t));
    }

    // Save indicies as the second element in the array
    // (so that we can reconstruct triangle order after sorting)
    for (size_t i=0; i < tri_count*3; ++i)
    {
        verts[i].second = i;
    }

    // Sort the set of vertices (to deduplicate)
    std::sort(verts.begin(), verts.end());

    // This vector will store triangles as sets of 3 indices
    std::vector<GLuint> indices(tri_count*3);

    // Go through the sorted vertex list, deduplicating and creating
    // an indexed geometry representation for the triangles.
    // Unique vertices are moved so that they occupy the first vertex_count
    // positions in the verts array.
    size_t vertex_count = 0;
    for (auto v : verts)
    {
        if (!vertex_count || v.first != verts[vertex_count-1].first)
        {
            verts[vertex_count++] = v;
        }
        indices[v.second] = vertex_count - 1;
    }
    verts.resize(vertex_count);

    std::vector<GLfloat> flat_verts;
    flat_verts.reserve(vertex_count*3);
    for (auto v : verts)
    {
        flat_verts.push_back(v.first.x);
        flat_verts.push_back(v.first.y);
        flat_verts.push_back(v.first.z);
    }

    Polyhedron P_init, P_final;
    polyhedron_builder<HalfedgeDS> builder(flat_verts, indices );
    P_init.delegate( builder );

    if (P_init.is_closed()) {
        Nef_polyhedron_3 Np (P_init);
        Nef_polyhedron_3 Plane (Plane_3(0,0,1,0));
        Nef_polyhedron_3 Intersect = Np * Plane;

        if (Intersect.is_simple()) {
            Intersect.convert_to_polyhedron(P_final);
            // write the polyhedron out as a .OFF file
            std::ofstream os("dump.off");
            os << P_final;
            os.close();

            // Container holding last line read
            std::string readLine;
            // Open file for reading
            std::ifstream in("dump.off");

            // Check if file is in OFF format
            getline(in,readLine);
            if (readLine == "OFF") {
                int MAX_BUFFER_SIZE = 1000;
                char buffer[MAX_BUFFER_SIZE];

                // Read values for Nv and Nf
                int nv, nf;
                in.getline(buffer, MAX_BUFFER_SIZE);
                std::stringstream ss(buffer);
                ss >> nv >> nf;

                // skip blank line
                getline(in,readLine);

                // read the vertices
                std::vector<GLfloat> new_verts;
                new_verts.reserve(nv);
                for (int n=0; n<nv; n++) {
                   in.getline(buffer, MAX_BUFFER_SIZE);
                   std::stringstream ss(buffer);
                   GLfloat x, y, z;
                   ss >> x >> y >> z;
                   new_verts.push_back(x);
                   new_verts.push_back(y);
                   new_verts.push_back(z);
                }

                // read the faces
                std::vector<GLuint> new_faces;
                new_faces.reserve(nf * 3);
                for (int n=0; n<nf; n++) {
                   in.getline(buffer, MAX_BUFFER_SIZE);
                   std::stringstream ss(buffer);
                   GLuint f0, f1, f2, f3;
                   ss >> f0 >> f1 >> f2 >> f3;
                   new_faces.push_back(f1);
                   new_faces.push_back(f2);
                   new_faces.push_back(f3);
                }

                return new Mesh(new_verts, new_faces);
            } else {
                std::cout << "The file to read is not in OFF format." << std::endl;
            }

        } else {
            std::cout << "Error: Nef polyhedron is not simple." << std::endl;
        }
    } else {
        std::cout << "Error: Initial mesh is not closed." << std::endl;
    }
    return new Mesh(flat_verts, indices);
}

