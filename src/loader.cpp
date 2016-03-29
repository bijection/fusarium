#include "loader.h"

#include <igl/readSTL.h>

#define IGL_NO_CORK

#include <fstream>
#include <iostream>

#include <igl/copyleft/boolean/mesh_boolean.h>
#include <Eigen/Core>

#include <QtDebug>
#include <string>

#include "clipper.hpp"


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


void clipper_test(){
    using namespace ClipperLib;

    Path subj;
    Paths solution;
    
    subj << IntPoint(348,257) << IntPoint(364,148) << IntPoint(362,148) << 
        IntPoint(326,241) << IntPoint(295,219) << IntPoint(258,88) << 
        IntPoint(440,129) << IntPoint(370,196) << IntPoint(372,275);
    
    ClipperOffset co;
    co.AddPath(subj, jtRound, etClosedPolygon);
    co.Execute(solution, -7.0);

    for (size_t i = 0; i < solution.size(); ++i) {
    for (size_t j = 0; j < solution[i].size(); ++j) {
        qDebug() << solution[i][j].X << " " << solution[i][j].Y;
    }}
}


Mesh* Loader::load_stl()
{
    Eigen::MatrixXd VA,VB,VC,NA;
    Eigen::MatrixXi FA,FB,FC;

    igl::readSTL(filename.toUtf8().constData(), VA, FA, NA);

    // igl::readSTL("../gl/cube.stl", VB, FB, NA);
    clipper_test();


    VB= (Eigen::MatrixXd(8,3)<<
        0,0,0,
        0,0,1,
        0,1,0,
        0,1,1,
        1,0,0,
        1,0,1,
        1,1,0,
        1,1,1).finished().array()-0.5;

    FB = (Eigen::MatrixXi(12,3)<<
        1,7,5,
        1,3,7,
        1,4,3,
        1,2,4,
        3,8,7,
        3,4,8,
        5,7,8,
        5,8,6,
        1,5,6,
        1,6,2,
        2,6,8,
        2,8,4).finished().array() - 1;

    igl::copyleft::boolean::mesh_boolean(VA,FA,VB,FB,
        igl::copyleft::boolean::MESH_BOOLEAN_TYPE_UNION,
        VC,FC);

    std::vector<GLfloat> flat_verts;
    std::vector<GLuint> indices;

    for (size_t i = 0; i < VC.rows(); ++i)
    {
        flat_verts.push_back(VC(i,0));
        flat_verts.push_back(VC(i,1));
        flat_verts.push_back(VC(i,2));
    }

    for (size_t i = 0; i < FC.rows(); ++i)
    {
        indices.push_back(FC(i,0));
        indices.push_back(FC(i,1));
        indices.push_back(FC(i,2));
    }

    return new Mesh(flat_verts, indices);
}

