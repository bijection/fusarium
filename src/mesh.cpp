#include <QFile>
#include <QDataStream>
#include <QVector3D>

#include <iostream>
#include <cmath>
#include <limits>
#include "../vecmath/vecmath.h"
#include "mesh.h"
#include "clipper.hpp"

#include <igl/copyleft/boolean/mesh_boolean.h>
#include <Eigen/Core>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polygon_2.h>

struct vertexInfo
{
  GLuint index;
};

struct FaceInfo2
{
  FaceInfo2(){}
  int nesting_level;
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<vertexInfo, K> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                       Tds;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag>  CDT;
typedef K::Point_2                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;
////////////////////////////////////////////////////////////////////////////////

Mesh::Mesh(std::vector<GLfloat> v, std::vector<GLuint> i)
    : vertices(v), indices(i)
{
    // Nothing to do here
}

void Mesh::setTransform(float rotateX, float rotateY, float rotateZ)
{
    QVector3D v0, v1;
    QMatrix4x4 transform = QMatrix4x4();
    transform.rotate(rotateX, QVector3D(1, 0, 0));
    transform.rotate(rotateY,  QVector3D(0, 1, 0));
    transform.rotate(rotateZ,  QVector3D(0, 0, 1));
    for (size_t i=0; i<vertices.size(); i+= 3) {
        v0 = QVector3D(vertices[i], vertices[i+1], vertices[i+2]);
        v1 = transform * v0;
        if (i == 0) {
            bbox.xmin = v1.x();
            bbox.ymin = v1.y();
            bbox.zmin = v1.z();
            bbox.xmax = v1.x();
            bbox.ymax = v1.y();
            bbox.zmax = v1.z();
        } else {
            bbox.xmin = fmin(bbox.xmin, v1.x());
            bbox.ymin = fmin(bbox.ymin, v1.y());
            bbox.zmin = fmin(bbox.zmin, v1.z());
            bbox.xmax = fmax(bbox.xmax, v1.x());
            bbox.ymax = fmax(bbox.ymax, v1.y());
            bbox.zmax = fmax(bbox.zmax, v1.z());
        }
    }
}

float Mesh::calculateProjectedArea(QVector3D norm) {
    QVector3D x,y,z,p,q,r;
    float norm2 = norm.length();
    float totalArea = 0;
    for (size_t i=0; i<indices.size(); i+= 3) {
        GLuint xIdx = 3 * indices[i];
        GLuint yIdx = 3 * indices[i+1];
        GLuint zIdx = 3 * indices[i+2];
        // x y z are the coordinates of the triangle
        x = QVector3D(vertices[xIdx], vertices[xIdx+1], vertices[xIdx+2]);
        y = QVector3D(vertices[yIdx], vertices[yIdx+1], vertices[yIdx+2]);
        z = QVector3D(vertices[zIdx], vertices[zIdx+1], vertices[zIdx+2]);

        // p q r are the coordinates of the triangle after rotation and projected onto the xz plane
        p = x - QVector3D::dotProduct(x, norm) / norm2 * norm;
        q = y - QVector3D::dotProduct(y, norm) / norm2 * norm;
        r = z - QVector3D::dotProduct(z, norm) / norm2 * norm;

        float area = QVector3D::crossProduct((p-q), (p-r)).length() / 2;
        totalArea += area;
    }
    return totalArea;
}

bool Mesh::checkForOverhangs(QVector3D dir, QVector3D &qStart, QVector3D &qEnd) const {

    Vector3f direction = Vector3f(dir.x(), dir.y(), dir.z());
    for (GLuint i = 0; i < vertices.size(); i += 3) {
        Vector3f origin = Vector3f(vertices[i], vertices[i+1], vertices[i+2]);
        int hits = 0;
        float minT = -1;
        float maxT = 1;
        for (GLuint j = 0; j < indices.size(); j += 3) {
            if (i != 3*indices[j] && i != 3*indices[j+1] && i != 3*indices[j+2]) {
                Vector3f triA = Vector3f(vertices[3*indices[j]], vertices[3*indices[j]+1], vertices[3*indices[j]+2]);
                Vector3f triB = Vector3f(vertices[3*indices[j+1]], vertices[3*indices[j+1]+1], vertices[3*indices[j+1]+2]);
                Vector3f triC = Vector3f(vertices[3*indices[j+2]], vertices[3*indices[j+2]+1], vertices[3*indices[j+2]+2]);

                Matrix3f A = Matrix3f();
                A.setCol(0, triA - triB);
                A.setCol(1, triA - triC);
                A.setCol(2, direction);
                float detA = A.determinant();

                Matrix3f betaMatrix = Matrix3f(A);
                betaMatrix.setCol(0, triA - origin);
                float beta = betaMatrix.determinant() / detA;

                Matrix3f gammaMatrix = Matrix3f(A);
                gammaMatrix.setCol(1, triA - origin);
                float gamma = gammaMatrix.determinant() / detA;

                Matrix3f tMatrix = Matrix3f(A);
                tMatrix.setCol(2, triA - origin);
                float t = tMatrix.determinant() / detA;

                if ((beta + gamma) < 1 && beta > 0 && gamma > 0 && std::abs(t) >= 1e-16) {
                    hits += 1;
                    minT = fmin(minT, t);
                    maxT = fmax(maxT, t);
                    if (hits > 1) {
                        Vector3f start = origin + (minT-1) * direction;
                        Vector3f end = origin + (maxT+1) * direction;
                        qStart.setX(start.x());
                        qStart.setY(start.y());
                        qStart.setZ(start.z());
                        qEnd.setX(end.x());
                        qEnd.setY(end.y());
                        qEnd.setZ(end.z());
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

std::vector<Vector3f>* expandContour(std::vector<Vector3f>* innerContour, Matrix3f m, float dist){
    using namespace ClipperLib;
    float precision = 1e8;

    Path subj;
    Paths solution;

    for (size_t i = 0; i < innerContour->size(); ++i) {
        Vector3f point = m * (*innerContour)[i];
        subj << IntPoint(point.x()*precision, point.z()*precision);
    }

    ClipperOffset co;
    co.AddPath(subj, jtSquare, etClosedPolygon);
    co.Execute(solution, dist * precision);

    std::vector<Vector3f> *expanded = new std::vector<Vector3f>();
    for (size_t i = 0; i < solution[0].size(); ++i) {

        Vector2f point = Vector2f(solution[0][i].X / precision, solution[0][i].Y / precision);

        float closestDist = std::numeric_limits<float>::infinity();
        float closestY;
        float dist;
        for (size_t j = 0; j < innerContour->size(); ++j) {
            Vector3f cpoint = m * (*innerContour)[j];
            dist = (cpoint.xz() - point).abs();

            if(dist < closestDist){
                closestDist = dist;
                closestY = cpoint.y();
            }
        }

        expanded->push_back(m.inverse() * Vector3f(solution[0][i].X / precision, closestY, solution[0][i].Y/ precision));
    }

    return expanded;
}

void mark_domains(CDT& ct, CDT::Face_handle start, int index, std::list<CDT::Edge>& border) {
    if (start->info().nesting_level != -1) {
        return;
    }
    std::list<CDT::Face_handle> queue;
    queue.push_back(start);
    while (!queue.empty()) {
        CDT::Face_handle fh = queue.front();
        queue.pop_front();
        if(fh->info().nesting_level == -1){
            fh->info().nesting_level = index;
            for(int i = 0; i < 3; i++){
                CDT::Edge e(fh,i);
                CDT::Face_handle n = fh->neighbor(i);
                if (n->info().nesting_level == -1) {
                    if (ct.is_constrained(e)) {
                        border.push_back(e);
                    } else {
                        queue.push_back(n);
                    }
                }
            }
        }
    }
}

GLuint addNewVertex(std::vector<Vector3f> *outerContour, std::vector<Vector3f> *innerContour,
    std::vector<Vector3f>& allPoints, Point& p, Matrix3f m){
    float minDist = 10e10;
    GLuint closestIndex = 0;
    for (GLuint i = 0; i < innerContour->size() + outerContour->size(); i++) {
        Vector3f v;
        if(i >= innerContour->size()){
            v = (*outerContour)[i - innerContour->size()];
        } else {
            v = (*innerContour)[i];
        }

        float dist = ((m * v).xz() - Vector2f(p.x(), p.y())).abs();
        if (dist < minDist) {
            minDist = dist;
            closestIndex = i;
        }
    }
    Vector3f wau;
    if(closestIndex >= innerContour->size()){
        wau = (*outerContour)[closestIndex - innerContour->size()];
    } else {
        wau = (*innerContour)[closestIndex];
    }
    Vector3f newPt = m.inverse() * Vector3f(p.x(), (m * wau).y(), p.y());
    allPoints.push_back(newPt);
    return allPoints.size() - 1;
}


void Mesh::fillInContour(std::vector<Vector3f> *outerContour, std::vector<Vector3f> *innerContour, Matrix3f m,
    std::vector<Vector3f> &allPoints, std::vector<GLuint> &allFaces) {
    std::vector<Point> inner2dPoints;
    std::vector<Point> outer2dPoints;

    for (GLuint i = 0; i < innerContour->size(); i++) {
        Vector3f v = (*innerContour)[i];
        Vector2f v2d = (m * v).xz();
        inner2dPoints.push_back(Point(v2d[0], v2d[1]));
        allPoints.push_back(v);
    }

    for (GLuint i = 0; i < outerContour->size(); i++) {
        Vector3f v = (*outerContour)[i];
        Vector2f v2d = (m * v).xz();
        outer2dPoints.push_back(Point(v2d[0], v2d[1]));
        // optimize this
        allPoints.push_back(v);
    }

    //Insert the polygons into a constrained triangulation
    CDT cdt;
    cdt.insert_constraint(inner2dPoints.begin(), inner2dPoints.end(), true);
    cdt.insert_constraint(outer2dPoints.begin(), outer2dPoints.end(), true);

    for(CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
        it->info().nesting_level = -1;
    }
    std::list<CDT::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while(!border.empty()) {
        CDT::Edge e = border.front();
        border.pop_front();
        CDT::Face_handle n = e.first->neighbor(e.second);
        if(n->info().nesting_level == -1){
            mark_domains(cdt, n, e.first->info().nesting_level+1, border);
        }
    }

    for (CDT::Finite_faces_iterator ft = cdt.finite_faces_begin();
        ft != cdt.finite_faces_end(); ++ft) {
        if (ft->info().nesting_level >= 1) {
            CDT::Face f = *ft;
            Point& p0 = ft->vertex(0)->point();
            Point& p1 = ft->vertex(1)->point();
            Point& p2 = ft->vertex(2)->point();

            bool found0 = false;
            bool found1 = false;
            bool found2 = false;
            int i0 = 0;
            int i1 = 0;
            int i2 = 0;

            for (GLuint i = 0; i < inner2dPoints.size(); i++) {
                if (inner2dPoints[i] == p0) {
                    found0 = true;
                    i0 = i;
                }
                if (inner2dPoints[i] == p1) {
                    found1 = true;
                    i1 = i;
                }
                if (inner2dPoints[i] == p2) {
                    found2 = true;
                    i2 = i;
                }
            }
            int offset = inner2dPoints.size();
            for (GLuint i = 0; i < outer2dPoints.size(); i++) {
                if (outer2dPoints[i] == p0) {
                    found0 = true;
                    i0 = i+offset;
                }
                if (outer2dPoints[i] == p1) {
                    found1 = true;
                    i1 = i+offset;
                }
                if (outer2dPoints[i] == p2) {
                    found2 = true;
                    i2 = i+offset;
                }
            }

            if (!found0) {
                i0 = addNewVertex(outerContour, innerContour, allPoints, p0, m);
            }
            if (!found1) {
                i1 = addNewVertex(outerContour, innerContour, allPoints, p1, m);
            }
            if (!found2) {
                i2 = addNewVertex(outerContour, innerContour, allPoints, p2, m);
            }

            // need to push in opposite direction
            allFaces.push_back(i0);
            allFaces.push_back(i2);
            allFaces.push_back(i1);
        }
    }
}

void populateEigen(std::vector<Vector3f>& verts, std::vector<GLuint>& faces,
    Eigen::MatrixXd& vm, Eigen::MatrixXi& fm) {
    vm = Eigen::MatrixXd(verts.size(), 3);
    for (size_t i = 0; i < verts.size(); i++) {
        vm(i, 0) = verts[i][0];
        vm(i, 1) = verts[i][1];
        vm(i, 2) = verts[i][2];
    }

    fm = Eigen::MatrixXi(faces.size() / 3, 3);
    for (size_t i = 0; i < faces.size()/3; i++) {
        fm(i, 0) = faces[3*i];
        fm(i, 1) = faces[3*i+1];
        fm(i, 2) = faces[3*i+2];
    }
}

Mesh* Mesh::generateMold(QVector3D n, float meshScale) {
    Vector3f norm = Vector3f(n.x(), n.y(), n.z());
    Vector3f a, b;
    if (n.x() <= n.y() && n.x() <= n.z()) {
        a = Vector3f::cross(Vector3f(1,0,0), norm).normalized();
    } else if (n.y() <= n.x() && n.x() <= n.z()) {
        a = Vector3f::cross(Vector3f(0,1,0), norm).normalized();
    } else {
        a = -Vector3f::cross(Vector3f(0,0,1), norm).normalized();
    }
    b = Vector3f::cross(a, norm);
    Matrix3f minv = Matrix3f();
    minv.setCol(0, a);
    minv.setCol(1, norm);
    minv.setCol(2, b);
    Matrix3f m = minv.inverse();

    std::vector<Vector3f>* edges = getEdges(norm);
    std::vector<Vector3f>* innerContour = sortIntoContour(edges, m);
    std::vector<Vector3f>* middleContour = expandContour(innerContour, m, 0.5/meshScale);
    std::vector<Vector3f>* outerContour = expandContour(innerContour, m, 1/meshScale);

    std::vector<Vector3f> cutSurfaceVerts;
    std::vector<GLuint> cutSurfaceFaces;

    std::vector<Vector3f> innerSurfaceVerts;
    std::vector<GLuint> innerSurfaceFaces;

    fillInContour(outerContour, innerContour, m, cutSurfaceVerts, cutSurfaceFaces);
    fillInContour(middleContour, new std::vector<Vector3f>(), m, innerSurfaceVerts, innerSurfaceFaces);

    std::pair<std::vector<Vector3f>, std::vector<GLuint>> topBlock
        = generateBlock(cutSurfaceVerts, cutSurfaceFaces,
        innerContour->size(), outerContour->size(), m, true);

    std::pair<std::vector<Vector3f>, std::vector<GLuint>> topCutout
        = generateBlock(innerSurfaceVerts, innerSurfaceFaces,
        0, middleContour->size(), m, true);

    std::vector<Vector3f> topBlockVerts = topBlock.first;
    std::vector<GLuint> topBlockFaces = topBlock.second;

    std::vector<Vector3f> topCutoutVerts = topCutout.first;
    std::vector<GLuint> topCutoutFaces = topCutout.second;

    // std::vector<Vector3f> botBlockVerts;
    // std::vector<GLuint> botBlockFaces;

    // generateOuterBlock(cutSurfaceVerts, cutSurfaceFaces, botBlockVerts, botBlockFaces, false,
    //     innerContour->size(), outerContour->size(), m);

    // std::ofstream f;
    // f.open ("test.stl");
    // qDebug() << " opening ";
    // if (f.is_open()) {
    //     qDebug() << " writing ";

    //     f << "solid test" << std::endl;

    //     for (GLuint i = 0; i < topBlockFaces.size(); i+= 3) {
    //         f << "facet normal 0 0 0" << std::endl;
    //         f << "  outer loop" << std::endl;
    //         for (GLuint j = 0; j < 3; j++) {
    //             GLuint idx = topBlockFaces[i + j];
    //             Vector3f v = topBlockVerts[idx];
    //             f << "    vertex " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    //         }
    //         f << "  endloop" << std::endl;
    //         f << "endfacet" << std::endl;
    //     }
    // }

    // f.close();

    Eigen::MatrixXd modelVM, topBlockVM, topCutoutVM, topMinus1VM, topMinus2VM, topFinalVM;
    Eigen::MatrixXi modelFM, topBlockFM, topCutoutFM, topMinus1FM, topMinus2FM, topFinalFM;

    modelVM = Eigen::MatrixXd(vertices.size() / 3, 3);
    for (size_t i = 0; i < vertices.size()/3; i++) {
        modelVM(i, 0) = vertices[3*i];
        modelVM(i, 1) = vertices[3*i+1];
        modelVM(i, 2) = vertices[3*i+2];
    }

    modelFM = Eigen::MatrixXi(indices.size() / 3, 3);
    for (size_t i = 0; i < indices.size()/3; i++) {
        modelFM(i, 0) = indices[3*i];
        modelFM(i, 1) = indices[3*i+2];
        modelFM(i, 2) = indices[3*i+1];
    }

    populateEigen(topBlockVerts, topBlockFaces, topBlockVM, topBlockFM);
    populateEigen(topCutoutVerts, topCutoutFaces, topCutoutVM, topCutoutFM);

    qDebug() << "first one";
    igl::copyleft::boolean::mesh_boolean(topBlockVM, topBlockFM, modelVM, modelFM,
        igl::copyleft::boolean::MESH_BOOLEAN_TYPE_MINUS, topMinus1VM, topMinus1FM);

    qDebug() << "second one";
    igl::copyleft::boolean::mesh_boolean(topCutoutVM, topCutoutFM, topMinus1VM, topMinus1FM,
        igl::copyleft::boolean::MESH_BOOLEAN_TYPE_INTERSECT, topMinus2VM, topMinus2FM);

    for (size_t i = 0; i < topMinus2VM.size()/3; i++) {
        Vector3f v = m * Vector3f(topMinus2VM(i, 0), topMinus2VM(i, 1), topMinus2VM(i, 2));
        Vector3f w = minv * Vector3f(v[0], v[1]-0.05/meshScale, v[2]);
        topMinus2VM(i,0) = w[0];
        topMinus2VM(i,1) = w[1];
        topMinus2VM(i,2) = w[2];
    }

    igl::copyleft::boolean::mesh_boolean(topMinus1VM, topMinus1FM, topMinus2VM, topMinus2FM,
        igl::copyleft::boolean::MESH_BOOLEAN_TYPE_MINUS, topFinalVM, topFinalFM);

    std::vector<GLfloat> mold_verts = std::vector<GLfloat>();
    std::vector<GLuint> mold_faces = std::vector<GLuint>();
    for (long i = 0; i < topFinalVM.rows(); ++i)
    {
        mold_verts.push_back(topFinalVM(i,0) * meshScale);
        mold_verts.push_back(topFinalVM(i,1) * meshScale);
        mold_verts.push_back(topFinalVM(i,2) * meshScale);
    }

    for (long i = 0; i < topFinalFM.rows(); ++i)
    {
        mold_faces.push_back(topFinalFM(i,0));
        mold_faces.push_back(topFinalFM(i,1));
        mold_faces.push_back(topFinalFM(i,2));
    }

    // std::ofstream f;
    // f.open ("test.stl");
    // qDebug() << " opening ";
    // if (f.is_open()) {
    //     qDebug() << " writing ";

    //     f << "solid test" << std::endl;

    //     for (GLuint i = 0; i < mold_faces.size(); i+= 3) {
    //         f << "facet normal 0 0 0" << std::endl;
    //         f << "  outer loop" << std::endl;
    //         for (GLuint j = 0; j < 3; j++) {
    //             GLuint idx = mold_faces[i + j];
    //             GLfloat x = mold_verts[3*idx];
    //             GLfloat y = mold_verts[3*idx+1];
    //             GLfloat z = mold_verts[3*idx+2];
    //             f << "    vertex " << x << " " << y << " " << z << std::endl;
    //         }
    //         f << "  endloop" << std::endl;
    //         f << "endfacet" << std::endl;
    //     }
    // }

    // f.close();

    return new Mesh(mold_verts, mold_faces);
}

std::pair<std::vector<Vector3f>, std::vector<GLuint>> Mesh::generateBlock(
    std::vector<Vector3f> &cutSurfaceVerts, std::vector<GLuint> &cutSurfaceFaces,
    size_t innerContourSize, size_t outerContourSize, Matrix3f m, bool isUpper) {

    // to orient vertices correctly
    int i1 = isUpper ? 1 : 2;
    int i2 = isUpper ? 2 : 1;

    size_t cutSurfaceVerts_size = cutSurfaceVerts.size();
    size_t cutSurfaceFaces_size = cutSurfaceFaces.size();

    std::vector<Vector3f> blockVerts = std::vector<Vector3f>();
    std::vector<GLuint> blockFaces = std::vector<GLuint>();

    for(size_t i = 0; i < cutSurfaceVerts_size; i++){
        blockVerts.push_back(cutSurfaceVerts[i]);
    }

    for(size_t i = 0; i < cutSurfaceFaces_size; i+=3){
        blockFaces.push_back(cutSurfaceFaces[i]);
        blockFaces.push_back(cutSurfaceFaces[i+i1]);
        blockFaces.push_back(cutSurfaceFaces[i+i2]);
    }

    Matrix3f minv = m.inverse();

    float zOffset = 0.2;

    // fill in flat bottom
    for(size_t i = 0; i < cutSurfaceVerts_size; i++){
        Vector2f v2d = (m * cutSurfaceVerts[i]).xz();
        Vector3f flatCoords = minv * Vector3f(v2d[0],
            (isUpper ? bbox.zmin - zOffset : bbox.zmax + zOffset) ,
            v2d[1]);
        blockVerts.push_back(flatCoords);
    }

    for(size_t i = 0; i < cutSurfaceFaces_size; i+=3){
        blockFaces.push_back(cutSurfaceFaces[i]+cutSurfaceVerts_size);
        blockFaces.push_back(cutSurfaceFaces[i+i2]+cutSurfaceVerts_size);
        blockFaces.push_back(cutSurfaceFaces[i+i1]+cutSurfaceVerts_size);
    }

    // Extrusion code
    for (GLuint i = innerContourSize; i < outerContourSize + innerContourSize; i++) {
        if (i - innerContourSize < outerContourSize - 1) {
            if (isUpper) {
                blockFaces.push_back(i);
                blockFaces.push_back(i+1);
                blockFaces.push_back(i+cutSurfaceVerts_size);
                blockFaces.push_back(i+1);
                blockFaces.push_back(i+1+cutSurfaceVerts_size);
                blockFaces.push_back(i+cutSurfaceVerts_size);
            } else {
                blockFaces.push_back(i);
                blockFaces.push_back(i+cutSurfaceVerts_size);
                blockFaces.push_back(i+1);
                blockFaces.push_back(i+1);
                blockFaces.push_back(i+cutSurfaceVerts_size);
                blockFaces.push_back(i+1+cutSurfaceVerts_size);
            }
        } else {
            if (isUpper) {
                blockFaces.push_back(i);
                blockFaces.push_back(innerContourSize);
                blockFaces.push_back(i+cutSurfaceVerts_size);
                blockFaces.push_back(innerContourSize);
                blockFaces.push_back(innerContourSize+cutSurfaceVerts_size);
                blockFaces.push_back(i+cutSurfaceVerts_size);
            } else {
                blockFaces.push_back(i);
                blockFaces.push_back(i+cutSurfaceVerts_size);
                blockFaces.push_back(innerContourSize);
                blockFaces.push_back(innerContourSize);
                blockFaces.push_back(i+cutSurfaceVerts_size);
                blockFaces.push_back(innerContourSize+cutSurfaceVerts_size);
            }
        }
    }

    return std::make_pair(blockVerts, blockFaces);
}

float get_line_intersection(Vector2f p0, Vector2f p1, Vector2f p2, Vector2f p3) {

    Vector2f s1 = p1 - p0;
    Vector2f s2 = p3 - p2;
    float eps = 10e-8;

    float s, t;
    s = (-s1[1] * (p0[0] - p2[0]) + s1[0] * (p0[1] - p2[1])) / (-s2[0] * s1[1] + s1[0] * s2[1]);
    t = ( s2[0] * (p0[1] - p2[1]) - s2[1] * (p0[0] - p2[0])) / (-s2[0] * s1[1] + s1[0] * s2[1]);

    if (s > eps && (1 - s) > eps && t > eps && (1 - t) > eps) {
        // Collision detected
        return t;
    }

    return -1; // No collision
}

// test if e1-e2 are an edge on triangle t1-t2-t3
bool edgeInTriangle(Vector3f t1, Vector3f t2, Vector3f t3, Vector3f e1, Vector3f e2) {
    return (((e1 == t1) && ((e2 == t2) || (e2 == t3))) ||
            ((e1 == t2) && ((e2 == t1) || (e2 == t3))) ||
            ((e1 == t3) && ((e2 == t1) || (e2 == t2))));
};

std::vector<Vector3f>* Mesh::getEdges(Vector3f n) {
    std::vector<bool> isCulled = std::vector<bool>(indices.size()/3);
    std::fill(isCulled.begin(), isCulled.end(), false);

    for (GLuint i = 0; i < indices.size(); i+=3) {
        GLuint xIdx = 3 * indices[i];
        GLuint yIdx = 3 * indices[i+1];
        GLuint zIdx = 3 * indices[i+2];

        Vector3f x = Vector3f(vertices[xIdx], vertices[xIdx+1], vertices[xIdx+2]);
        Vector3f y = Vector3f(vertices[yIdx], vertices[yIdx+1], vertices[yIdx+2]);
        Vector3f z = Vector3f(vertices[zIdx], vertices[zIdx+1], vertices[zIdx+2]);

        if (Vector3f::dot(Vector3f::cross(y-x, z-x).normalized(), n) <= 0) {
            isCulled[i/3] = true;
        }
    }

    std::vector<Vector3f> *edges = new std::vector<Vector3f>();

    for (GLuint i = 0; i < indices.size(); i+=3) {
        if (isCulled[i/3]) {
            for (GLuint j = 0; j < indices.size(); j+=3) {
                if (!isCulled[j/3]) {
                    GLuint xiIdx = 3 * indices[i];
                    GLuint yiIdx = 3 * indices[i+1];
                    GLuint ziIdx = 3 * indices[i+2];

                    Vector3f xi = Vector3f(vertices[xiIdx], vertices[xiIdx+1], vertices[xiIdx+2]);
                    Vector3f yi = Vector3f(vertices[yiIdx], vertices[yiIdx+1], vertices[yiIdx+2]);
                    Vector3f zi = Vector3f(vertices[ziIdx], vertices[ziIdx+1], vertices[ziIdx+2]);

                    GLuint xjIdx = 3 * indices[j];
                    GLuint yjIdx = 3 * indices[j+1];
                    GLuint zjIdx = 3 * indices[j+2];

                    Vector3f xj = Vector3f(vertices[xjIdx], vertices[xjIdx+1], vertices[xjIdx+2]);
                    Vector3f yj = Vector3f(vertices[yjIdx], vertices[yjIdx+1], vertices[yjIdx+2]);
                    Vector3f zj = Vector3f(vertices[zjIdx], vertices[zjIdx+1], vertices[zjIdx+2]);

                    if (edgeInTriangle(xj,yj,zj,xi,yi)) {
                        edges->push_back(xi);
                        edges->push_back(yi);
                    }
                    if (edgeInTriangle(xj,yj,zj,yi,zi)) {
                        edges->push_back(yi);
                        edges->push_back(zi);
                    }
                    if (edgeInTriangle(xj,yj,zj,xi,zi)) {
                        edges->push_back(xi);
                        edges->push_back(zi);
                    }
                }
            }
        }
    }

    return edges;
}

std::vector<Vector3f>* Mesh::sortIntoContour(std::vector<Vector3f>* edges, Matrix3f m) {
    // find point with largest x-coord when projected

    float runningMax = -10e10;
    GLuint rightmostIdx = 0;
    for (GLuint i = 0; i < edges->size(); i++) {
        Vector3f e = m * (*edges)[i];
        if (e.x() > runningMax) {
            runningMax = e.x();
            rightmostIdx = i;
        }
    }

    Vector3f rightmostPt = (*edges)[rightmostIdx];
    Vector3f currentPt = rightmostPt;
    GLuint currentIdx = rightmostIdx;

    std::vector<GLuint> innerContourIdxs = std::vector<GLuint>();
    std::vector<GLuint> uniqueIdxs = std::vector<GLuint>();
    do {
        GLuint neighborIdx = currentIdx + ((currentIdx % 2 == 0) ? 1 : -1);
        Vector3f neighbor = (*edges)[neighborIdx];

        float runningMin = 10e10;
        GLuint closestToNeighbor = -1;

        // find other point closest to that edge
        for (GLuint i = 0; i < edges->size(); i++) {
            // search only among points not already in the innerContour
            if (i != neighborIdx &&
                std::find(innerContourIdxs.begin(), innerContourIdxs.end(), i) == std::end(innerContourIdxs)) {
                float distToNeighbor = (neighbor - (*edges)[i]).abs();
                if (distToNeighbor < runningMin) {
                    closestToNeighbor = i;
                    runningMin = distToNeighbor;
                }
            }
        }

        if (runningMin > 1e10) {
            qDebug() << "ERROR";
            // all points have already been added to the innerContour
            break;
        }

        // add both to set of points to not take again
        innerContourIdxs.push_back(neighborIdx);
        innerContourIdxs.push_back(closestToNeighbor);

        // but if both points are approx the same, just add one of them to uniqueIdxs
        uniqueIdxs.push_back(neighborIdx);
        if (runningMin > 0.01) {
            uniqueIdxs.push_back(closestToNeighbor);
        }

        currentIdx = closestToNeighbor;
        currentPt = (*edges)[closestToNeighbor];

    } while (currentPt != rightmostPt);

    std::vector<Vector3f> *innerContour = new std::vector<Vector3f>();
    for (GLuint i = 0; i < uniqueIdxs.size(); i++) {
        innerContour->push_back((*edges)[uniqueIdxs[i]]);
    }

    return innerContour;
}
