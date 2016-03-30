#include <QFile>
#include <QDataStream>
#include <QVector3D>

#include <iostream>
#include <cmath>
#include "../vecmath/vecmath.h"
#include "mesh.h"

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

std::vector<Vector3f> *Mesh::fillInContour(std::vector<Vector3f> *contour, Matrix3f m,
    std::vector<GLfloat> &new_verts, std::vector<GLuint> &new_faces) {
    std::vector<Point> points;
    std::vector<Vector3f> *newContour = new std::vector<Vector3f>();
    Matrix3f minv = m.inverse();

    float aveY = 0;
    for (GLuint i = 0; i < contour->size(); i++) {
        Vector3f v = (*contour)[i];
        aveY += v.y();
        Vector2f v2d = (m * v).xz();
        points.push_back(Point(v2d[0], v2d[1]));
        // qDebug() << i << v2d[0] << v2d[1];
        new_verts.push_back(v[0]);
        new_verts.push_back(v[1]);
        new_verts.push_back(v[2]);
    }
    aveY = aveY / contour->size();

    // aveY = 100;

    qDebug() << contour->size() << " " << new_verts.size();

    //Insert the polygons into a constrained triangulation
    CDT cdt;
    cdt.insert_constraint(points.begin(), points.end(), true);

    // for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
    //     vit != cdt.finite_vertices_end(); ++vit) {
    //     Point& p0 = vit->point();
    //     vit->info().index = -5;
    //     // qDebug() << p0.x() << p0.y();
    // }

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

            for (GLuint i = 0; i < points.size(); i++) {
                if (points[i] == p0) {
                    found0 = true;
                    i0 = i;
                }
                if (points[i] == p1) {
                    found1 = true;
                    i1 = i;
                }
                if (points[i] == p2) {
                    found2 = true;
                    i2 = i;
                }
            }

            if (!found0) {
                float minDist = 10e10;
                GLuint closestPt = 0;
                for (GLuint i = 0; i < contour->size(); i++) {
                    float dist = ((m * (*contour)[i]).xz() - Vector2f(p0.x(), p0.y())).abs();
                    if (dist < minDist) {
                        minDist = dist;
                        closestPt = i;
                    }
                }
                Vector3f newPt = minv * Vector3f(p0.x(), (m * (*contour)[closestPt]).y(), p0.y());
                new_verts.push_back(newPt[0]);
                new_verts.push_back(newPt[1]);
                new_verts.push_back(newPt[2]);
                i0 = (new_verts.size() - 3)/3;
                newPt.print();
            }
            if (!found1) {
                float minDist = 10e10;
                GLuint closestPt = 0;
                for (GLuint i = 0; i < contour->size(); i++) {
                    float dist = ((m * (*contour)[i]).xz() - Vector2f(p1.x(), p1.y())).abs();
                    if (dist < minDist) {
                        minDist = dist;
                        closestPt = i;
                    }
                }
                Vector3f newPt = minv * Vector3f(p1.x(), (m * (*contour)[closestPt]).y(), p1.y());
                new_verts.push_back(newPt[0]);
                new_verts.push_back(newPt[1]);
                new_verts.push_back(newPt[2]);
                i1 = (new_verts.size() - 3)/3;
                newPt.print();
            }
            if (!found2) {
                float minDist = 10e10;
                GLuint closestPt = 0;
                for (GLuint i = 0; i < contour->size(); i++) {
                    float dist = ((m * (*contour)[i]).xz() - Vector2f(p2.x(), p2.y())).abs();
                    if (dist < minDist) {
                        minDist = dist;
                        closestPt = i;
                    }
                }
                Vector3f newPt = minv * Vector3f(p2.x(), (m * (*contour)[closestPt]).y(), p2.y());
                new_verts.push_back(newPt[0]);
                new_verts.push_back(newPt[1]);
                new_verts.push_back(newPt[2]);
                i2 = (new_verts.size() - 3)/3;
                newPt.print();
            }

            qDebug() << i0 << i1 << i2;
            // if (found0 && found1 && found2) {
                new_faces.push_back(i0);
                new_faces.push_back(i1);
                new_faces.push_back(i2);
            // }
        }
    }


    // for (CDT::Finite_edges_iterator eit = cdt.finite_edges_begin();
    //     eit != cdt.finite_edges_end(); ++eit){

    //     CDT::Edge f = *eit;
    //     CDT::Segment s = cdt.segment(eit);
    //     const Point& p0 = s.point(0);
    //     const Point& p1 = s.point(1);

    //     Vector3f v0 = minv * Vector3f(p0.x(), 0, p0.y());
    //     Vector3f v1 = minv * Vector3f(p1.x(), 0, p1.y());
    //     newContour->push_back(v0);
    //     newContour->push_back(v1);

    // }
    return newContour;
}

Mesh* Mesh::getExtrudedOutline(QVector3D n) {
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
    std::vector<Vector3f>* contour = sortIntoContour(edges, m);

    // std::vector<Vector3f>* oldContour = sortIntoContour(edges, m);
    // std::vector<Vector3f>* contour = fillInContour(oldContour, m, new_verts, new_faces);

    std::vector<GLfloat> new_verts;
    std::vector<GLuint> new_faces;

    // for (GLuint i = 0; i < contour->size(); i++) {
    //     Vector2f point2d = (m * (*contour)[i]).xz();
    //     Vector3f projTopCoords = minv * Vector3f(point2d[0], 3, point2d[1]);
    //     Vector3f projBotCoords = minv * Vector3f(point2d[0], -3, point2d[1]);

    //     new_verts.push_back(projTopCoords.x());
    //     new_verts.push_back(projTopCoords.y());
    //     new_verts.push_back(projTopCoords.z());
    //     new_verts.push_back(projBotCoords.x());
    //     new_verts.push_back(projBotCoords.y());
    //     new_verts.push_back(projBotCoords.z());

    //     if (i %2 == 0) {
    //         new_faces.push_back(2*i);
    //         new_faces.push_back(2*i + 2);
    //         new_faces.push_back(2*i + 3);
    //         new_faces.push_back(2*i);
    //         new_faces.push_back(2*i + 3);
    //         new_faces.push_back(2*i + 1);
    //     }


    //     if (i < contour->size()-1) {
    //         new_faces.push_back(2*i);
    //         new_faces.push_back(2*i + 2);
    //         new_faces.push_back(2*i + 3);
    //         new_faces.push_back(2*i);
    //         new_faces.push_back(2*i + 3);
    //         new_faces.push_back(2*i + 1);
    //     } else {
    //         new_faces.push_back(2*i);
    //         new_faces.push_back(0);
    //         new_faces.push_back(1);
    //         new_faces.push_back(2*i);
    //         new_faces.push_back(1);
    //         new_faces.push_back(2*i + 1);
    //     }
    // }

    fillInContour(contour, m, new_verts, new_faces);
    qDebug() << " Here: " << new_verts.size();

    return new Mesh(new_verts, new_faces);
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

    std::vector<GLuint> contourIdxs = std::vector<GLuint>();
    std::vector<GLuint> uniqueIdxs = std::vector<GLuint>();
    do {
        GLuint neighborIdx = currentIdx + ((currentIdx % 2 == 0) ? 1 : -1);
        Vector3f neighbor = (*edges)[neighborIdx];

        float runningMin = 10e10;
        GLuint closestToNeighbor = -1;

        // find other point closest to that edge
        for (GLuint i = 0; i < edges->size(); i++) {
            // search only among points not already in the contour
            if (i != neighborIdx &&
                std::find(contourIdxs.begin(), contourIdxs.end(), i) == std::end(contourIdxs)) {
                float distToNeighbor = (neighbor - (*edges)[i]).abs();
                if (distToNeighbor < runningMin) {
                    closestToNeighbor = i;
                    runningMin = distToNeighbor;
                }
            }
        }

        if (runningMin > 1e10) {
            qDebug() << "ERROR";
            // all points have already been added to the contour
            break;
        }

        // add both to set of points to not take again
        contourIdxs.push_back(neighborIdx);
        contourIdxs.push_back(closestToNeighbor);

        // but if both points are approx the same, just add one of them to uniqueIdxs
        uniqueIdxs.push_back(neighborIdx);
        if (runningMin > 0.01) {
            uniqueIdxs.push_back(closestToNeighbor);
        }

        currentIdx = closestToNeighbor;
        currentPt = (*edges)[closestToNeighbor];

    } while (currentPt != rightmostPt);

    std::vector<Vector3f> *contour = new std::vector<Vector3f>();
    for (GLuint i = 0; i < uniqueIdxs.size(); i++) {
        contour->push_back((*edges)[uniqueIdxs[i]]);
    }

    // qDebug() << contourIdxs.size() << " vs " << uniqueIdxs.size() << " vs " << edges->size();

    return contour;
}
