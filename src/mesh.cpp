#include <QFile>
#include <QDataStream>
#include <QVector3D>

#include <iostream>
#include <cmath>
#include <limits>
#include "../vecmath/vecmath.h"
#include "mesh.h"
#include "clipper.hpp"


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

std::vector<Vector3f>* expandContour(std::vector<Vector3f>* contour, Matrix3f m){
    using namespace ClipperLib;

    float precision = 1e8;

    Path subj;
    Paths solution;
    
    qDebug() << "setting up path, contour size "<< contour->size();

    for (size_t i = 0; i < contour->size(); ++i) {
        Vector3f point = m * (*contour)[i];
        qDebug() << "creating point" << point.x() << " " << point.z();
        subj << IntPoint(point.x()*precision, point.z()*precision);
    }

    qDebug() << "expanding";

    ClipperOffset co;
    co.AddPath(subj, jtMiter, etClosedPolygon);
    co.Execute(solution, 1e8);

    qDebug() << "converting expanded path to 4d " << solution.size();


    // for (size_t i = 0; i < solution.size(); ++i) {
    // for (size_t j = 0; j < solution[i].size(); ++j) {
    //     qDebug() << solution[i][j].X << " " << solution[i][j].Y;
    // }}

    std::vector<Vector3f> *expanded = new std::vector<Vector3f>();
    for (size_t i = 0; i < solution[0].size(); ++i) {
        
        Vector2f point = Vector2f(solution[0][i].X / precision, solution[0][i].Y / precision);

        float closestDist = std::numeric_limits<float>::infinity();
        float closestY;
        float dist;
        for (size_t j = 0; j < contour->size(); ++j) {
            Vector3f cpoint = m * (*contour)[j];
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

    std::vector<Vector3f>* edges = getEdges(norm, m);
    std::vector<Vector3f>* contour = sortIntoContour(edges, m);
    qDebug() << "Preparing to expand";

    std::vector<Vector3f>* bigContour = expandContour(contour, m);

    contour = bigContour;

    std::vector<GLfloat> new_verts;
    std::vector<GLuint> new_faces;


    for (GLuint i = 0; i < contour->size(); i++) {
        Vector2f point2d = (m * (*contour)[i]).xz();
        Vector3f projTopCoords = minv * Vector3f(point2d[0], 3, point2d[1]);
        Vector3f projBotCoords = (*contour)[i];
        // Vector3f projBotCoords = minv * Vector3f(point2d[0], -3, point2d[1]);

        new_verts.push_back(projTopCoords.x());
        new_verts.push_back(projTopCoords.y());
        new_verts.push_back(projTopCoords.z());
        new_verts.push_back(projBotCoords.x());
        new_verts.push_back(projBotCoords.y());
        new_verts.push_back(projBotCoords.z());

        if (i < contour->size()-1) {
            new_faces.push_back(2*i);
            new_faces.push_back(2*i + 2);
            new_faces.push_back(2*i + 3);
            new_faces.push_back(2*i);
            new_faces.push_back(2*i + 3);
            new_faces.push_back(2*i + 1);
        } else {
            new_faces.push_back(2*i);
            new_faces.push_back(0);
            new_faces.push_back(1);
            new_faces.push_back(2*i);
            new_faces.push_back(1);
            new_faces.push_back(2*i + 1);

        }
    }


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

std::vector<Vector3f>* Mesh::getEdges(Vector3f n, Matrix3f m) {
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
