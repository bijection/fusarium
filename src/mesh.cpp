#include <QFile>
#include <QDataStream>
#include <QVector3D>

#include <iostream>
#include <cmath>
#include "../vecmath/vecmath.h"
#include "mesh.h"

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

    std::vector<Vector2f>* edges2d = get2DGraph(m);
    std::vector<Vector2f>* contour = calculate2DContour(edges2d);
    std::vector<GLfloat> new_verts;
    std::vector<GLuint> new_faces;

    for (GLuint i = 0; i < contour->size()-1; i++) {
        Vector3f projTopCoords = minv * Vector3f((*contour)[i][0],  5, (*contour)[i][1]);
        Vector3f projBotCoords = minv * Vector3f((*contour)[i][0],  -5, (*contour)[i][1]);

        new_verts.push_back(projTopCoords.x());
        new_verts.push_back(projTopCoords.y());
        new_verts.push_back(projTopCoords.z());
        new_verts.push_back(projBotCoords.x());
        new_verts.push_back(projBotCoords.y());
        new_verts.push_back(projBotCoords.z());

        if (i != contour->size() - 2) {
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

    GLuint offset = new_faces.size();
    for (GLuint i = 0; i < edges2d->size()-1; i++) {
        Vector3f projTopCoords = minv * Vector3f((*edges2d)[i][0],  1, (*edges2d)[i][1]);
        Vector3f projBotCoords = minv * Vector3f((*edges2d)[i][0],  -1, (*edges2d)[i][1]);

        new_verts.push_back(projTopCoords.x());
        new_verts.push_back(projTopCoords.y());
        new_verts.push_back(projTopCoords.z());
        new_verts.push_back(projBotCoords.x());
        new_verts.push_back(projBotCoords.y());
        new_verts.push_back(projBotCoords.z());

        if (i%2 == 0) {
            new_faces.push_back(offset + 2*i);
            new_faces.push_back(offset + 2*i + 2);
            new_faces.push_back(offset + 2*i + 3);
            new_faces.push_back(offset + 2*i);
            new_faces.push_back(offset + 2*i + 3);
            new_faces.push_back(offset + 2*i + 1);
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

std::vector<Vector2f>* Mesh::get2DGraph(Matrix3f m) {
    std::vector<Vector2f> orig_edges2d = std::vector<Vector2f>(0);
    std::vector<Vector2f>* edges2d = new std::vector<Vector2f>(0);
    float eps = 10e-16;

    // project all edges down to xz plane
    for (GLuint i = 0; i < indices.size(); i += 3) {
        GLuint xIdx = 3 * indices[i];
        GLuint yIdx = 3 * indices[i+1];
        GLuint zIdx = 3 * indices[i+2];

        Vector2f x2d = (m * Vector3f(vertices[xIdx], vertices[xIdx+1], vertices[xIdx+2])).xz();
        Vector2f y2d = (m * Vector3f(vertices[yIdx], vertices[yIdx+1], vertices[yIdx+2])).xz();
        Vector2f z2d = (m * Vector3f(vertices[zIdx], vertices[zIdx+1], vertices[zIdx+2])).xz();

        orig_edges2d.push_back(x2d);
        orig_edges2d.push_back(y2d);
        orig_edges2d.push_back(x2d);
        orig_edges2d.push_back(z2d);
        orig_edges2d.push_back(y2d);
        orig_edges2d.push_back(z2d);
    }

    // account for intersections
    for (GLuint i = 0; i < orig_edges2d.size(); i += 2) {
        if (i % (orig_edges2d.size()/10) == 0) {
            qDebug() << i << " out of " << orig_edges2d.size();
        }
        std::vector<float> ts = std::vector<float>(0);
        std::vector<float> sort_ts = std::vector<float>(0);
        Vector2f p0 = orig_edges2d[i];
        Vector2f p1 = orig_edges2d[i+1];
        // if this edge is not vertical
        if ((p0 - p1).abs() > eps) {
            for (GLuint j = 0; j < orig_edges2d.size(); j += 2) {
                if (i != j) {
                    Vector2f p2 = orig_edges2d[j];
                    Vector2f p3 = orig_edges2d[j+1];
                    float t = get_line_intersection(p0, p1, p2, p3);
                    if (t > 0) {
                        ts.push_back(t);
                    }
                }
            }
            // deduplicate
            ts.insert(ts.begin(), 0);
            ts.push_back(1);

            std::sort(ts.begin(), ts.end());
            std::unique_copy(ts.begin(), ts.end(), std::back_inserter(sort_ts),
                [](double l, double r) { return std::abs(l - r) < 0.01; });

            // for (auto t: sort_ts) {
            //     std::cout << t << " ";
            // }
            // std::cout << std::endl;

            for (GLuint i = 0; i < sort_ts.size() - 1; i++) {
                edges2d->push_back(p0 + sort_ts[i] * (p1-p0));
                edges2d->push_back(p0 + sort_ts[i+1] * (p1-p0));
            }
        }
    }

    qDebug() << orig_edges2d.size() << " vs " << edges2d->size();
    return edges2d;
}

std::vector<Vector2f>* Mesh::calculate2DContour(std::vector<Vector2f>* edges2dPtr) {
    float eps = 10e-16;
    std::vector<Vector2f> edges2d = *edges2dPtr;
    // find point with largest x coordinate
    Vector2f rightmostPt = edges2d[0];

    float runningMax = -10e10;
    for (GLuint i = 0; i < edges2d.size(); i++) {
        Vector2f v = edges2d[i];
        if (v.x() > runningMax) {
            runningMax = v.x();
            rightmostPt = v;
        }
    };

    qDebug() << "rightmostPt Is :";
    rightmostPt.print();

    std::vector<Vector2f> neighbors;
    std::vector<Vector2f>* contour = new std::vector<Vector2f>(0);
    float mostRecentAngle = 0;
    contour->push_back(rightmostPt);

    Vector2f currentPt = rightmostPt;
    float cumulativeAngle = 0;

    int ctr = 0;
    do {
        // loop through list of edges
        for (GLuint i = 0; i < edges2d.size(); i+= 2) {
            Vector2f x = edges2d[i];
            Vector2f y = edges2d[i+1];
            float eps = 2;

            if ((x - currentPt).abs() < eps) {
                neighbors.push_back(y);
            } else if ((y - currentPt).abs() < eps) {
                neighbors.push_back(x);
            }

        }

        float runningMin = 2 * M_PI;
        float angleDelta;
        if (neighbors.size() == 0) {
            qDebug() << "Error: no neighbors found!";
        }

        Vector2f nextPt = neighbors[0];
        float nextMostRecentAngle;


        for (GLuint i = 0; i < neighbors.size(); i++) {
            Vector2f neighborPt = neighbors[i];

            Vector2f displacement = neighborPt - currentPt;
            float absoluteAngle = atan2(displacement[1], displacement[0]);

            angleDelta = mostRecentAngle - absoluteAngle;
            if (angleDelta < 0 ) angleDelta += 2 * M_PI;

            // qDebug() << "Neighbor " << " " << mostRecentAngle << " " << absoluteAngle << " " << angleDelta << " " << runningMin;

            if (angleDelta <= runningMin && angleDelta > 0) {
                runningMin = angleDelta;
                nextMostRecentAngle = absoluteAngle;
                nextPt = neighborPt;
            }
        }

        // check all other edges if when projected down to XZ they cross this line
        bool hasIntersection = false;
        float tIntersect = 1;
        for (GLuint i = 0; i < edges2d.size(); i+= 3) {
            Vector2f cross1 = edges2d[i];
            Vector2f cross2 = edges2d[i+1];

            if ((cross1 - currentPt).absSquared() > eps && (cross2 - currentPt).absSquared() > eps) {

                float xyt = get_line_intersection(currentPt, nextPt, cross1, cross2);

                if (xyt > 0) {
                    qDebug() << "INTERSECTION FOUND: " << xyt;
                    tIntersect = fmin(tIntersect, xyt);
                }
            }
        }

        if (tIntersect < 1) {
            nextPt = currentPt + tIntersect * (nextPt - currentPt);
        }


        // add this neighbor to the contour
        cumulativeAngle += angleDelta - M_PI;
        mostRecentAngle = nextMostRecentAngle;
        currentPt = nextPt;
        contour->push_back(currentPt);

        qDebug() << (angleDelta - M_PI) << cumulativeAngle;

        neighbors.clear();
        ctr ++;

    } while ((currentPt - rightmostPt).absSquared() > 1e-16 && ctr < 50);
        //&& !(std::find(contour->begin(), contour->end(), nextIdx) != contour->end() - 1));

    qDebug() << contour->size();

    return contour;
}