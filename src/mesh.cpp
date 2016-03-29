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

Mesh* Mesh::getExtrudedOutline() {
    std::vector<GLuint>* contour = calculate2DContour();
    std::vector<GLfloat> new_verts;
    std::vector<GLuint> new_faces;

    for (GLuint i = 0; i < contour->size()-1; i++) {
        float x = vertices[(*contour)[i]];
        // float y = vertices[(*contour)[i+1]];
        float z = vertices[(*contour)[i]+2];

        qDebug() << (*contour)[i] << " " << x << " " << z;

        new_verts.push_back(x);
        new_verts.push_back(-1);
        new_verts.push_back(z);
        new_verts.push_back(x);
        new_verts.push_back(1);
        new_verts.push_back(z);

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

    return new Mesh(new_verts, new_faces);
}

std::vector<GLuint>* Mesh::calculate2DContour() {

    // find point with largest x coordinate
    GLuint rightMostIdx = 0;
    GLuint currentIdx, nextIdx;
    float runningMax = vertices[indices[0]];
    for (GLuint i = 0; i < vertices.size(); i += 3) {
        currentIdx = i; // x index of the ith vertex
        if (vertices[currentIdx] > runningMax) {
            runningMax = vertices[currentIdx];
            rightMostIdx = currentIdx;
        }
    };

    qDebug() << "rightMostIdx Is " << rightMostIdx;

    std::vector<GLuint> neighbors;
    std::vector<GLuint>* contour = new std::vector<GLuint>(0);
    float mostRecentAngle = 0;
    contour->push_back(rightMostIdx);
    currentIdx = rightMostIdx;
    QVector3D rightMostCoords = QVector3D(vertices[rightMostIdx], vertices[rightMostIdx+1], vertices[rightMostIdx+2]);

    QVector3D x,y,z, currentCoords;
    int ctr = 0;
    do {
        currentCoords = QVector3D(vertices[currentIdx], vertices[currentIdx+1], vertices[currentIdx+2]);
        qDebug() << currentIdx << " " << mostRecentAngle << " " << vertices[currentIdx+1];
        // get all neighbors to currentIdx
        // loop through list of faces
        for (GLuint i = 0; i < indices.size(); i+= 3) {
            GLuint xIdx = 3 * indices[i];
            GLuint yIdx = 3 * indices[i+1];
            GLuint zIdx = 3 * indices[i+2];
            // qDebug() << xIdx << " " << yIdx << " " << zIdx;

            x = QVector3D(vertices[xIdx], vertices[xIdx+1], vertices[xIdx+2]);
            y = QVector3D(vertices[yIdx], vertices[yIdx+1], vertices[yIdx+2]);
            z = QVector3D(vertices[zIdx], vertices[zIdx+1], vertices[zIdx+2]);
            float eps = 1e-16;
            if ((currentCoords - x).length() < eps) {
                neighbors.push_back(yIdx);
                neighbors.push_back(zIdx);
            } else if ((currentCoords - y).lengthSquared() < eps) {
                neighbors.push_back(xIdx);
                neighbors.push_back(zIdx);
            } else if ((currentCoords - z).lengthSquared() < eps) {
                neighbors.push_back(xIdx);
                neighbors.push_back(yIdx);
            }
        }

        // std::cout << "Neighbors: ";
        // for (size_t i = 0; i < neighbors.size(); i++) {
        //     std::cout << neighbors[i] << " ";
        // }
        // std::cout << std::endl;

        // find neighbor with smallest angle
        QVector2D currentPt = QVector2D(vertices[currentIdx], vertices[currentIdx+2]);
        qDebug() << "In 2D: " << currentIdx << " " << currentPt[0] << " " << currentPt[1] << atan2(currentPt[1], currentPt[0]);
        float runningMin = 2 * M_PI;
        float angleDelta;
        if (neighbors.size() == 0) {
            qDebug() << "Error: no neighbors found!";
        }

        nextIdx = neighbors[0];

        float nextMostRecentAngle;

        for (GLuint i = 0; i < neighbors.size(); i++) {
            int neighborIdx = neighbors[i];
            QVector2D neighborPt = QVector2D(vertices[neighborIdx], vertices[neighborIdx+2]);
            QVector2D displacement = neighborPt - currentPt;
            float absoluteAngle = atan2(displacement[1], displacement[0]);

            angleDelta = mostRecentAngle - absoluteAngle;
            if (angleDelta < 0 ) angleDelta += 2 * M_PI;

            // if (mostRecentAngle > absoluteAngle) {
            //     angleDelta = mostRecentAngle - absoluteAngle;
            // } else {
            //     angleDelta = mostRecentAngle - absoluteAngle + 2 * M_PI;
            // }

            // qDebug() << "Neighbor " << neighborIdx << " " << mostRecentAngle << " " << absoluteAngle << " " << angleDelta << " " << runningMin;

            if (angleDelta <= runningMin && angleDelta > 0) {
                runningMin = angleDelta;
                nextMostRecentAngle = absoluteAngle;
                nextIdx = neighborIdx;
            }
        }

        // add this neighbor to the contour
        contour->push_back(nextIdx);
        mostRecentAngle = nextMostRecentAngle;
        currentIdx = nextIdx;
        currentCoords = QVector3D(vertices[currentIdx], vertices[currentIdx+1], vertices[currentIdx+2]);
        neighbors.clear();
        ctr ++;
        // qDebug() << currentIdx << " " << mostRecentAngle;


    } while ((currentCoords - rightMostCoords).length() > 1e-16);
        //&& !(std::find(contour->begin(), contour->end(), nextIdx) != contour->end() - 1));

    return contour;
}