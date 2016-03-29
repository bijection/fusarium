#ifndef MESH_H
#define MESH_H

#include <QString>
#include <QtOpenGL/QtOpenGL>
#include <QMatrix4x4>
#include <vector>
#include <../vecmath/vecmath.h>


struct BoundingBox {
    float xmin;
    float xmax;
    float ymin;
    float ymax;
    float zmin;
    float zmax;
};

class Mesh
{
public:
    Mesh(std::vector<GLfloat> vertices, std::vector<GLuint> indices);
    void setTransform(float rotateX, float rotateY, float rotateZ);
    float calculateProjectedArea(QVector3D norm);
    bool checkForOverhangs(QVector3D norm, QVector3D &start, QVector3D &end) const;
    Mesh* getExtrudedOutline(QVector3D norm);
    BoundingBox bbox;

private:
    std::vector<Vector3f>* getEdges(Vector3f n, Matrix3f m);
    std::vector<Vector3f>* sortIntoContour(std::vector<Vector3f>* edges, Matrix3f m);
    std::vector<GLfloat> vertices;
    std::vector<GLuint> indices;
    QMatrix4x4 transform;

    friend class GLMesh;
};

#endif // MESH_H
