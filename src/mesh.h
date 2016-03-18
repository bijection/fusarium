#ifndef MESH_H
#define MESH_H

#include <QString>
#include <QtOpenGL/QtOpenGL>
#include <QMatrix4x4>
#include <vector>


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
    BoundingBox bbox;

private:
    std::vector<GLfloat> vertices;
    std::vector<GLuint> indices;
    QMatrix4x4 transform;

    friend class GLMesh;
};

#endif // MESH_H
