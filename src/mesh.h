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
    std::vector<Vector2f>* get2DGraph(Matrix3f m);
    std::vector<Vector2f>* calculate2DContour(std::vector<Vector2f>* edges2d);
    std::vector<GLfloat> vertices;
    std::vector<GLuint> indices;
    QMatrix4x4 transform;

    friend class GLMesh;
};

#endif // MESH_H
