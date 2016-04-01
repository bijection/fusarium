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
    Mesh* generateMold(QVector3D norm, float meshScale);
    BoundingBox bbox;

private:
    std::vector<Vector3f>* getEdges(Vector3f n);
    std::vector<Vector3f>* sortIntoContour(std::vector<Vector3f>* edges, Matrix3f m);
    std::vector<GLfloat> vertices;
    std::vector<GLuint> indices;
    void fillInContour(std::vector<Vector3f> *contour, std::vector<Vector3f> *bigContour, Matrix3f m,
        std::vector<Vector3f> &new_verts, std::vector<GLuint> &new_faces);
    std::pair<std::vector<Vector3f>, std::vector<GLuint>> generateBlock(
        std::vector<Vector3f> &cutSurfaceVerts, std::vector<GLuint> &cutSurfaceFaces,
        size_t innerContourSize, size_t outerContourSize, Matrix3f m, bool isUpper);

    QMatrix4x4 transform;

    friend class GLMesh;
};

#endif // MESH_H
