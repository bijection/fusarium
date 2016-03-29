#ifndef GLMESH_H
#define GLMESH_H

#include <QtOpenGL/QGLBuffer>
#include <QtOpenGL/QGLFunctions>
#include <QVector3D>

class Mesh;

class GLMesh : protected QGLFunctions
{
public:
    GLMesh(const Mesh* const mesh);
    void draw(GLuint vp);
    void drawBoundingBox();
    void drawOverhangLine();
	bool checkForOverhangs(QVector3D norm);
private:
	const Mesh* mesh;
    QGLBuffer vertices;
    QGLBuffer indices;
    QVector3D start, end;
    bool hasOverhang = false;
};

#endif // GLMESH_H
