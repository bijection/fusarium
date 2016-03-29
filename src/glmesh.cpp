#include <iostream>

#include "glmesh.h"
#include "mesh.h"

GLMesh::GLMesh(const Mesh* const mesh)
    : mesh(mesh), vertices(QGLBuffer::VertexBuffer), indices(QGLBuffer::IndexBuffer)
{
    initializeGLFunctions();

    vertices.create();
    indices.create();

    vertices.setUsagePattern(QGLBuffer::StaticDraw);
    indices.setUsagePattern(QGLBuffer::StaticDraw);

    vertices.bind();
    vertices.allocate(mesh->vertices.data(),
                      mesh->vertices.size() * sizeof(float));
    vertices.release();

    indices.bind();
    indices.allocate(mesh->indices.data(),
                     mesh->indices.size() * sizeof(uint32_t));
    indices.release();
}

void GLMesh::draw(GLuint vp)
{
    vertices.bind();
    indices.bind();

    glVertexAttribPointer(vp, 3, GL_FLOAT, false, 3*sizeof(float), NULL);
    glDrawElements(GL_TRIANGLES, indices.size() / sizeof(uint32_t),
                   GL_UNSIGNED_INT, NULL);

    vertices.release();
    indices.release();
}

bool GLMesh::checkForOverhangs(QVector3D norm) {
    hasOverhang = mesh->checkForOverhangs(norm, start, end);
    return hasOverhang;
}

void GLMesh::drawOverhangLine() {
    if (hasOverhang) {
        glColor3f(1.0, 1.0, 1.0);
        glLineWidth(5);
        glBegin(GL_LINES);
        glVertex3f(start.x(), start.y(), start.z());
        glVertex3f(end.x(), end.y(), end.z());
        glEnd();
    }
}

void GLMesh::drawBoundingBox()
{
    BoundingBox bbox = mesh->bbox;
    float x0 = bbox.xmin;
    float x1 = bbox.xmax;
    float y0 = bbox.ymin;
    float y1 = bbox.ymax;
    float z0 = bbox.zmin;
    float z1 = bbox.zmax;

    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(2);
    glBegin(GL_LINES);
    glVertex3f(x0, y0, z0);
    glVertex3f(x1, y0, z0);
    glVertex3f(x0, y0, z0);
    glVertex3f(x0, y1, z0);
    glVertex3f(x0, y0, z0);
    glVertex3f(x0, y0, z1);
    glVertex3f(x1, y1, z1);
    glVertex3f(x0, y1, z1);
    glVertex3f(x1, y1, z1);
    glVertex3f(x1, y0, z1);
    glVertex3f(x1, y1, z1);
    glVertex3f(x1, y1, z0);
    glVertex3f(x1, y0, z0);
    glVertex3f(x1, y1, z0);
    glVertex3f(x1, y0, z0);
    glVertex3f(x1, y0, z1);
    glVertex3f(x0, y1, z0);
    glVertex3f(x1, y1, z0);
    glVertex3f(x0, y1, z0);
    glVertex3f(x0, y1, z1);
    glVertex3f(x0, y0, z1);
    glVertex3f(x1, y0, z1);
    glVertex3f(x0, y0, z1);
    glVertex3f(x0, y1, z1);
    glEnd();
    glLineStipple(4, 0xAAAA);
    glEnable(GL_LINE_STIPPLE);
    glLineWidth(0.1);
    glBegin(GL_LINES);
    glVertex3f(x0, y0, (z0+z1)/2);
    glVertex3f(x0, y1, (z0+z1)/2);
    glVertex3f(x0, y0, (z0+z1)/2);
    glVertex3f(x1, y0, (z0+z1)/2);
    glVertex3f(x1, y1, (z0+z1)/2);
    glVertex3f(x0, y1, (z0+z1)/2);
    glVertex3f(x1, y1, (z0+z1)/2);
    glVertex3f(x1, y0, (z0+z1)/2);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
}