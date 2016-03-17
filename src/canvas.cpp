#include <QMouseEvent>
#include <QDebug>

#include <cmath>
#include <QtDebug>
#include <iostream>

#include "canvas.h"
#include "backdrop.h"
#include "glmesh.h"
#include "mesh.h"

Canvas::Canvas(const QGLFormat& format, QWidget *parent)
    : QGLWidget(format, parent), mesh(NULL),
      scale(1), zoom(1), tilt(90), yaw(0), status(" ")
{
    // Nothing to do here
}

Canvas::~Canvas()
{
    delete glmesh;
}

void Canvas::load_mesh(Mesh* m)
{
    mesh = m;
    glmesh = new GLMesh(m);

    meshRotateX = 0;
    meshRotateY = 0;
    meshRotateZ = 0;
    updateMeshBbox();
    BoundingBox bbox = mesh->bbox;

    center = QVector3D(bbox.xmin + bbox.xmax,
                       bbox.ymin + bbox.ymax,
                       bbox.zmin + bbox.zmax) / 2;
    scale = 2 / sqrt(
                pow(bbox.xmax - bbox.xmin, 2) +
                pow(bbox.ymax - bbox.ymin, 2) +
                pow(bbox.zmax - bbox.zmin, 2));

    // Reset other camera parameters
    zoom = 1;
    yaw = 0;
    tilt = 90;

    update();
}

void Canvas::set_status(const QString &s)
{
    status = s;
    update();
}

void Canvas::clear_status()
{
    status = "";
    update();
}

void Canvas::setMeshRotateX(int degrees)
{
    meshRotateX = degrees + 0.0;
    updateMeshBbox();
    update();
}

void Canvas::setMeshRotateY(int degrees)
{
    meshRotateY = degrees + 0.0;
    updateMeshBbox();
    update();
}

void Canvas::setMeshRotateZ(int degrees)
{
    meshRotateZ = degrees + 0.0;
    updateMeshBbox();
    update();
}

void Canvas::setMeshScale(int factor)
{
    if (factor > 0) {
        meshScale = factor/10.0f + 1.0f;
    } else if (factor < 0) {
        meshScale = 1 / (-factor/10.0f + 1.0f);
    } else {
        meshScale = 1;
    }
    updateMeshBbox();
    update();
}

void Canvas::updateMeshBbox() {
    mesh->setTransform(meshRotateX, meshRotateY, meshRotateZ);
    BoundingBox bbox = mesh->bbox;
    center = QVector3D(bbox.xmin + bbox.xmax,
                       bbox.ymin + bbox.ymax,
                       bbox.zmin + bbox.zmax) / 2;
    emit updatedBbox(
        (bbox.xmax - bbox.xmin) * meshScale,
        (bbox.ymax - bbox.ymin) * meshScale,
        (bbox.zmax - bbox.zmin) * meshScale);
}

void Canvas::initializeGL()
{
    initializeGLFunctions();

    mesh_shader.addShaderFromSourceFile(QGLShader::Vertex, ":/gl/mesh.vert");
    mesh_shader.addShaderFromSourceFile(QGLShader::Fragment, ":/gl/mesh.frag");
    mesh_shader.link();

    line_shader.addShaderFromSourceFile(QGLShader::Vertex, ":/gl/mesh.vert");
    line_shader.addShaderFromSourceFile(QGLShader::Fragment, ":/gl/lines.frag");
    line_shader.link();

    backdrop = new Backdrop();
}

void Canvas::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event);

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glHint(GL_CLIP_VOLUME_CLIPPING_HINT_EXT,GL_FASTEST);

    backdrop->draw();
    if (glmesh) draw_mesh();

    if (status.isNull())    return;

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.drawText(10, height() - 10, status);
}


void Canvas::draw_mesh()
{
    mesh_shader.bind();

    // Load the transform and view matrices into the shader
    glUniformMatrix4fv(
                mesh_shader.uniformLocation("transform_matrix"),
                1, GL_FALSE, transform_matrix().data());
    glUniformMatrix4fv(
                mesh_shader.uniformLocation("view_matrix"),
                1, GL_FALSE, view_matrix().data());

    // Compensate for z-flattening when zooming
    glUniform1f(mesh_shader.uniformLocation("zoom"), 1/(zoom * meshScale));

    // Find and enable the attribute location for vertex position
    const GLuint vp = mesh_shader.attributeLocation("vertex_position");
    glEnableVertexAttribArray(vp);

    // Then draw the mesh with that vertex position
    glmesh->draw(vp);
    mesh_shader.release();

    line_shader.bind();
    glUniformMatrix4fv(
                line_shader.uniformLocation("transform_matrix"),
                1, GL_FALSE, bbox_transform_matrix().data());
    glUniformMatrix4fv(
                line_shader.uniformLocation("view_matrix"),
                1, GL_FALSE, view_matrix().data());
    glmesh->drawBoundingBox();

    // Clean up state machine
    glDisableVertexAttribArray(vp);
    line_shader.release();
}

QMatrix4x4 Canvas::transform_matrix() const
{
    QMatrix4x4 m;
    m.rotate(tilt, QVector3D(1, 0, 0));
    m.rotate(yaw,  QVector3D(0, 0, 1));
    m.scale(scale);
    m.translate(-center);
    m.rotate(meshRotateX, QVector3D(1, 0, 0));
    m.rotate(meshRotateY,  QVector3D(0, 1, 0));
    m.rotate(meshRotateZ,  QVector3D(0, 0, 1));
    return m;
}

QMatrix4x4 Canvas::bbox_transform_matrix() const
{
    QMatrix4x4 m;
    m.rotate(tilt, QVector3D(1, 0, 0));
    m.rotate(yaw,  QVector3D(0, 0, 1));
    m.scale(scale);
    m.translate(-center);
    return m;
}

QMatrix4x4 Canvas::view_matrix() const
{
    QMatrix4x4 m;
    if (width() > height())
    {
        m.scale(-height() / float(width()), 1, 0.5);
    }
    else
    {
        m.scale(-1, width() / float(height()), 0.5);
    }
    m.scale(zoom * meshScale, zoom * meshScale, 1);
    return m;
}

void Canvas::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton ||
        event->button() == Qt::RightButton)
    {
        mouse_pos = event->pos();
        setCursor(Qt::ClosedHandCursor);
    }
}

void Canvas::mouseReleaseEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton ||
        event->button() == Qt::RightButton)
    {
        unsetCursor();
    }
}

void Canvas::mouseMoveEvent(QMouseEvent* event)
{
    auto p = event->pos();
    auto d = p - mouse_pos;
    if (event->buttons() & Qt::LeftButton)
    {
        yaw = fmod(yaw - d.x(), 360);
        tilt = tilt - d.y();
        update();
    }
    else if (event->buttons() & Qt::RightButton)
    {
        center = transform_matrix().inverted() *
                 view_matrix().inverted() *
                 QVector3D(-d.x() / (0.5*width()),
                            d.y() / (0.5*height()), 0);
        update();
    }
    mouse_pos = p;
}

void Canvas::wheelEvent(QWheelEvent *event)
{
    // Find GL position before the zoom operation
    // (to zoom about mouse cursor)
    auto p = event->pos();
    QVector3D v(1 - p.x() / (0.5*width()),
                p.y() / (0.5*height()) - 1, 0);
    QVector3D a = transform_matrix().inverted() *
                  view_matrix().inverted() * v;

    if (event->delta() < 0)
    {
        for (int i=0; i > event->delta(); --i)
            zoom *= 1.001;
    }
    else if (event->delta() > 0)
    {
        for (int i=0; i < event->delta(); ++i)
            zoom /= 1.001;
    }

    // Then find the cursor's GL position post-zoom and adjust center.
    QVector3D b = transform_matrix().inverted() *
                  view_matrix().inverted() * v;
    center += b - a;
    update();
}
