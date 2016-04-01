#include <QMouseEvent>
#include <QDebug>

#include <cmath>
#include <QtDebug>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "canvas.h"
#include "backdrop.h"
#include "glmesh.h"
#include "mesh.h"

Canvas::Canvas(const QGLFormat& format, EditorPanel* editorPanel, QWidget *parent)
    : QGLWidget(format, parent), mesh(NULL),
      scale(1), zoom(1), tilt(90), yaw(0), status(" ")
{
    panel = editorPanel;
}

Canvas::~Canvas()
{
    delete glmesh;
}

void Canvas::load_mesh(Mesh* m)
{
    mesh = m;
    glmesh = new GLMesh(mesh);
    glmoldGenerated = false;

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
    QMatrix4x4 m;
    m.rotate(meshRotateX, QVector3D(1, 0, 0));
    m.rotate(meshRotateY, QVector3D(0, 1, 0));
    m.rotate(meshRotateZ, QVector3D(0, 0, 1));
    mesh->calculateProjectedArea(QVector3D(0, 0, 1) * m);

    emit updatedBbox(
        (bbox.xmax - bbox.xmin) * meshScale,
        (bbox.ymax - bbox.ymin) * meshScale,
        (bbox.zmax - bbox.zmin) * meshScale);
}

void Canvas::updateModelView(int state) {
    modelView = (state > 0);
    update();
}


void Canvas::updateBboxView(int state) {
    bboxView = (state > 0);
    update();
}


void Canvas::updateMoldView(int state) {
    moldView = (state > 0);
    update();
}

void Canvas::generateMold() {
    glmoldGenerated = false;
    update();

    QMatrix4x4 m;
    m.rotate(meshRotateX, QVector3D(1, 0, 0));
    m.rotate(meshRotateY, QVector3D(0, 1, 0));
    m.rotate(meshRotateZ, QVector3D(0, 0, 1));

    glMold = new GLMesh(mesh->generateMold(QVector3D(0,0,1) * m, meshScale,
        panel->zThickness, panel->moldWidth, panel->connectors,
        (panel->moldCombo->currentIndex() == 0)));
    glmoldGenerated = true;
    update();
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
    if (modelView) {
        mesh_shader.bind();

        // Load the transform and view matrices into the shader
        glUniformMatrix4fv(
                    mesh_shader.uniformLocation("transform_matrix"),
                    1, GL_FALSE, transform_matrix().data());
        glUniformMatrix4fv(
                    mesh_shader.uniformLocation("view_matrix"),
                    1, GL_FALSE, view_matrix().data());

        // Compensate for z-flattening when zooming
        glUniform1f(mesh_shader.uniformLocation("zoom"), 1/(zoom));

        // Find and enable the attribute location for vertex position
        const GLuint vp = mesh_shader.attributeLocation("vertex_position");
        glEnableVertexAttribArray(vp);

        // Then draw the mesh with that vertex position
        glmesh->draw(vp);

        glDisableVertexAttribArray(vp);
        mesh_shader.release();
    }

    if (moldView && glmoldGenerated) {
        mesh_shader.bind();
        glUniformMatrix4fv(
                    mesh_shader.uniformLocation("transform_matrix"),
                    1, GL_FALSE, mold_transform_matrix().data());
        glUniformMatrix4fv(
                    mesh_shader.uniformLocation("view_matrix"),
                    1, GL_FALSE, view_matrix().data());
        glUniform1f(mesh_shader.uniformLocation("zoom"), 1/(zoom));
        const GLuint vp = mesh_shader.attributeLocation("vertex_position");
        glEnableVertexAttribArray(vp);
        glMold->draw(vp);

        glDisableVertexAttribArray(vp);
        mesh_shader.release();
    }

    if (bboxView) {
        line_shader.bind();
        glUniformMatrix4fv(
                    line_shader.uniformLocation("transform_matrix"),
                    1, GL_FALSE, bbox_transform_matrix().data());
        glUniformMatrix4fv(
                    line_shader.uniformLocation("view_matrix"),
                    1, GL_FALSE, view_matrix().data());
        glmesh->drawBoundingBox();

        line_shader.release();
        line_shader.bind();
        glUniformMatrix4fv(
                    line_shader.uniformLocation("transform_matrix"),
                    1, GL_FALSE, transform_matrix().data());
        glUniformMatrix4fv(
                    line_shader.uniformLocation("view_matrix"),
                    1, GL_FALSE, view_matrix().data());

        glmesh->drawOverhangLine();

        line_shader.release();
    }
}

QMatrix4x4 Canvas::transform_matrix() const
{
    QMatrix4x4 m;
    m.rotate(tilt, QVector3D(1, 0, 0));
    m.rotate(yaw,  QVector3D(0, 0, 1));
    m.scale(scale);
    m.translate(-center);
    m.scale(meshScale);
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
    m.scale(meshScale);
    return m;
}


QMatrix4x4 Canvas::mold_transform_matrix() const
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
    m.scale(zoom, zoom, 1);
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

double calcArea (const gsl_vector *v, void *params) {
    double meshRotateX, meshRotateY, meshRotateZ;
    Mesh *p = (Mesh *) params;

    meshRotateX = gsl_vector_get(v, 0);
    meshRotateY = gsl_vector_get(v, 1);
    meshRotateZ = gsl_vector_get(v, 2);

    QMatrix4x4 m;
    m.rotate(meshRotateX, QVector3D(1, 0, 0));
    m.rotate(meshRotateY, QVector3D(0, 1, 0));
    m.rotate(meshRotateZ, QVector3D(0, 0, 1));
    return -(p[0]).calculateProjectedArea(QVector3D(0, 0, 1) * m);
}

void Canvas::optimizeMesh() {
    Mesh par[1] = {*mesh};

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point */
    x = gsl_vector_alloc (3);
    gsl_vector_set (x, 0, meshRotateX);
    gsl_vector_set (x, 1, meshRotateY);
    gsl_vector_set (x, 2, meshRotateZ);

    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (3);
    gsl_vector_set_all (ss, 1.0);

    /* Initialize method and iterate */
    minex_func.n = 3;
    minex_func.f = calcArea;
    minex_func.params = par;

    s = gsl_multimin_fminimizer_alloc (T, 3);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status) break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);

        // uncomment this out to see the gory details of the optimization
        // if (status == GSL_SUCCESS) {
        //     printf ("converged to minimum at\n");
        // }

        // printf ("%5d %3.3f %3.3f %3.3f f() = %7.3f size = %.3f\n",
        //       iter,
        //       gsl_vector_get (s->x, 0),
        //       gsl_vector_get (s->x, 1),
        //       gsl_vector_get (s->x, 2),
        //       s->fval, size);
    }
    while (status == GSL_CONTINUE && iter < 100);

    setMeshRotateX(gsl_vector_get (s->x, 0));
    setMeshRotateY(gsl_vector_get (s->x, 1));
    setMeshRotateZ(gsl_vector_get (s->x, 2));

    emit updatedOrientation(meshRotateX, meshRotateY, meshRotateZ);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
}