#ifndef CANVAS_H
#define CANVAS_H

#include <QWidget>
#include <QtOpenGL/QGLWidget>
#include <QtOpenGL/QGLFunctions>
#include <QtOpenGL/QGLShaderProgram>
#include <QMatrix4x4>

#include "editorpanel.h"

class GLMesh;
class Mesh;
class Backdrop;

class Canvas : public QGLWidget, protected QGLFunctions
{
    Q_OBJECT

public:
    Canvas(const QGLFormat& format, EditorPanel* editorPanel, QWidget* parent=0);

    void initializeGL();
    void paintEvent(QPaintEvent* event);
    ~Canvas();

signals:
    void updatedBbox(float x, float y, float z);
    void updatedOrientation(float x, float y, float z);

public slots:
    void set_status(const QString& s);
    void clear_status();
    void load_mesh(Mesh* m);
    void setMeshRotateX(int degrees);
    void setMeshRotateY(int degrees);
    void setMeshRotateZ(int degrees);
    void setMeshScale(int factor);
    void optimizeMesh();
    void generateMold();
    void exportSTL();
    void updateModelView(int state);
    void updateBboxView(int state);
    void updateMoldView(int state);

protected:
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);


private:
    void draw_mesh();
    void updateMeshBbox();

    QMatrix4x4 transform_matrix() const;
    QMatrix4x4 view_matrix() const;
    QMatrix4x4 mold_transform_matrix() const;
    QMatrix4x4 bbox_transform_matrix() const;

    QGLShaderProgram mesh_shader;
    QGLShaderProgram line_shader;
    QGLShaderProgram quad_shader;

    GLMesh* glmesh;
    GLMesh* glMold;
    Mesh* mesh;
    Mesh* mold_mesh;
    Backdrop* backdrop;
    EditorPanel* panel;

    QVector3D center;
    float scale;
    float zoom;
    float tilt;
    float yaw;

    float meshRotateX, meshRotateY, meshRotateZ;
    float meshScale = 1;
    bool glmoldGenerated = false;
    bool modelView = true;
    bool bboxView = true;
    bool moldView = true;

    QPoint mouse_pos;
    QString status;
};

#endif // CANVAS_H
