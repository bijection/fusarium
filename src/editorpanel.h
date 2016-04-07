#include <QWidget>
#include <QtWidgets>
#include <QGridLayout>
#include <QGroupBox>

class EditorPanel : public QWidget
{

public:
    explicit EditorPanel(QWidget *parent);
    QSlider *xRotateSlider;
    QSlider *yRotateSlider;
    QSlider *zRotateSlider;
    QSlider *scaleSlider;
    QSlider *zThicknessSlider;
    QSlider *moldWidthSlider;
    QSlider *connectorsSlider;
    QPushButton *optimizeBtn;
    QPushButton *generateMoldBtn;
    QComboBox *moldCombo;
    QCheckBox *modelView;
    QCheckBox *bboxView;
    QCheckBox *moldView;

    float zThickness;
    float moldWidth;
    int connectorSpacing;

public slots:
    void updateBboxLabels(float x, float y, float z);
    void updateBboxUnits(int index);
    void updateOrientation(float x, float y, float z);
    void updateXRotate(int deg);
    void updateYRotate(int deg);
    void updateZRotate(int deg);
    void updateMeshScale(int factor);
    void updateZThickness(int thickness);
    void updateMoldWidth(int width);
    void updateConnectors(int num);
    void resetControls();

private:
    QGroupBox *createOrientationGroup();
    QGroupBox *createScaleGroup();
    QGroupBox *createParametersGroup();
    QGroupBox *createViewGroup();
    QLabel *width;
    QLabel *height;
    QLabel *depth;
    QLabel *xRotateLabel;
    QLabel *yRotateLabel;
    QLabel *zRotateLabel;
    QLabel *scaleLabel;
    QLabel *zThicknessLabel;
    QLabel *moldWidthLabel;
    QLabel *connectorsLabel;
    float x, y, z;
    QString unit = tr(" mm");
};