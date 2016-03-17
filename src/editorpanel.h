#include <QWidget>
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

public slots:
    void updateBboxLabels(float x, float y, float z);
    void updateBboxUnits(int index);
    void resetControls();

private:
    QGroupBox *createOrientationGroup();
    QGroupBox *createScaleGroup();
    QGroupBox *createParametersGroup();
    QGroupBox *createViewGroup();
    QLabel *width;
    QLabel *height;
    QLabel *depth;
    float x, y, z;
    QString unit = tr(" inches");
    QComboBox* unitCombo;
};