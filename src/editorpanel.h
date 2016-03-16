#include <QWidget>
#include <QGridLayout>
#include <QGroupBox>

class EditorPanel : public QWidget
{

public:
    explicit EditorPanel(QWidget *parent);

private:
    QGroupBox *createOrientationGroup();
    QGroupBox *createScaleGroup();
    QGroupBox *createParametersGroup();
    QGroupBox *createViewGroup();
};