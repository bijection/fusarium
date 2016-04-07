#include <QtWidgets>
#include <iostream>

#include "editorpanel.h"

EditorPanel::EditorPanel(QWidget *parent)
    : QWidget(parent)
{
    QFormLayout *grid = new QFormLayout;
    grid->addRow(createOrientationGroup());
    grid->addRow(createScaleGroup());
    grid->addRow(createParametersGroup());
    grid->addRow(createViewGroup());
    setLayout(grid);

    QObject::connect(xRotateSlider, &QSlider::valueChanged, this, &EditorPanel::updateXRotate);
    QObject::connect(yRotateSlider, &QSlider::valueChanged, this, &EditorPanel::updateYRotate);
    QObject::connect(zRotateSlider, &QSlider::valueChanged, this, &EditorPanel::updateZRotate);
    QObject::connect(scaleSlider, &QSlider::valueChanged, this, &EditorPanel::updateMeshScale);
    QObject::connect(zThicknessSlider, &QSlider::valueChanged, this, &EditorPanel::updateZThickness);
    QObject::connect(moldWidthSlider, &QSlider::valueChanged, this, &EditorPanel::updateMoldWidth);
    QObject::connect(connectorsSlider, &QSlider::valueChanged, this, &EditorPanel::updateConnectors);
}

QGroupBox *EditorPanel::createOrientationGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Orientation"));
    QFormLayout *layout = new QFormLayout;
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

    xRotateSlider = new QSlider(Qt::Horizontal);
    yRotateSlider = new QSlider(Qt::Horizontal);
    zRotateSlider = new QSlider(Qt::Horizontal);

    xRotateSlider->setRange(-180, 180);
    xRotateSlider->setTickPosition(QSlider::TicksBelow);
    xRotateSlider->setTickInterval(30);
    yRotateSlider->setRange(-180, 180);
    yRotateSlider->setTickPosition(QSlider::TicksBelow);
    yRotateSlider->setTickInterval(30);
    zRotateSlider->setRange(-180, 180);
    zRotateSlider->setTickPosition(QSlider::TicksBelow);
    zRotateSlider->setTickInterval(30);

    optimizeBtn = new QPushButton(tr("Optimize"));

    xRotateLabel = new QLabel(tr("0°"));
    xRotateLabel->setMinimumSize(35,1);
    yRotateLabel = new QLabel(tr("0°"));
    yRotateLabel->setMinimumSize(35,1);
    zRotateLabel = new QLabel(tr("0°"));
    zRotateLabel->setMinimumSize(35,1);

    QBoxLayout *xRotate = new QBoxLayout(QBoxLayout::RightToLeft);
    xRotate->addWidget(xRotateLabel);
    xRotate->addWidget(xRotateSlider);

    QBoxLayout *yRotate = new QBoxLayout(QBoxLayout::RightToLeft);
    yRotate->addWidget(yRotateLabel);
    yRotate->addWidget(yRotateSlider);

    QBoxLayout *zRotate = new QBoxLayout(QBoxLayout::RightToLeft);
    zRotate->addWidget(zRotateLabel);
    zRotate->addWidget(zRotateSlider);

    layout->addRow(new QLabel(tr("X rotation")), xRotate);
    layout->addRow(new QLabel(tr("Y rotation")), yRotate);
    layout->addRow(new QLabel(tr("Z rotation")), zRotate);
    layout->addRow(optimizeBtn);

    groupBox->setLayout(layout);

    return groupBox;
}

QGroupBox *EditorPanel::createScaleGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Dimensions"));
    QFormLayout *layout = new QFormLayout;
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

	scaleSlider = new QSlider(Qt::Horizontal);
    scaleSlider->setRange(-90, 90);
    scaleSlider->setTickPosition(QSlider::TicksBelow);
    scaleSlider->setTickInterval(20);

    scaleLabel = new QLabel(tr("1x"));
    scaleLabel->setMinimumSize(35,1);

    QBoxLayout *scaleRow = new QBoxLayout(QBoxLayout::RightToLeft);
    scaleRow->addWidget(scaleLabel);
    scaleRow->addWidget(scaleSlider);

    width = new QLabel();
    height = new QLabel();
    depth = new QLabel();

    layout->addRow(new QLabel(tr("Width")), width);
    layout->addRow(new QLabel(tr("Height")), height);
    layout->addRow(new QLabel(tr("Depth")), depth);
    layout->addRow(new QLabel(tr("Scale by")), scaleRow);

    groupBox->setLayout(layout);
    return groupBox;
}

QGroupBox *EditorPanel::createParametersGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Mold Parameters"));
    QFormLayout *layout = new QFormLayout;
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

    moldCombo = new QComboBox();
    moldCombo->addItem(tr("Top"));
    moldCombo->addItem(tr("Bottom"));

    zThickness = 2;
	zThicknessSlider = new QSlider(Qt::Horizontal);
    zThicknessSlider->setRange(10, 30);
    zThicknessSlider->setValue(zThickness * 10);
    zThicknessSlider->setTickPosition(QSlider::TicksBelow);
    zThicknessSlider->setTickInterval(5);

    zThicknessLabel = new QLabel(tr("2.0") + unit);
    zThicknessLabel->setMinimumSize(50,1);

    QBoxLayout *zThicknessRow = new QBoxLayout(QBoxLayout::RightToLeft);
    zThicknessRow->addWidget(zThicknessLabel);
    zThicknessRow->addWidget(zThicknessSlider);

    moldWidth = 19;
    moldWidthSlider = new QSlider(Qt::Horizontal);
    moldWidthSlider->setRange(0, 400);
    moldWidthSlider->setValue(moldWidth * 10);
    moldWidthSlider->setTickPosition(QSlider::TicksBelow);
    moldWidthSlider->setTickInterval(40);

    moldWidthLabel = new QLabel(tr("19.0") + unit);
    moldWidthLabel->setMinimumSize(50,1);

    QBoxLayout *moldWidthRow = new QBoxLayout(QBoxLayout::RightToLeft);
    moldWidthRow->addWidget(moldWidthLabel);
    moldWidthRow->addWidget(moldWidthSlider);

    connectorSpacing = 25;
    connectorsSlider = new QSlider(Qt::Horizontal);
    connectorsSlider->setRange(10, 50);
    connectorsSlider->setValue(connectorSpacing);
    connectorsSlider->setTickPosition(QSlider::TicksBelow);
    connectorsSlider->setTickInterval(5);

    connectorsLabel = new QLabel(tr("25") + unit);
    connectorsLabel->setMinimumSize(50,1);

    QBoxLayout *connectorsRow = new QBoxLayout(QBoxLayout::RightToLeft);
    connectorsRow->addWidget(connectorsLabel);
    connectorsRow->addWidget(connectorsSlider);

    generateMoldBtn = new QPushButton(tr("Generate Mold"));

    layout->addRow(new QLabel(tr("Mold")), moldCombo);
    layout->addRow(new QLabel(tr("Z thickness")), zThicknessRow);
    layout->addRow(new QLabel(tr("Mold width")), moldWidthRow);
    layout->addRow(new QLabel(tr("Spacing")), connectorsRow);
    layout->addRow(generateMoldBtn);

    groupBox->setLayout(layout);
    return groupBox;
}

QGroupBox *EditorPanel::createViewGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Display Options"));
    QFormLayout *layout = new QFormLayout;
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

    modelView = new QCheckBox(tr("Model"));
    bboxView = new QCheckBox(tr("Bounding Box"));
    moldView = new QCheckBox(tr("Mold"));

    modelView->setChecked(true);
    bboxView->setChecked(true);
    moldView->setChecked(true);

    layout->addRow(modelView);
    layout->addRow(bboxView);
    layout->addRow(moldView);
    groupBox->setLayout(layout);
    return groupBox;
}

void EditorPanel::updateBboxLabels(float xNew, float yNew, float zNew) {
    x = xNew;
    y = yNew;
    z = zNew;
    width->setText(QString::number(x, 'f', 2) + unit);
    depth->setText(QString::number(y, 'f', 2) + unit);
    height->setText(QString::number(z, 'f', 2) + unit);
}

void EditorPanel::updateOrientation(float xRotate, float yRotate, float zRotate) {
    // convert floats to closest ints between -180 and 180
    int xRotateInt = (xRotate) >= 0 ? (float)(xRotate+0.5) : (float)(xRotate-0.5);
    xRotateInt = ((xRotateInt + 180) % 360) - 180;
    xRotateSlider->setValue(xRotateInt);
    updateXRotate(xRotateInt);

    int yRotateInt = (yRotate) >= 0 ? (float)(yRotate+0.5) : (float)(yRotate-0.5);
    yRotateInt = ((yRotateInt + 180) % 360) - 180;
    yRotateSlider->setValue(yRotateInt);
    updateYRotate(yRotateInt);

    int zRotateInt = (zRotate) >= 0 ? (float)(zRotate+0.5) : (float)(zRotate-0.5);
    zRotateInt = ((zRotateInt + 180) % 360) - 180;
    zRotateSlider->setValue(zRotateInt);
    updateZRotate(zRotateInt);
}


void EditorPanel::updateXRotate(int deg) {
    xRotateLabel->setText(QString::number(deg)+"°");
}

void EditorPanel::updateYRotate(int deg) {
    yRotateLabel->setText(QString::number(deg)+"°");
}

void EditorPanel::updateZRotate(int deg) {
    zRotateLabel->setText(QString::number(deg)+"°");
}

void EditorPanel::updateMeshScale(int factor) {
    float meshScale;
    if (factor > 0) {
        meshScale = factor/10.0f + 1.0f;
    } else if (factor < 0) {
        meshScale = 1 / (-factor/10.0f + 1.0f);
    } else {
        meshScale = 1;
    }
    scaleLabel->setText(QString::number(meshScale, 'f' ,2) + "x");
}

void EditorPanel::updateZThickness(int thickness) {
    zThickness = thickness / 10.0;
    zThicknessLabel->setText(QString::number(zThickness) + unit);
}

void EditorPanel::updateMoldWidth(int width) {
    moldWidth = width / 10.0;
    moldWidthLabel->setText(QString::number(moldWidth) + unit);
}

void EditorPanel::updateConnectors(int num) {
    connectorSpacing = num;
    connectorsLabel->setText(QString::number(connectorSpacing) + unit);
}

void EditorPanel::updateBboxUnits(int index) {
    unit = (index == 0) ? tr(" mm") : tr(" inches");
    width->setText(QString::number(x, 'f', 2) + unit);
    depth->setText(QString::number(y, 'f', 2) + unit);
    height->setText(QString::number(z, 'f', 2) + unit);
}

void EditorPanel::resetControls() {
    xRotateSlider->setValue(0);
    yRotateSlider->setValue(0);
    zRotateSlider->setValue(0);
    scaleSlider->setValue(0);
}
