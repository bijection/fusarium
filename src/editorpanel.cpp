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

    layout->addRow(new QLabel(tr("X rotation")), xRotateSlider);
    layout->addRow(new QLabel(tr("Y rotation")), yRotateSlider);
    layout->addRow(new QLabel(tr("Z Rotation")), zRotateSlider);
    layout->addRow(new QPushButton(tr("Optimize")));

    groupBox->setLayout(layout);

    return groupBox;
}

QGroupBox *EditorPanel::createScaleGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Dimensions"));
    QFormLayout *layout = new QFormLayout;
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

	QSlider* scale = new QSlider(Qt::Horizontal);

    width = new QLabel();
    height = new QLabel();
    depth = new QLabel();

    layout->addRow(new QLabel(tr("Width")), width);
    layout->addRow(new QLabel(tr("Height")), height);
    layout->addRow(new QLabel(tr("Depth")), depth);
    layout->addRow(new QLabel(tr("Scale by")), scale);

    groupBox->setLayout(layout);
    return groupBox;
}

QGroupBox *EditorPanel::createParametersGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Mold Parameters"));
    QFormLayout *layout = new QFormLayout;
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

	QSlider* zThickness = new QSlider(Qt::Horizontal);
	QSlider* moldWidth = new QSlider(Qt::Horizontal);
	QSlider* connectors = new QSlider(Qt::Horizontal);
	QPushButton* generateMold = new QPushButton(tr("Generate Mold"));
	// generateMold->setMinimumHeight(10);

    layout->addRow(new QLabel(tr("Z Thickness")), zThickness);
    layout->addRow(new QLabel(tr("Mold Width")), moldWidth);
    layout->addRow(new QLabel(tr("Connectors")), connectors);
    layout->addRow(generateMold);

    groupBox->setLayout(layout);
    return groupBox;
}

QGroupBox *EditorPanel::createViewGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Display Options"));
    QFormLayout *layout = new QFormLayout;
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

    QCheckBox* modelView = new QCheckBox(tr("Model"));
    QCheckBox* axesView = new QCheckBox(tr("Axes"));
    QCheckBox* bboxView = new QCheckBox(tr("Bounding Box"));
    QCheckBox* moldView = new QCheckBox(tr("Mold"));

    modelView->setChecked(true);
    bboxView->setChecked(true);

    layout->addRow(modelView);
    layout->addRow(axesView);
    layout->addRow(bboxView);
    layout->addRow(moldView);
    groupBox->setLayout(layout);
    return groupBox;
}

void EditorPanel::updateBboxLabels(float x, float y, float z) {
    width->setText(QString("%1 mm").arg(x));
    height->setText(QString("%1 mm").arg(y));
    depth->setText(QString("%1 mm").arg(z));
}

void EditorPanel::resetControls() {
    xRotateSlider->setValue(0);
    yRotateSlider->setValue(0);
    zRotateSlider->setValue(0);
}
