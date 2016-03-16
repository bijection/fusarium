#include <QtWidgets>

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

    QSlider* xSpin = new QSlider(Qt::Horizontal);
    QSlider* ySpin = new QSlider(Qt::Horizontal);
    QSlider* zSpin = new QSlider(Qt::Horizontal);

    layout->addRow(new QLabel(tr("x°")), xSpin);
    layout->addRow(new QLabel(tr("y°")), ySpin);
    layout->addRow(new QLabel(tr("z°")), zSpin);
    layout->addRow(new QPushButton(tr("Optimize")));

    groupBox->setLayout(layout);
    return groupBox;
}

QGroupBox *EditorPanel::createScaleGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Scale"));
    QFormLayout *layout = new QFormLayout;
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

	QSlider* scale = new QSlider(Qt::Horizontal);

    layout->addRow(new QLabel(tr("Width")), new QLabel(tr("5.5\"")));
    layout->addRow(new QLabel(tr("Height")), new QLabel(tr("6.1\"")));
    layout->addRow(new QLabel(tr("Depth")), new QLabel(tr("2.3\"")));
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

    layout->addRow(new QCheckBox(tr("Model")));
    layout->addRow(new QCheckBox(tr("Axes")));
    layout->addRow(new QCheckBox(tr("Bounding Box")));
    layout->addRow(new QCheckBox(tr("Mold")));
    groupBox->setLayout(layout);
    return groupBox;
}
