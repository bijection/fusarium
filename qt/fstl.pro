QT += core gui opengl widgets

TARGET = fstl
TEMPLATE = app

# Bump optimization up to -O3 in release builds
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

INCLUDEPATH += /usr/local/include ../libigl/include /usr/local/include/eigen3

LIBS += -L"/usr/local/lib" -lCGAL -lboost_thread-mt -lgmp -lmpfr -lgsl -lgslcblas -lm -lboost_system

SOURCES += \
    ../src/app.cpp\
    ../src/main.cpp\
    ../src/editorpanel.cpp\
    ../src/canvas.cpp \
    ../src/mesh.cpp \
    ../src/glmesh.cpp \
    ../src/loader.cpp \
    ../src/window.cpp \
    ../src/backdrop.cpp \
    ../vecmath/Matrix2f.cpp \
    ../vecmath/Matrix3f.cpp \
    ../vecmath/Matrix4f.cpp \
    ../vecmath/Quat4f.cpp \
    ../vecmath/Vector2f.cpp \
    ../vecmath/Vector3f.cpp \
    ../vecmath/Vector4f.cpp

HEADERS  += \
    ../vecmath/Matrix2f.h \
    ../vecmath/Matrix3f.h \
    ../vecmath/Matrix4f.h \
    ../vecmath/Quat4f.h \
    ../vecmath/vecmath.h \
    ../vecmath/Vector2f.h \
    ../vecmath/Vector3f.h \
    ../vecmath/Vector4f.h \
    ../src/app.h\
    ../src/canvas.h \
    ../src/editorpanel.h \
    ../src/mesh.h \
    ../src/glmesh.h \
    ../src/loader.h \
    ../src/window.h \
    ../src/backdrop.h

CONFIG += c++11

RESOURCES += \
    qt.qrc \
    ../gl/gl.qrc

macx {
    QMAKE_INFO_PLIST = ../app/Info.plist
    ICON = ../app/fstl.icns
}

win32 {
    RC_FILE = ../exe/fstl.rc
}

static {
    CONFIG += static
}
