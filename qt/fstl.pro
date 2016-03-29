QT += core gui opengl widgets

TARGET = fstl
TEMPLATE = app

# Bump optimization up to -O3 in release builds
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

INCLUDEPATH += /usr/local/include /usr/local/include/eigen3 ../lib/libigl/include

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
    ../src/clipper.cpp \
    ../src/backdrop.cpp

HEADERS  += \
    ../src/app.h\
    ../src/canvas.h \
    ../src/editorpanel.h \
    ../src/mesh.h \
    ../src/glmesh.h \
    ../src/loader.h \
    ../src/window.h \
    ../src/clipper.hpp \
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
