QT       += core gui
QT += datavisualization
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = LiGu_Wave

TEMPLATE = app

DEFINES += QT_DEPRECATED_WARNINGS


CONFIG += c++11

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    surfacegraph.cpp \
    wave.cpp \
    controlwindow.cpp

HEADERS += \
        mainwindow.h \
    surfacegraph.h \
    wave.h \
    controlwindow.h

qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
