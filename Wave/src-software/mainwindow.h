#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "surfacegraph.h"
#include <QApplication>
#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>
#include <QtGui/QPainter>
#include <QtGui/QScreen>
#include <QMainWindow>
#include <QTimer>
#include <QTimer>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include "wave.h"
#include "controlwindow.h"
class MainWindow : public QMainWindow
{
    Q_OBJECT
private:
    int epoch=0;
    QTimer *timer_time= new QTimer(this);
    bool enable=1;
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
    Q3DSurface *graph = new Q3DSurface();
    QWidget *container = QWidget::createWindowContainer(graph,this);
    controlwindow *ControlWindow = new controlwindow();
    SurfaceGraph *modifier = new SurfaceGraph(graph);
    wave wave_t;

public slots:
    void getView();
    void suspendView();
    void initView();
    void setBoundCond();
};

#endif // MAINWINDOW_H
