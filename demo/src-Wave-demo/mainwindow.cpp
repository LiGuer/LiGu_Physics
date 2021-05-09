#include "mainwindow.h"
MainWindow::MainWindow(QWidget *parent): QMainWindow(parent)
{
    QSize screenSize = graph->screen()->size();
    screenSize = QSize(screenSize.width() / 2, screenSize.height() / 2);
    this->setMinimumSize(screenSize);
    this->setMaximumSize(screenSize);
    this->setWindowTitle("LiguのWave");
    this->setWindowOpacity(0.9);
    this->setWindowFlags(Qt::FramelessWindowHint);
    container->setGeometry(0,0,screenSize.width(), screenSize.height());
    wave_t.init();
    ControlWindow->BoundCondPointEdit->setText("250,125,1\n250,375,1");
    setBoundCond();
    getView();
    modifier->initSurface();

    connect(timer_time, SIGNAL(timeout()), this, SLOT(getView()));
    connect(ControlWindow->suspendButton, SIGNAL(clicked()), this, SLOT(suspendView()));
    connect(ControlWindow->initButton, SIGNAL(clicked()), this, SLOT(initView()));
    connect(ControlWindow->BoundCondPointButton, SIGNAL(clicked()), this, SLOT(setBoundCond()));

    timer_time->start(100);
}
MainWindow::~MainWindow()
{}

void MainWindow::getView()
{
    if(enable){
        wave_t.boundCond(epoch++);
        wave_t.simulate();
        modifier->ViewSurface(wave_t.N, &wave_t.Map[0][0]);
    }
    timer_time->start(10);
}
void MainWindow::suspendView()
{
    enable = enable == 1 ? 0 : 1;
}
void MainWindow::initView()
{
    wave_t.init();
    epoch = 0;
}
void MainWindow::setBoundCond()
{
    wave_t.clearBoundCondPoint();
    QString str = ControlWindow->BoundCondPointEdit->toPlainText();
    QStringList strlist = str.split("\n");
    for(int i=0;i<strlist.count();i++){
        int x = strlist[i].section(",", 0, 0).toInt();
        int y = strlist[i].section(",", 1, 1).toInt();
        double A = strlist[i].section(",", 2, 2).toDouble();
        wave_t.addBoundCondPoint(x,y,A,1);
    }
}

