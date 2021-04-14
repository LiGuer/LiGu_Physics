#include "controlwindow.h"
QString ButtonStyle1 = "QPushButton{font-family:'times new roman';font-size:20px; "
"background:rgba(253,253,253,255);border:2px rgb(200,200,200);border-style:inset;}"
"QPushButton:pressed{background:rgba(85,170,255,100);border:2px rgb(85,170,255);border-style:inset;}";
controlwindow::controlwindow()
{
    QList<QScreen *> list_screen =  QGuiApplication::screens();  //多显示器
    QRect screenSize = list_screen.at(0)->geometry();

    ControlWidget->resize(QSize(screenSize.width() / 4, screenSize.height()/2 ));
    ControlWidget->setWindowTitle("ControlWindow");
    ControlWidget->setStyleSheet("background:rgba(240,240,240,255);");
    ControlWidget->show();

    const int ButtonSize = 30;
    suspendButton->setText("Suspend");
    initButton->setText("Init");
    suspendButton->setGeometry(0,0,ButtonSize*3,ButtonSize);
    initButton->setGeometry(0,ButtonSize+2,ButtonSize*3,ButtonSize);
    suspendButton->setStyleSheet(ButtonStyle1);
    initButton->setStyleSheet(ButtonStyle1);

    BoundCondPointButton->setText("enable");
    BoundCondPointButton->setGeometry(0,12*ButtonSize+3*2,ButtonSize*3,ButtonSize);
    BoundCondPointButton->setStyleSheet(ButtonStyle1);
    BoundCondPointEdit->setGeometry(0,2*ButtonSize+2*2,ButtonSize*10,ButtonSize*10);
    BoundCondPointEdit->setStyleSheet("font-family:'times new roman';font-size:20px;border:2px rgb(85,170,255);background:rgba(255,255,255,255);");
}
