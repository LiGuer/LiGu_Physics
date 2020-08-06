#include "controlwindow.h"

controlwindow::controlwindow()
{
    QList<QScreen *> list_screen =  QGuiApplication::screens();  //多显示器
    QRect screenSize = list_screen.at(0)->geometry();

    ControlWidget->resize(QSize(screenSize.width() / 4, screenSize.height()/2 ));
    ControlWidget->setWindowTitle("ControlWindow");
    ControlWidget->show();

    const int ButtonSize = 20;
    suspendButton->setText("Suspend");
    initButton->setText("Init");
    suspendButton->setGeometry(0,0,ButtonSize*5,ButtonSize);
    initButton->setGeometry(0,ButtonSize,ButtonSize*5,ButtonSize);

    BoundCondPointButton->setText("确定");
    BoundCondPointButton->setGeometry(0,22*ButtonSize,ButtonSize*5,ButtonSize);
    BoundCondPointEdit->setGeometry(0,2*ButtonSize,ButtonSize*20,ButtonSize*20);
}
