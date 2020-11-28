#ifndef CONTROLWINDOW_H
#define CONTROLWINDOW_H
#include <QGuiApplication>
#include <QScreen>
#include <QRect>
#include <QWidget>
#include <QPushButton>
#include <QTextEdit>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QVBoxLayout>
class controlwindow
{
public:
    controlwindow();
    QWidget* ControlWidget = new QWidget();
    QHBoxLayout *hLayout = new QHBoxLayout(ControlWidget);
    QPushButton* suspendButton = new QPushButton(ControlWidget);
    QPushButton* initButton = new QPushButton(ControlWidget);
    QPushButton* BoundCondPointButton = new QPushButton(ControlWidget);
    QTextEdit* BoundCondPointEdit = new QTextEdit(ControlWidget);
};

#endif // CONTROLWINDOW_H
