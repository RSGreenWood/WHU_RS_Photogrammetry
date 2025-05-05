#pragma once
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>  
#include <QtWidgets/QWidget>
#include <QFileDialog>  
#include <QImage> 
#include <QDebug>
#include <QPainter>
#include <qmessagebox.h>
#include <qmainwindow.h>
#include "ui_photogrammetry.h"
#include <string>
#include <qlabel.h>
#include"result_show.h"

using namespace std;
using namespace cv;


class photogrammetry : public QWidget
{
    Q_OBJECT

public:
    photogrammetry(QWidget *parent = nullptr);
    ~photogrammetry();
private slots:
    void on_Load_Left_Button_clicked();//载入左影像
    void on_Load_Right_Button_clicked();//载入右影像
    void on_Caculation_clicked();//匹配计算
private:
    Ui::photogrammetryClass *ui;

    //QPixmap与Mat互相转换
    Mat QPix2Mat(QPixmap &ori);
    Mat QImage2Mat(const QImage& image);
    QPixmap Mat2QPix(Mat &ori);
    //QLineEdit* GetDelta;
    QLabel* img_depart;
//protected:
public:
    Mat left_img, right_img;
    int delta;
};

