#pragma once
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>  
#include <QDialog>
#include "ui_result_show.h"
#include "calculation.h"
#include <vector>
#include<fstream>
#include<QFile>

using namespace cv;


class Result_Show : public QDialog
{
	Q_OBJECT

public:
	Result_Show(const int d, const Mat &l = Mat(), const cv::Mat &r = cv::Mat(),  QWidget *parent = nullptr);
	~Result_Show();
	cv::Mat left_img, right_img;
	QImage mat2QImage(const cv::Mat& mat);
	QImage rst_q;
	vector<Match_Pts> pts;
	int delta_x;
	void on_save_points_clicked();//±£¥Ê∆•≈‰µ„
	void on_save_result_clicked();
	void write_points();
private:
	Ui::Result_ShowClass *ui;

};
