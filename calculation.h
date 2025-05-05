#pragma once
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>  
#include <opencv2/imgproc/imgproc.hpp>  
#include <vector>
#include <cstdlib> // 包含rand()和srand()  
#include <ctime>   // 包含time() 
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace cv;
using namespace std;

struct Match_Pts//同名点
{
	Point left;
	Point right;
	float r;//相关系数
};

class Calculation
{
	//包含功能：寻找特征点、匹配、最小二乘优化、输出结果
	//寻找特征点：给定一张图像，输出特征点坐标集
	//匹配：给定两张图像和各自特征点坐标，输出匹配结果（连线）和转换矩阵
	//最小二乘优化：给定转换矩阵和匹配点坐标集，输出优化后的匹配结果
	//输出结果：包括同名点文件、两张连线结果、转换后的图像

	//特征点部分
	vector<Point> feature_points;//特征点集
	int F_Threshold=35000;//兴趣阈值
	int K_Size=3;//窗口半径
	Mat Mora_Resultimg;//检测结果图
	Mat Interest_Mat;//兴趣矩阵
	
	//匹配部分
	float M_Threshold = 0.92;
	int M_Size = 3;
public:
	Mat LeftImg;
	Mat RightImg;
	Mat MatchImg;//彩色
	//优化部分
	int L_Size = 7;
	float L_Threshold = 0.9;
	vector<Match_Pts> Match_Points;//同名点
	int deltax;
	//算法功能
private:
	double Get_NCC(int lr,int lc, int rr, int rc);
	double Get_NCC(const Mat& leftWindow, const Mat& rightWindow);
	void lqm(Match_Pts& match);
	void CreatResult();
	double a0 = 0, a1 = 1, a2 = 0;
	double b0 = 0, b1 = 0, b2 = 0;
	double h0 = 0, h1 = 1;
	
public:
	void Moravec(Mat srcImage);
	void Match();
	Mat DrawRst();
	Calculation();
	~Calculation();

};

