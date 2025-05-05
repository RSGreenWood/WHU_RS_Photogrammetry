#pragma once
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>  
#include <opencv2/imgproc/imgproc.hpp>  
#include <vector>
#include <cstdlib> // ����rand()��srand()  
#include <ctime>   // ����time() 
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace cv;
using namespace std;

struct Match_Pts//ͬ����
{
	Point left;
	Point right;
	float r;//���ϵ��
};

class Calculation
{
	//�������ܣ�Ѱ�������㡢ƥ�䡢��С�����Ż���������
	//Ѱ�������㣺����һ��ͼ��������������꼯
	//ƥ�䣺��������ͼ��͸������������꣬���ƥ���������ߣ���ת������
	//��С�����Ż�������ת�������ƥ������꼯������Ż����ƥ����
	//������������ͬ�����ļ����������߽����ת�����ͼ��

	//�����㲿��
	vector<Point> feature_points;//�����㼯
	int F_Threshold=35000;//��Ȥ��ֵ
	int K_Size=3;//���ڰ뾶
	Mat Mora_Resultimg;//�����ͼ
	Mat Interest_Mat;//��Ȥ����
	
	//ƥ�䲿��
	float M_Threshold = 0.92;
	int M_Size = 3;
public:
	Mat LeftImg;
	Mat RightImg;
	Mat MatchImg;//��ɫ
	//�Ż�����
	int L_Size = 7;
	float L_Threshold = 0.9;
	vector<Match_Pts> Match_Points;//ͬ����
	int deltax;
	//�㷨����
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

