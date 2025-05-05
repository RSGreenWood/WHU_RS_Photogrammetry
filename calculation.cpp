#include "calculation.h"

using namespace cv;
using namespace std;
void Calculation::Moravec(Mat srcImage)
{
	Mora_Resultimg = srcImage.clone();
	const int nRows = srcImage.rows;
	const int nCols = srcImage.cols;
	Interest_Mat = Mat::ones(nRows, nCols, CV_8U);
	for (int i = K_Size; i < nRows - K_Size; i++) {
		for (int j = K_Size; j < nCols - K_Size; j++) {
			int wV1, wV2, wV3, wV4;
			wV1 = wV2 = wV3 = wV4 = 0;
			//计算水平方向的兴趣值
			for (int k = -K_Size; k < K_Size + 1; k++) {
				for (int l = -K_Size; l < K_Size + 1; l++) {
					if (j + l+ 1 >= nCols||i + k + 1 >= nRows|| j + l - 1 < 0) {//出界
						continue;
					}
					wV1 += (srcImage.at<uchar>(i + k, j + l + 1) - srcImage.at<uchar>(i + k, j + l))
						* (srcImage.at<uchar>(i + k, j + l + 1) - srcImage.at<uchar>(i + k, j + l));//计算水平方向的兴趣值
					wV2+= (srcImage.at<uchar>(i + k + 1, j + l) - srcImage.at<uchar>(i + k, j + l))
						* (srcImage.at<uchar>(i + k + 1, j + l) - srcImage.at<uchar>(i + k, j + l));//计算竖直方向的兴趣值
					wV3 += (srcImage.at<uchar>(i + k + 1, j + l + 1) - srcImage.at<uchar>(i + k, j + l))
						* (srcImage.at<uchar>(i + k + 1, j + l + 1) - srcImage.at<uchar>(i + k, j + l));//计算/方向的兴趣值
					wV4 += (srcImage.at<uchar>(i + k + 1, j + l - 1) - srcImage.at<uchar>(i + k, j + l))
						* (srcImage.at<uchar>(i + k + 1, j + l - 1) - srcImage.at<uchar>(i + k, j + l));//计算\方向的兴趣值
				}
			}
			int Interest_Value = min(min(wV1, wV2), min(wV3, wV4));
			if (Interest_Value > F_Threshold) {
				Point the_point(j,i);
				feature_points.push_back(the_point);
				Interest_Mat.at<uchar>(i, j) = 0;
			}
		}
	}
	cvtColor(Mora_Resultimg, Mora_Resultimg, COLOR_GRAY2BGR);
	for (int i = 0; i < feature_points.size(); i++) {
		Scalar randomColour(
			std::rand() % 256, // 随机红色分量  
			std::rand() % 256, // 随机绿色分量  
			std::rand() % 256  // 随机蓝色分量
		);
		circle(Mora_Resultimg, feature_points[i], 5, randomColour, 1);
	}
	namedWindow("Moravec Result", WINDOW_NORMAL);
	imshow("Moravec Result", Mora_Resultimg);
	waitKey(-1);

}

double Calculation::Get_NCC(int lr, int lc, int rr, int rc)
{
	double gLeftAverage = 0;//左影像窗口灰度平均值  
	double gRightAverage = 0;//右影像窗口灰度平均值;
	for (int i = -M_Size; i < 1 + M_Size; i++)
	{
		for (int j = -M_Size; j < 1 + M_Size; j++)
		{ 
			gLeftAverage += LeftImg.at<uchar>(lr + i, lc + j);
			gRightAverage += RightImg.at<uchar>(rr + i, rc + j);
		}
	}
	gLeftAverage /= (1 + M_Size) * (1 + M_Size);
	gRightAverage /= (1 + M_Size) * (1 + M_Size);
	double a = 0;
	double b = 0;
	double c = 0;
	for (int i = -M_Size; i < 1 + M_Size; i++)
	{
		for (int j = -M_Size; j < 1 + M_Size; j++)
		{
			double left_av = LeftImg.at<uchar>(lr + i, lc + j) - gLeftAverage;
			double right_av = RightImg.at<uchar>(rr + i, rc + j) - gRightAverage;
			a += left_av * right_av;
			b += left_av * left_av;
			c += right_av * right_av;
		}
	}
	return a / sqrt(b * c);
}

double Calculation::Get_NCC(const Mat& leftWindow, const Mat& rightWindow)
{
	double gLeftAverage = 0;//左影像窗口灰度平均值  
	double gRightAverage = 0;//右影像窗口灰度平均值;
	for (int i = 0; i < leftWindow.rows; i++)
	{
		for (int j = 0; j < leftWindow.cols; j++)
		{
			gLeftAverage += leftWindow.at<uchar>(i, j);
			gRightAverage += rightWindow.at<uchar>(i, j);
		}
	}
	gLeftAverage /= (leftWindow.rows) * (leftWindow.cols);
	gRightAverage /= (rightWindow.rows) * (rightWindow.cols);
	double a = 0;
	double b = 0;
	double c = 0;
	for (int i = 0; i < leftWindow.rows; i++)
	{
		for (int j = 0; j < leftWindow.cols; j++)
		{
			double left_av = leftWindow.at<uchar>(i, j) - gLeftAverage;
			double right_av = rightWindow.at<uchar>(i, j) - gRightAverage;
			a += left_av * right_av;
			b += left_av * left_av;
			c += right_av * right_av;
		}
	}
	return a / sqrt(b * c);
}

void Calculation::CreatResult()
{
	// 确保至少有一个图像有效
	if (LeftImg.empty() && RightImg.empty())
	{
		cerr << "Error: One or both images are empty." << std::endl;
		return;
	}
	int depth = LeftImg.depth();
	int channels = LeftImg.channels();
	int result_row = max(LeftImg.rows, RightImg.rows);
	MatchImg.create(result_row, LeftImg.cols + RightImg.cols, LeftImg.type());

	// 复制LeftImg到MatchImg的左边
	for (int i = 0; i < LeftImg.rows; i++)
	{
		for (int j = 0; j < LeftImg.cols; j++)
		{
			MatchImg.at<uchar>(i, j) = LeftImg.at<uchar>(i, j);
		}
	}

	// 复制RightImg到MatchImg的右边
	for (int i = 0; i < RightImg.rows; i++)
	{
		for (int j = 0; j < RightImg.cols; j++)
		{
			MatchImg.at<uchar>(i, LeftImg.cols + j) = RightImg.at<uchar>(i, j);
		}
		i;
	}
	// 如果LeftImg和RightImg的高度不同，需要处理未填充的部分
	if (LeftImg.rows < result_row)
	{
		for (int i = LeftImg.rows; i < result_row; i++)
		{
			for (int j = 0; j < LeftImg.cols; j++)
			{
				// 可以选择将这部分填充为黑色或其他颜色
				MatchImg.at<uchar>(i, j) = uchar(0);  // 填充黑色
			}
		}
	}

	if (RightImg.rows < result_row)
	{
		for (int i = RightImg.rows; i < result_row; i++)
		{
			for (int j = 0; j < RightImg.cols; j++)
			{
				// 可以选择将这部分填充为黑色或其他颜色
				MatchImg.at<uchar>(i, LeftImg.cols + j) = uchar(0);  // 填充黑色
			}
		}
	}
}

void Calculation::Match()
{
	CreatResult();

	int Lr = LeftImg.rows;
	int Lc = LeftImg.cols;
	int Rr = RightImg.rows;
	int Rc = RightImg.cols;
	const int window_radius = 15;
	const int offset_x = deltax;
	const int offset_y = -20;
	const int step = 5;
	for (int i = M_Size+window_radius; i < Lr - M_Size-window_radius;  i = i+ step )
	{
		for (int j = M_Size+window_radius; j < Lc - window_radius - M_Size; j = j + step)
		{
			if (i<0 || i>Lr || j<0 || j>Lc) {
				continue;
			}
			if (Interest_Mat.at<uchar>(i, j) == 0)
				//特征点作为模板中心  
			{
				double maxscore = 0;
				Match_Pts temp;
				// 定义左图窗口范围
				Rect left_window_rect(j - window_radius, i - window_radius, M_Size, M_Size);
				Mat leftWindow = LeftImg(left_window_rect);

				// 计算右图搜索范围
				int search_r_start = max(0, i + offset_y - window_radius);
				int search_r_end = min(Rr - M_Size, i + offset_y + window_radius);
				int search_c_start = max(offset_x, j + offset_x - window_radius);
				int search_c_end = min(Rc - M_Size, j + offset_x + window_radius);

				// 在右图中搜索匹配点
				if (search_r_start<0 || search_r_start>Rr || search_r_end<0 || search_r_end>Rr ||
					search_c_start<0 || search_c_start>Rc || search_c_start<0 || search_c_start>Rc) {
					continue;
				}
				for (int r = search_r_start; r <= search_r_end; r += step)
				{
					for (int c = search_c_start; c <= search_c_end; c += step)
					{
						// 定义右图窗口范围
						// 确保右图窗口不会越界
						int right_window_x = std::max(0, c - window_radius);
						int right_window_y = std::max(0, r - window_radius);
						int right_window_width = std::min(M_Size, RightImg.cols - right_window_x);
						int right_window_height = std::min(M_Size, RightImg.rows - right_window_y);

						// 如果窗口大小不足，则跳过
						if (right_window_width != M_Size || right_window_height != M_Size) {
							continue;
						}

						// 定义右图窗口范围
						Rect right_window_rect(right_window_x, right_window_y, right_window_width, right_window_height);
						Mat rightWindow = RightImg(right_window_rect);
						double score = Get_NCC(leftWindow, rightWindow); // 计算相关系数
						if (score > maxscore)
						{
							maxscore = score;//计算相关系数的最大值
							temp.left = Point(j, i);//顺序反着的
							temp.right = Point(c, r);
							temp.r = score;
						}
					}
				}
				if (maxscore > this->M_Threshold) {
					Match_Points.push_back(temp);
				}

				}
			}
		}
	for (auto& mtchpts : Match_Points) {
		lqm(mtchpts);
	}
}

void Calculation::lqm(Match_Pts &match)
{
	int windowsize = L_Size * 2 + 1;
	Mat left_copy, right_copy;
	left_copy = LeftImg.clone();
	right_copy = RightImg.clone();
	//x y相反
	double x1 = match.left.y;
	double y1 = match.left.x;
	double x2 = match.right.y;
	double y2 = match.right.x;
	Rect rectL, rectR;
	Mat winL, winR;
	rectL = Rect(y1 - L_Size, x1 - L_Size, windowsize, windowsize);
	winL = left_copy(rectL);
	winR.create(Size(windowsize, windowsize), CV_8UC1);

	// 设定几何畸变初值
	a0 = x2 - x1;
	a1 = 1;
	a2 = 0;
	b0 = y2 - y1;
	b1 = 0;
	b2 = 1;

	// 设定灰度畸变初值
	h0 = 0;
	h1 = 1;
	double xs = 0.0, ys = 0.0;
	double tempCorIdx, bestCorIdx = 0.0;
	Point2d bestPt;
	for (int iter = 0; iter < 50; iter++) // 设定最大迭代次数不超过50次
	{
		Eigen::MatrixXd A(windowsize* windowsize, 8), L(windowsize* windowsize, 1), x;

		int num = 0;
		double xNumerator = 0.0, yNumerator = 0.0, xDenominator = 0.0, yDenominator = 0.0;

		for (int i = x1 - L_Size; i <= x1 + L_Size; i++)
			for (int j = y1 - L_Size; j <= y1 + L_Size; j++)
			{
				// 几何变形改正
				double m = a0 + a1 * i + a2 * j;
				double n = b0 + b1 * i + b2 * j;

				int I = floor(m);
				int J = floor(n);

				// 如果当前的点在图像的边界附近，就舍弃当前点，因为后面求导会出现问题
				if (I < 1 || I >= right_copy.rows-1 || J < 1 || J >= right_copy.cols-1)
					continue;

				// 重采样：双线性内插
				double pixelValue = (J + 1 - n) * ((I + 1 - m) * right_copy.at<uchar>(I, J) 
					+ (m - I) * right_copy.at<uchar>(I + 1, J)) + (n - J) * ((I + 1 - m) * right_copy.at<uchar>(I, J + 1) 
						+ (m - I) * right_copy.at<uchar>(I + 1, J + 1));

				// 辐射畸变改正
				pixelValue = h0 + h1 * pixelValue;
				winR.at<uchar>(i - x1 + L_Size, j - y1 + L_Size) = static_cast<uchar>(pixelValue);;

				// 构建误差方程
				double gxDst = 0.5 * (right_copy.at<uchar>(I + 1, J) - right_copy.at<uchar>(I - 1, J));
				double gyDst = 0.5 * (right_copy.at<uchar>(I, J + 1) - right_copy.at<uchar>(I, J - 1));
				A(num, 0) = 1;
				A(num, 1) = pixelValue;
				A(num, 2) = gxDst;
				A(num, 3) = m * gxDst;
				A(num, 4) = n * gxDst;
				A(num, 5) = gyDst;
				A(num, 6) = m * gyDst;
				A(num, 7) = n * gyDst;

				L(num, 0) = left_copy.at<uchar>(i, j) - pixelValue;

				// 计算最佳匹配点位
				double gxSrc = 0.5 * (left_copy.at<uchar>(i + 1, j) - left_copy.at<uchar>(i - 1, j));
				double gySrc = 0.5 * (left_copy.at<uchar>(i, j + 1) - left_copy.at<uchar>(i, j - 1));

				xNumerator += i * gxSrc * gxSrc;
				xDenominator += gxSrc * gxSrc;
				yNumerator += j * gySrc * gySrc;
				yDenominator += gySrc * gySrc;

				num++;
			}
		if (num < 8) // 无法求解法方程
			return;

		tempCorIdx = Get_NCC(winL, winR);

		// 计算变形参数
		x = (A.transpose() * A).inverse() * (A.transpose() * L);
		double a0_old = a0;
		double a1_old = a1;
		double a2_old = a2;
		double b0_old = b0;
		double b1_old = b1;
		double b2_old = b2;
		double h0_old = h0;
		double h1_old = h1;

		a0 = a0_old + x(2, 0) + a0_old * x(3, 0) + b0_old * x(4, 0);
		a1 = a1_old + a1_old * x(3, 0) + b1_old * x(4, 0);
		a2 = a2_old + a2_old * x(3, 0) + b2_old * x(4, 0);
		b0 = b0_old + x(5, 0) + a0_old * x(6, 0) + b0_old * x(7, 0);
		b1 = b1_old + a1_old * x(6, 0) + b1_old * x(7, 0);
		b2 = b2_old + a2_old * x(6, 0) + b2_old * x(7, 0);
		h0 = h0_old + x(0, 0) + h0_old * x(1, 0);
		h1 = h1_old + h1_old * x(1, 0);

		// 计算最佳匹配点位
		double xt = xNumerator / xDenominator;
		double yt = yNumerator / yDenominator;

		xs = a0 + a1 * xt + a2 * yt;
		ys = b0 + b1 * xt + b2 * yt;

		if (tempCorIdx > bestCorIdx)
		{
			bestPt.x = ys;
			bestPt.y = xs;
			bestCorIdx = tempCorIdx;
		}

		if (bestCorIdx > L_Threshold)
		{
			match.right.x = bestPt.x;
			match.right.y = bestPt.y;
			match.r = bestCorIdx;
			return;
		}
	}
	match.right.x = bestPt.x;
	match.right.y = bestPt.y;
	match.r = bestCorIdx;
	return;
}

Mat Calculation::DrawRst()
{
	if (MatchImg.empty()) {
		Mat mat1 = cv::Mat::zeros(3, 3, CV_8UC1);
		return mat1;
	}
	cvtColor(MatchImg, MatchImg, COLOR_GRAY2BGR);
	for (Match_Pts& pts : Match_Points) {
		if (pts.r < 0.85) {
			continue;
		}
		Point pl = pts.left;
		Point pr = pts.right;
		pr.x += LeftImg.cols;//偏移
		Scalar randomColour(
			std::rand() % 256, // 随机红色分量  
			std::rand() % 256, // 随机绿色分量  
			std::rand() % 256  // 随机蓝色分量
		);
		circle(MatchImg, pl, 5, randomColour, 1);
		circle(MatchImg, pr, 5, randomColour, 1);
		line(MatchImg, pl, pr, randomColour, 2, LINE_AA);
	}
	return MatchImg;
}



Calculation::Calculation()
{
}

Calculation::~Calculation()
{
}
