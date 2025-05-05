#include "result_show.h"
using namespace cv;
Result_Show::Result_Show(const int d, const Mat &l, const Mat &r,  QWidget* parent)
	: QDialog(parent), ui(new Ui::Result_ShowClass), left_img(l), right_img(r),delta_x(d)
{
    QPushButton* pushbutton = findChild<QPushButton*>("save_reslut");
    connect(pushbutton, &QPushButton::clicked, this, &Result_Show::on_save_points_clicked);
    pushbutton = findChild<QPushButton*>("save_points");
    connect(pushbutton, &QPushButton::clicked, this, &Result_Show::on_save_result_clicked);
	ui->setupUi(this);
	Calculation test;
	Mat temp_img;
	cvtColor(left_img, temp_img, COLOR_BGR2GRAY);
	test.LeftImg = temp_img.clone();
	cvtColor(right_img, temp_img, COLOR_BGR2GRAY);
	test.RightImg = temp_img.clone();
	test.Moravec(test.LeftImg);
    test.deltax = delta_x*-1;
	test.Match();
	Mat rst_mat = test.DrawRst();
    namedWindow("result", WINDOW_NORMAL);
    imshow("result", rst_mat);
    waitKey(-1);
    rst_q = mat2QImage(rst_mat);
    pts = test.Match_Points;

    write_points();
}

QImage Result_Show::mat2QImage(const cv::Mat& mat)
{
    if (mat.empty()) {
        return QImage();
    }

    // 获取矩阵的尺寸和通道数
    int width = mat.cols;
    int height = mat.rows;
    int bytesPerLine = mat.step;
    int channels = mat.channels();

    // 根据不同的颜色空间进行转换
    QImage img;
    if (mat.type() == CV_8UC1) { // 灰度图像
        img = QImage(mat.data, width, height, bytesPerLine, QImage::Format_Grayscale8);
    }
    else if (mat.type() == CV_8UC3) { // BGR图像
        // OpenCV 默认使用 BGR 颜色空间，而 Qt 使用 RGB，所以需要转换颜色空间
        // OpenCV 默认使用 BGR 颜色空间，而 Qt 使用 RGB，所以需要转换颜色空间
        Mat rgb;
        cvtColor(mat, rgb, cv::COLOR_BGR2RGB);
        img = QImage(rgb.data, width, height, width * channels, QImage::Format_RGB888).copy();
    }
    else if (mat.type() == CV_8UC4) { // BGRA图像
        Mat rgba;
        cvtColor(mat, rgba, cv::COLOR_BGRA2RGBA);
        img = QImage(rgba.data, width, height, bytesPerLine, QImage::Format_RGBA8888);
    }
    else {
        // 对于其他类型的图像，这里只做简单的处理，实际应用中可能需要更复杂的转换逻辑
        return QImage();
    }

    // 返回转换后的 QImage
    return img; // 保证返回的是深拷贝，避免原始数据被释放后导致的问题
}

Result_Show::~Result_Show()
{
	delete ui;
}


void Result_Show::on_save_points_clicked() {
    if (rst_q.isNull()) {
        return ;
    }
    rst_q.save(":/phottgrammetry/ResultImg", "PNG");
    return ;
}

void Result_Show::on_save_result_clicked()
{
    String filepath = ":/phottgrammetry/ResultPoints.txt";
    ofstream outFile(filepath);
    if (!outFile.is_open()) {
        return;
    }
    for (int i = 0; i < pts.size();i++) {
        outFile << "Points" << i << ":(" << pts[i].left.x << "," << pts[i].left.y << ")  ("
            << pts[i].right.x << "," << pts[i].right.y << ")" << pts[i].r << endl;
    }
    outFile.close();
    return ;
}

void Result_Show::write_points()
{
    QPixmap pixmap = QPixmap::fromImage(rst_q);
    if (!pixmap.isNull()) {
        ui->result_image->setPixmap(pixmap.scaled(ui->result_image->size(), Qt::KeepAspectRatio));//显示到窗口
    }
    else {
        ui->result_image->setText("error");
    }
    QString appDirPath = QCoreApplication::applicationDirPath();
    QString filepath = appDirPath + "/ResultPoints.txt";
    QFile file(filepath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        return;
    }

    QTextStream out(&file);

    // 写入 Match_Pts 数据
    for (const auto& pt : pts) {
        if (pt.r < 0.85) {
            continue;
        }
        out << "Left: (" << pt.left.x << ", " << pt.left.y << "), "
            << "Right: (" << pt.right.x << ", " << pt.right.y << "), "
            << "Correlation: " << pt.r << "\n";
    }

    file.close();
    if (rst_q.isNull()) {
        return;
    }
    QString imgpath = appDirPath + "/ResultImg.png";
    rst_q.save(imgpath);
}