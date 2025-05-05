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

    // ��ȡ����ĳߴ��ͨ����
    int width = mat.cols;
    int height = mat.rows;
    int bytesPerLine = mat.step;
    int channels = mat.channels();

    // ���ݲ�ͬ����ɫ�ռ����ת��
    QImage img;
    if (mat.type() == CV_8UC1) { // �Ҷ�ͼ��
        img = QImage(mat.data, width, height, bytesPerLine, QImage::Format_Grayscale8);
    }
    else if (mat.type() == CV_8UC3) { // BGRͼ��
        // OpenCV Ĭ��ʹ�� BGR ��ɫ�ռ䣬�� Qt ʹ�� RGB��������Ҫת����ɫ�ռ�
        // OpenCV Ĭ��ʹ�� BGR ��ɫ�ռ䣬�� Qt ʹ�� RGB��������Ҫת����ɫ�ռ�
        Mat rgb;
        cvtColor(mat, rgb, cv::COLOR_BGR2RGB);
        img = QImage(rgb.data, width, height, width * channels, QImage::Format_RGB888).copy();
    }
    else if (mat.type() == CV_8UC4) { // BGRAͼ��
        Mat rgba;
        cvtColor(mat, rgba, cv::COLOR_BGRA2RGBA);
        img = QImage(rgba.data, width, height, bytesPerLine, QImage::Format_RGBA8888);
    }
    else {
        // �����������͵�ͼ������ֻ���򵥵Ĵ���ʵ��Ӧ���п�����Ҫ�����ӵ�ת���߼�
        return QImage();
    }

    // ����ת����� QImage
    return img; // ��֤���ص������������ԭʼ���ݱ��ͷź��µ�����
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
        ui->result_image->setPixmap(pixmap.scaled(ui->result_image->size(), Qt::KeepAspectRatio));//��ʾ������
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

    // д�� Match_Pts ����
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