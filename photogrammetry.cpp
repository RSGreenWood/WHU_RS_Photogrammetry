#include "photogrammetry.h"
#include "ui_photogrammetry.h"
#include <QVBoxLayout> 

photogrammetry::photogrammetry(QWidget *parent)
    : QWidget(parent), ui(new Ui::photogrammetryClass)
{
    ui->setupUi(this);
    //Mat left_img, right_img;

    //载入学院图标
    img_depart = new QLabel(this);
//GetDelta = new QLineEdit(this);
    QString imgPath=":/photogrammetry/department";
    QPixmap pixmap(imgPath);
    QSize originalSize = pixmap.size();
    QSize NewSize = 0.1 * originalSize;
    QPixmap scaledPixmap = pixmap.scaled(NewSize, Qt::KeepAspectRatio, Qt::SmoothTransformation);
    if (!pixmap.isNull()) {
        img_depart->setPixmap(scaledPixmap);
    }
    else {
        img_depart->setText("error");
    }

    // 确保 QLineEdit 是可见的和启用的
    ui->GetDelta->setEnabled(true);
    ui->GetDelta->setVisible(true);

    // 打印 QLineEdit 的可见性和启用状态以进行调试
    qDebug() << "QLineEdit visible:" << ui->GetDelta->isVisible();
    qDebug() << "QLineEdit enabled:" << ui->GetDelta->isEnabled();
    //设置窗口布局
    
    //无需创建按钮响应，qt会自动连接
    //connect(ui->Load_Left_Button, SIGNAL(clicked()), this, SLOT(on_Load_Left_Button_clicked()));
    //connect(ui->Load_Right_Button, SIGNAL(clicked()), this, SLOT(on_Load_Right_Button_clicked()));
}

photogrammetry::~photogrammetry()
{
    delete ui;
    delete img_depart;
}


Mat photogrammetry::QImage2Mat(const QImage& image)
{
    // 检查输入的 QImage 是否是有效的
    if (image.isNull())
    {
        return Mat();
    }

    // 根据 QImage 的格式选择合适的 Mat 类型
    switch (image.format())
    {
    case QImage::Format_RGB32: {
        Mat mat(image.height(), image.width(), CV_8UC4, const_cast<uchar*>(image.bits()), image.bytesPerLine());

        // 创建一个深拷贝，避免 QImage 被销毁后指向无效内存
        Mat matCopy;
        mat.copyTo(matCopy);

        // 转换颜色空间，从 RGB 到 BGR
        //cvtColor(matCopy, matCopy, COLOR_RGBA2BGR);

        return matCopy;
    }
    case QImage::Format_ARGB32:
    {
        Mat mat(image.height(), image.width(), CV_8UC4, const_cast<uchar*>(image.bits()), image.bytesPerLine());
        Mat matCopy;
        mat.copyTo(matCopy); 
        cvtColor(matCopy, matCopy, COLOR_RGBA2BGR); // 转换颜色空间
        return matCopy;
    }
    case QImage::Format_RGB888:
    {
        Mat mat(image.height(), image.width(), CV_8UC3, const_cast<uchar*>(image.bits()), image.bytesPerLine());
        Mat matCopy;
        mat.copyTo(matCopy); 
        cvtColor(matCopy, matCopy, cv::COLOR_BGRA2BGR);
        return matCopy;
    }
    case QImage::Format_Grayscale8:
    {
        Mat mat(image.height(), image.width(), CV_8UC1, const_cast<uchar*>(image.bits()), image.bytesPerLine());
        Mat matCopy;
        mat.copyTo(matCopy); 

        return matCopy;
    }
    default:
        qWarning("Unsupported QImage format");
        return Mat();
    }
}

QPixmap photogrammetry::Mat2QPix(Mat &ori)
{
    Mat rgb;
    cvtColor(ori, rgb, cv::COLOR_BGR2RGB); // 转换为RGB格式，因为QImage默认以RGB格式处理图像  
    QImage image(rgb.data, rgb.cols, rgb.rows, rgb.step, QImage::Format_RGB888);
    return QPixmap::fromImage(image);
}


void photogrammetry::on_Load_Left_Button_clicked() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Image"), "", tr("Images (*.png *.xpm *.jpg)"));
    if (!fileName.isEmpty()) {
        QImage qimg(fileName);
        if (qimg.isNull())
        {
            QMessageBox::information(this, "Image Viewer", "Cannot load " + fileName);
            return;
        }
        QPixmap pixmap=QPixmap::fromImage(qimg);
        if (!pixmap.isNull()) {
            ui->LeftImg->setPixmap(pixmap.scaled(ui->LeftImg->size(), Qt::KeepAspectRatio));//显示到窗口
            this->left_img = QImage2Mat(qimg);
            return;
        }
        else {
            QMessageBox::warning(this, tr("Image Viewer"), tr("Cannot load %1.").arg(fileName));
            return ;
        }
    }
    else {
        QMessageBox::warning(this, tr("File does not exit"),tr("!!!"));
        return;
    }
}

void photogrammetry::on_Load_Right_Button_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Image"), "", tr("Images (*.png *.xpm *.jpg)"));
    if (!fileName.isEmpty()) {
        QPixmap pixmap(fileName);
        QImage qimg(fileName);
        if (!pixmap.isNull()) {
            ui->RightImg->setPixmap(pixmap.scaled(ui->RightImg->size(), Qt::KeepAspectRatio));//显示到窗口
            this->right_img = QImage2Mat(qimg);
        }
        else {
            QMessageBox::warning(this, tr("Image Viewer"), tr("Cannot load %1.").arg(fileName));
            return;
        }
    }
    else {
        QMessageBox::warning(this, tr("File does not exit"), tr("!!!"));
        return;
    }
}

void photogrammetry::on_Caculation_clicked()
{
    qDebug() << "按钮触发";
    bool ok;
    QString qinput = ui->GetDelta->text();
    int num = qinput.toInt(&ok);
    Result_Show result(num,left_img,right_img,this);
    result.exec();
}