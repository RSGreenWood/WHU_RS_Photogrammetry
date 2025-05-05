#include "photogrammetry.h"
#include "ui_photogrammetry.h"
#include <QVBoxLayout> 

photogrammetry::photogrammetry(QWidget *parent)
    : QWidget(parent), ui(new Ui::photogrammetryClass)
{
    ui->setupUi(this);
    //Mat left_img, right_img;

    //����ѧԺͼ��
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

    // ȷ�� QLineEdit �ǿɼ��ĺ����õ�
    ui->GetDelta->setEnabled(true);
    ui->GetDelta->setVisible(true);

    // ��ӡ QLineEdit �Ŀɼ��Ժ�����״̬�Խ��е���
    qDebug() << "QLineEdit visible:" << ui->GetDelta->isVisible();
    qDebug() << "QLineEdit enabled:" << ui->GetDelta->isEnabled();
    //���ô��ڲ���
    
    //���贴����ť��Ӧ��qt���Զ�����
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
    // �������� QImage �Ƿ�����Ч��
    if (image.isNull())
    {
        return Mat();
    }

    // ���� QImage �ĸ�ʽѡ����ʵ� Mat ����
    switch (image.format())
    {
    case QImage::Format_RGB32: {
        Mat mat(image.height(), image.width(), CV_8UC4, const_cast<uchar*>(image.bits()), image.bytesPerLine());

        // ����һ����������� QImage �����ٺ�ָ����Ч�ڴ�
        Mat matCopy;
        mat.copyTo(matCopy);

        // ת����ɫ�ռ䣬�� RGB �� BGR
        //cvtColor(matCopy, matCopy, COLOR_RGBA2BGR);

        return matCopy;
    }
    case QImage::Format_ARGB32:
    {
        Mat mat(image.height(), image.width(), CV_8UC4, const_cast<uchar*>(image.bits()), image.bytesPerLine());
        Mat matCopy;
        mat.copyTo(matCopy); 
        cvtColor(matCopy, matCopy, COLOR_RGBA2BGR); // ת����ɫ�ռ�
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
    cvtColor(ori, rgb, cv::COLOR_BGR2RGB); // ת��ΪRGB��ʽ����ΪQImageĬ����RGB��ʽ����ͼ��  
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
            ui->LeftImg->setPixmap(pixmap.scaled(ui->LeftImg->size(), Qt::KeepAspectRatio));//��ʾ������
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
            ui->RightImg->setPixmap(pixmap.scaled(ui->RightImg->size(), Qt::KeepAspectRatio));//��ʾ������
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
    qDebug() << "��ť����";
    bool ok;
    QString qinput = ui->GetDelta->text();
    int num = qinput.toInt(&ok);
    Result_Show result(num,left_img,right_img,this);
    result.exec();
}