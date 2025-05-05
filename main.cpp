#include "photogrammetry.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    photogrammetry w;
    w.show();
    return a.exec();
}
