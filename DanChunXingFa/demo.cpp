#include <iostream>
#include "Eigen/Eigen"
#include "DanChunXingFa.h"
using namespace std;
using namespace Eigen;

int main() {
    DanChunXingFa liangzi1;

    RowVectorXd C1(5);
    C1<<5,2,3,-1,1;
    MatrixXd allA1(2,5);
    allA1<<1,2,2,1,0,
            3,4,1,0,1;
    VectorXd b1(2);
    b1<<8,7;
    liangzi1.SolveBiaoZhunXing(C1,allA1,b1);

    RowVectorXd C2(5);
    C2<<1,2,0,0,0;
    MatrixXd allA2(3,5);
    allA2<<1,0,1,0,0,
           0,1,0,1,0,
           1,2,0,0,1;
    VectorXd b2(3);
    b2<<4,3,8;
    liangzi1.SolveBiaoZhunXing(C2,allA2,b2);

    RowVectorXd C3(5);
    C3<<40,45,24,0,0;
    MatrixXd allA3(2,5);
    allA3<<2,3,1,1,0,
           3,3,2,0,1;
    VectorXd b3(2);
    b3<<100,120;
    liangzi1.SolveBiaoZhunXing(C3,allA3,b3);

    return 0;
}