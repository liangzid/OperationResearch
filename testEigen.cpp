#include<iostream>
#include "Eigen/Eigen"
using namespace std;
using namespace Eigen;

int main()
{
    MatrixXd testMatriXd;
    cout<<"=====\n"<<testMatriXd<<endl;
    MatrixXd test2(0,0);
    cout<<"----\n"<<test2<<endl;
    cout<<testMatriXd.cols()<<endl;
    cout<<test2.cols()<<test2.rows()<<endl;

    //test2(3,3);
    
    cout<<"this is a test.\n"<<endl;
    return 0;
}
