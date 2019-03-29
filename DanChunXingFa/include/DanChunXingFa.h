/*******************************************
 * 单纯形法C++实现.于2019.3.21开始
 * 目前已完成:
 * 0.单纯型表算法
 * 1.标准型问题的计算
 * 还需要完成:
 * +对于添加人工变量的问题的解决?
 * +对于b向量中存在小于0的量的预处理?
 * +对于min目标函数问题的求解?
 * +大M法,两个状态法的实现(即去除人工变量的部分)
 * +更多兼容性问题的思考?
 * +GUI用户界面?
 * +...(想到再补充
 *
 * --liangzia!
 *
 * ****************************************/
#ifndef DANCHUNXINGFA_DANCHUNXINGFA_H
#define DANCHUNXINGFA_DANCHUNXINGFA_H

#define MAXNUM 1000

#include "Eigen/Eigen"
#include <iostream>
#include <string>
using namespace Eigen;
using namespace std;


class DanChunXingFa{
public:
    DanChunXingFa()  //构造函数
    {
        this->ebuxinong=1e-3;
    }
    ~DanChunXingFa() //析构函数
    {
        ;
    }

    //用于求解标准型问题的方法,关于何为标准型请见定义部分
    VectorXd solveBiaoZhunXing(RowVectorXd &C, MatrixXd &allA, VectorXd &b);

    //通过交互输入数据来对模型问题进行求解.
    VectorXd solveWithString(const string &problem);

private:

    template<typename T>
    T getmin(T x,T y)
    {
        return x<=y?x:y;
    }

    template<typename T>
    T getmax(T x,T y)
    {
        return x>=y?x:y;
    }

    template<typename T1,typename T2>
    T2 get_max_and_max_index(T1 &storeMax,const T1 &compareValue)
    {
        if(storeMax<compareValue)
        {
            storeMax=compareValue;
            return 1;
        }
        else
        {
            return 0;
        }

    }

    VectorXi getXuLieImprove(long from,long length);

    void xuanZhuanJuZhenFa(MatrixXd &allA,VectorXd &b,VectorXi &XB);
    //求得当前最优解对应的目标最大值
    double getMaxNum(const VectorXi &XB, const RowVectorXd &C, const VectorXd &b);

    VectorXd getBestSolution(const VectorXd &b,const VectorXi &XB,int num_of_juecebianliang);

    int isTheBestResult(const RowVectorXd &C,const VectorXi &XB, const MatrixXd &allA);

    int getXuanChuIndex(const VectorXd &b,const MatrixXd &allA, const int &flag);

    void replace_XB_vector(const int &flag,const int &XuanChuIndex,VectorXi &XBB);

    VectorXd using_danchunxing_Table(const RowVectorXd &C,MatrixXd allA,VectorXi XB,VectorXd &b);

    void getDataFromString(const string &problem,RowVectorXd &C,MatrixXd &dayuA,MatrixXd &dengyuA,MatrixXd &xiaoyuA,
            VectorXd &dayub,VectorXd &dengyub,VectorXd &xiaoyub);

    //定义一个用来减小数值误差的变量
    double ebuxinong;
};

/***********************************************分界线******************************************************************/


/*******************************************
 * 一个标准型通常满足以下定义:
 * 1.最大化问题
 * 2.所有的约束均为等式
 * 3.b向量的值恒大于(特殊的,在退化时,等于)零.
 * *****************************************
 *
 * */
VectorXd DanChunXingFa::solveBiaoZhunXing(RowVectorXd &C, MatrixXd &allA, VectorXd &b)
{
    //使用最后几列作为基变量.即使不符合要求.
    VectorXi XB=getXuLieImprove(allA.cols()-allA.rows(),allA.cols());

    VectorXd x=using_danchunxing_Table(C,allA,XB,b);
    return x;
}

void DanChunXingFa::replace_XB_vector(const int &flag,const int &XuanChuIndex,VectorXi &XBB)
{
    XBB(XuanChuIndex)=flag;
}


/***********************************************************************************************************
     * 标准型单纯型表思路:
     * 1.对于得到的系数矩阵allA,使用初等行变换将其最后面的方阵化为标准型,此处涉及到一个化为标准型的函数XuanZhuanJuZhenFa
     * 2.计算出目标最大值
     * 3.得到最优解
     * 4.判断是否当前解为最优解,如果是,返回一个标志,反之,返回需要被选入的非基变量的索引
     * 5.如果4.中表明已经得到了最优解,退出程序,返回结果;否则进行第六步
     * 6.利用被选入的非基变量的索引计算出thita向量,并从中找到最小的值对应的索引作为基变量中的选出变量
     * 7.根据换入变量和换出变量修正向量XB及其对应的系数CB.
     * 8.旋转矩阵法,跳转至第二步.
     * *****************************************************************************************************
     * */
VectorXd DanChunXingFa::using_danchunxing_Table(const RowVectorXd &C,MatrixXd allA,VectorXi XB,VectorXd &b)
{
    xuanZhuanJuZhenFa(allA,b,XB);
    int i=0;
    while(true)
    {
        cout<<"===================正在进行第"<<i<<"次迭代===================="<<endl;
        cout<<"系数矩阵为:\n"<<allA<<"\nb向量为:\n"<<b<<"\nXB向量为:\n"<<XB<<endl;
        i++;
        if(i>100)
        {
            cout<<"程序执行此处过多!\n已退出..."<<endl;
            return b.setZero();
        }
        double max_value;
        max_value=DanChunXingFa::getMaxNum(XB,C,b);
        VectorXd bestsolution=getBestSolution(b,XB,(int)C.size());
        int flag=isTheBestResult(C,XB,allA);

        if(flag<0)
        {
            std::cout<<"已经找到了最优解,现在准备把它输出出来..."<<std::endl;
            std::cout<<"目标最优解为:\n"<<bestsolution<<std::endl;
            std::cout<<"最优解对应的目标函数最优值为"<<max_value<<std::endl;
            return bestsolution;
        }
        else
        {
            int XuanChuIndex=getXuanChuIndex(b,allA,flag);
            if(XuanChuIndex<0)
            {
                std::cout<<"程序出错,出错位置:getXuanChuIndex函数."<<std::endl;
            }
            //更新XBB,将换入变量替换掉换出变量.
            replace_XB_vector(flag,XuanChuIndex,XB);

            //不需要更新,因为没有将CB进行显式的定义
            //replace_XB_and_CB_vector(C,flag,XuanChuIndex,XBB,CB);
            xuanZhuanJuZhenFa(allA,b,XB);
            #ifdef debug
            cout<<"=======================debug================"<<endl;
            cout<<"换入变量为:x"<<flag<<endl;
            cout<<"换出变量的索引为:"<<XuanChuIndex<<endl;
            cout<<"XB向量的结果为:"<<XB<<endl;
            cout<<allA<<endl;
            cout<<"========================enddebug============"<<endl;
            #endif
        }
    }

}

int DanChunXingFa::getXuanChuIndex(const VectorXd &b,const MatrixXd &allA, const int &flag)
{
    VectorXd a_vector=allA.col(flag);
    VectorXd thita(b.size());
    thita.setZero();
    thita=thita*1000000;
    if(allA.any()>1000000)
    {
        std::cout<<"warning:现在正在使用的1000000作为最大量,如果您想系数有大于该数值的,则需要进行修改."<<std::endl;
        std::cout<<"修改地址为:DanChunXingFa.cpp -> class: DanChunXingFa -> function: getXuanChuIndex"<<std::endl;
    }
    double min=1000000;
    int index=-1;
    for(int i=0;i<b.size();i++)
    {
        if (a_vector(i)<=0)
        {
            continue;
        }
        thita(i)=b(i)/a_vector(i);
        if(min>thita(i))
        {
            min=thita(i);
            index=i;
        }
    }
    return index;
}

int DanChunXingFa::isTheBestResult(const RowVectorXd &C,const VectorXi &XB, const MatrixXd &allA)
{
    VectorXd CB(XB.size());
    for(int i=0;i<XB.size();i++)
    {
        CB(i)=C(XB(i));
    }
    RowVectorXd list_using_to_select_Which_in(allA.cols());
    int flag=-1;
    double value_select_in=0;
    for(int i=0;i<allA.cols();i++)
    {
        list_using_to_select_Which_in(i)=C(i)-CB.dot(allA.col(i));
        if(list_using_to_select_Which_in(i)>0)
        {
            int isDaYu=0;
            isDaYu=DanChunXingFa::get_max_and_max_index<double,int>(value_select_in,list_using_to_select_Which_in(i));
            if(isDaYu)
            {
                flag=i;
            }
        }
    }
    return flag;
}

double DanChunXingFa::getMaxNum(const VectorXi &XB, const RowVectorXd &C, const VectorXd &b)
{
    VectorXd CB(XB.size());
    for(int i=0;i<XB.size();i++)
    {
        CB(i)=C(XB(i));
    }
    double result;
    result=CB.dot(b);
    return result;
}

VectorXd DanChunXingFa::getBestSolution(const VectorXd &b,const VectorXi &XB,int num_of_juecebianliang)
{
    VectorXd result(num_of_juecebianliang);
    result.setZero();
    for(int i=0;i<XB.size();i++)
    {
        result(XB(i))=b(i);
    }
    return result;
}

void DanChunXingFa::xuanZhuanJuZhenFa(MatrixXd &allA,VectorXd &b,VectorXi &XB)
{
    for(int i=0;i<XB.size();i++)
    {
        //首先进行归一化
        if(allA(i,XB(i))!=1)
        {
            double temp1= allA(i,XB(i));
            for(int ii=0;ii<allA.cols();ii++)
            {
                allA(i,ii)/=temp1;
            }
            b(i)/=temp1;
        }
        //之后将该列里面非该行的所有元素化为0,实现目的.
        for(int j=0;j<XB.size();j++)
        {
            double temp2= allA(j,XB(i));
            if(j==i)
                continue;
            for(int jj=0;jj<allA.cols();jj++)
            {
                allA(j,jj)-=allA(i,jj)*temp2;
            }
            b(j)-=b(i)*temp2;
        }
    }
}


VectorXi DanChunXingFa::getXuLieImprove(long from,long length)
{
    long tanqi=from;
    VectorXi XB(length-from);
    for(int i=0;i<length-from;i++)
    {
        XB(i)=tanqi;
        tanqi++;
    }
    return XB;
}

/***************************************************************************************************************************
 * solve方法[已经废弃]:求解一般线性规划问题,在交互页面手动输入数据.(适用于小型求解问题)
 * 1.询问为最小优化问题还是最大优化问题
 * 2.请求输入C向量的数字
 * 3.询问小于的约束是否存在,如果是,请求输入
 * 4.询问大于的约束是否存在,如果是,请求输入
 * 5.询问等于的约束是否存在,如果是,请求输入
 * 6.结束询问.生成问题模型并打印到终端,询问是否正确,如果错误,全部清零,跳转到第一步
 * 7.将得到的小于等于大于约束转化为等式约束,更新C,并在等于与大于约束后面添加人工变量
 * 8.分类讨论:
 *      1)如果使用两个状态方法,则:对人工变量生成一个最大化问题目标函数,并对该问题使用单纯型表进行求解,直至求得最优解为止;之后进入第二个状态,
 *      对排除了非基变量的问题调用标准型法进行求解,直至问题解决;
 *      2)如果使用大M法,则在allA,b,C向量中找寻最大的数值,并将该数值的1000倍选取为大M,更新C向量,之后对形成的标准型采用单纯形法进行求解.
 * 8.将目标结果显示,跳转到第一步(是的,这是一个死循环)
 ***************************************************************************************************************************
 * 异常处理:
 * 1.如果人工变量无法被换出?
 * *************************************************************************************************************************/

//将字符串的数据读取出来.
void getDataFromString(const string &problem,RowVectorXd &C,MatrixXd &dayuA,MatrixXd &dengyuA,MatrixXd &xiaoyuA,
                       VectorXd &dayub,VectorXd &dengyub,VectorXd &xiaoyub)
{
    ;
}


/***********************************************************************************************************************
 * solveWithString:求解一般线性规划问题的最好办法.
 * 1.从字符串中提取出来所需要的所有信息,包括C,dayuA,dayub,xiaoyuA,xiaoyub,dengyuA,dengyub
 * 2.对于小于和等于的问题,分别添加决策变量,修正dayuA,xiaoyuA,C(添加0)
 * 3.对于大于和等于的情况,在其末尾添加人工变量,并跳转到第四步:
 * 4.分类讨论:
 *      1)如果使用两个状态方法,则:对人工变量生成一个最大化问题目标函数,并对该问题使用单纯型表进行求解,直至求得最优解为止;之后进入第二个状态,
 *      对排除了非基变量的问题调用标准型法进行求解,直至问题解决;
 *      2)如果使用大M法,则在allA,b,C向量中找寻最大的数值,并将该数值的1000倍选取为大M,更新C向量,之后对形成的标准型采用单纯形法进行求解.
 * 5.将结果导出,结束.
 ***********************************************************************************************************************/
VectorXd DanChunXingFa::solveWithString(const string &problem)
{
    RowVectorXd C;
    MatrixXd dayuA,xiaoyuA,dengyuA;
    VectorXd dayub,xiaoyub,dengyub;

    getDataFromString(problem,C,dayuA,dengyuA,xiaoyuA,dayub,dengyub,xiaoyub);

//-------------------------------------------------代码写到这里啦!-------------------------------------------------------------------------


}






#endif //DANCHUNXINGFA_DANCHUNXINGFA_H
