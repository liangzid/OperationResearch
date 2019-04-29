/**********************************************************
火柴棒小游戏！该游戏目前可实现以下功能：
1）任意确定火柴棒个数
2）任意确定敌方我方最多最少可以选多少火柴棒
仍无法实现的是：
1）无法再确定计算机完全胜利的时候进行启发式下棋
2）目前只能计算机走第一步


liangzid
2019年4月29日
**********************************************************/
# include<iostream>
# include<string>
# include<vector>

using namespace std;

vector<int> NumberOfMatchsticksNeedToYouMustNeedToGet(int numOfALlMatchSticks,int IMaxGet,
int IMinGet, int heMaxGet, int heMinGet);
void selfPause();

int main(int argc,char** argv)
{
    int allnum,iMaxnum,heMaxnum,iMinnum,heMinnum;
    int nownum=0;//已经被抽取的火柴棒个数
    vector<int> needGet;
    string yesorno;
    cout<<"欢迎来到火柴棒游戏,你是想用默认的数据[请输入: y]还是自定义数据[请输入: n]呢?"<<endl;
    cin>>yesorno;
    if(yesorno=="y"||yesorno=="yes")
    {
        //随便给一组数字开始尝试
        allnum=30;
        iMaxnum=3;iMinnum=1;
        heMaxnum=3;heMinnum=1;
    }
    else if(yesorno=="n"||yesorno=="no")
    {
        cout<<"请输入总共的火柴棒个数:";
        cin>>allnum;
        cout<<"\n请输入计算机最大所能抽取的火柴棒个数:";
        cin>>iMaxnum;
        cout<<"\n请输入计算机最少所能抽取的火柴棒个数:";
        cin>>iMinnum;
        cout<<"\n请输入你最大所能抽取的火柴棒个数:";
        cin>>heMaxnum;
        cout<<"\n请输入你最少所能抽取的火柴棒个数:";
        cin>>heMinnum;
        cout<<"\n";
    }

    //输出这些基本信息
    cout<<"**********游戏基本信息输出**********"<<endl;
    cout<<"目前总共拥有的火柴棒个数:"<<allnum<<endl;
    cout<<"计算机最多可以抽取"<<iMaxnum<<"根,最少可以抽取"<<iMinnum<<"根.\n";
    cout<<"而人类最多可以抽取"<<heMaxnum<<"根,最少可以抽取"<<heMinnum<<"根.\n";
    cout<<"**********游戏规则介绍**********"<<endl;
    cout<<"谁拿到最后一根火柴棒,谁就输了.\n";
    cout<<"**********那么,游戏开始!----------\n"<<endl;

    //获取必须要得到的火柴棒序列
    needGet=NumberOfMatchsticksNeedToYouMustNeedToGet(allnum,iMaxnum,iMinnum,heMaxnum,heMinnum);
    
    //定义每次迭代需要决策的抽取数目和对方选择进行抽取的数目
    int iNeedDo,heHasDone;
    //第一次是计算机抽取,计算机准备抽取到第一个必须抵达的点处
    iNeedDo=needGet[needGet.size()-1];
    if (iNeedDo==-1)
    {
        //开始随机瞎做状态...
        cout<<"该部分还未设计,即计算机没有把握稳赢的情形下的决策还未设计.\n待后来补充......";
    }
    else//开始获胜
    {
        nownum=iNeedDo;//计算机已经完成了第一次抓取
        cout<<"计算机选择抓取"<<iNeedDo<<"根,目前总共被抓去了"<<nownum<<"根,还剩"<<(allnum-nownum)<<
        "根,那么.\n";
        for(int i=1;nownum<=allnum;i++)
        {
            cout<<"这次你准备抽取几根?"<<endl;
            //此处需要判断狡猾的人类是否真的输入了正确的数字
            while(1)
            {
                cin>>heHasDone;
                if(heHasDone>=heMinnum&&heHasDone<=heMaxnum&&heHasDone<=(allnum-nownum))
                {
                    break;
                }
                else
                {
                    cout<<"+-+-+-+-+-+-+-错误!你的输入不符合规则!\n查看一下你的输入是否介于"<<heMinnum
                    <<"与"<<heMaxnum<<"之间,或者你抽取的数目是否大于了火柴棒的剩余数.\n";
                    cout<<"检查完毕后,请重新输入:\n";
                }
            }
            nownum+=heHasDone;
            cout<<"你选择抽取了"<<heHasDone<<"根,很不错,目前总共抽取掉了"<<nownum<<"根,还剩"
            <<(allnum-nownum)<<"根.\n";
            if(nownum==allnum)
            {
                cout<<"================================\n";
                cout<<"愚蠢的人类啊,你失败了 0_0"<<endl;
                cout<<"================================\n";
                break;
            }
            //计算机需要进行决策,显然需要补掉空缺,准时到达下一个点位
            iNeedDo=needGet[needGet.size()-1-i]-nownum;
            nownum+=iNeedDo;
            cout<<"----------.\n计算机决定抽取"<<iNeedDo<<"根,目前总共被抽去了"<<nownum<<"根,还剩"
            <<(allnum-nownum)<<"根,那么\n";
            if(nownum==allnum)
            {
                cout<<"虽然有这一行但是它永远不会被打印出来"<<endl;
            }

        }
    }
    selfPause();
    return 0;
}

//火柴棒算法--计算出需要抽取的火柴棒个数
vector<int> NumberOfMatchsticksNeedToYouMustNeedToGet(int numOfALlMatchSticks,int IMaxGet,
int IMinGet, int heMaxGet, int heMinGet)
{
    vector<int> numNeedGet;
    if(numOfALlMatchSticks<=IMaxGet)
    {
        //如果一次就能抓取尽的话,那么就抓的只剩下最后一根,留给对方
        numNeedGet.push_back(numOfALlMatchSticks-heMinGet);
        return numNeedGet;
    }
    else if(numOfALlMatchSticks==(min(IMaxGet+heMinGet,IMinGet+heMaxGet)+IMinGet))
	{
		cout<<"计算机拒绝进行本场游戏，请换用其他数字，游戏结束。\n";
		numNeedGet.push_back(-1);
        return numNeedGet;
	 } 
    
    else
    {
        //一次不能抓取尽,这时我就必须得想方设法抓到第n-1根.一定能抓到第n-1根的要求是,无论对方选择抓多少根(1->h),我都
        //必须抓到第n-1根.这种要求来自于我方所能抓取的火柴棒个数范围不能小于对方所能抓取的火柴棒个数范围.
        if((IMaxGet-IMinGet)>=(heMaxGet-heMinGet))
        {
            //可以获得绝对胜利
            cout<<"\n------\n评判完毕,计算机可以获得完全胜利\n------\n"<<endl;
            //获取距离胜利最近的一系列点的数值
            int needGet=numOfALlMatchSticks-heMinGet;
            while(1)
            {
                numNeedGet.push_back(needGet);

                needGet=needGet-(min(IMaxGet+heMinGet,IMinGet+heMaxGet));

                if(needGet<=0)
                {
                    //虚幻结束
                    return numNeedGet;
                }
            }
        }
        else
        {
            //如果不满足上述条件,我们只能根据人类玩家的聪明水平进行无胜算的博弈了.在这种博弈形式下,由于先天的劣势
            //很难保证最终的胜利
            cout<<"评判完毕,这个游戏对计算机没有把握能赢,计算机表示不想玩这个游戏."<<endl;
            numNeedGet.push_back(-1);//-1代表无返回,届时进行特殊处理
            return numNeedGet;
        }
    }
}

void selfPause()
{
    string liangzi;
    cout<<"程序运行结束,输入任意字符退出吧...... = =!\n";
    cin>>liangzi;
}
