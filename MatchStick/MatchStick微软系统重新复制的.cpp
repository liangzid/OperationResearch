/**********************************************************
����С��Ϸ������ϷĿǰ��ʵ�����¹��ܣ�
1������ȷ����������
2������ȷ���з��ҷ�������ٿ���ѡ���ٻ���
���޷�ʵ�ֵ��ǣ�
1���޷���ȷ���������ȫʤ����ʱ���������ʽ����
2��Ŀǰֻ�ܼ�����ߵ�һ��


liangzid
2019��4��29��
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
    int nownum=0;//�Ѿ�����ȡ�Ļ�������
    vector<int> needGet;
    string yesorno;
    cout<<"��ӭ����������Ϸ,��������Ĭ�ϵ�����[������: y]�����Զ�������[������: n]��?"<<endl;
    cin>>yesorno;
    if(yesorno=="y"||yesorno=="yes")
    {
        //����һ�����ֿ�ʼ����
        allnum=30;
        iMaxnum=3;iMinnum=1;
        heMaxnum=3;heMinnum=1;
    }
    else if(yesorno=="n"||yesorno=="no")
    {
        cout<<"�������ܹ��Ļ�������:";
        cin>>allnum;
        cout<<"\n����������������ܳ�ȡ�Ļ�������:";
        cin>>iMaxnum;
        cout<<"\n�����������������ܳ�ȡ�Ļ�������:";
        cin>>iMinnum;
        cout<<"\n��������������ܳ�ȡ�Ļ�������:";
        cin>>heMaxnum;
        cout<<"\n���������������ܳ�ȡ�Ļ�������:";
        cin>>heMinnum;
        cout<<"\n";
    }

    //�����Щ������Ϣ
    cout<<"**********��Ϸ������Ϣ���**********"<<endl;
    cout<<"Ŀǰ�ܹ�ӵ�еĻ�������:"<<allnum<<endl;
    cout<<"����������Գ�ȡ"<<iMaxnum<<"��,���ٿ��Գ�ȡ"<<iMinnum<<"��.\n";
    cout<<"�����������Գ�ȡ"<<heMaxnum<<"��,���ٿ��Գ�ȡ"<<heMinnum<<"��.\n";
    cout<<"**********��Ϸ�������**********"<<endl;
    cout<<"˭�õ����һ������,˭������.\n";
    cout<<"**********��ô,��Ϸ��ʼ!----------\n"<<endl;

    //��ȡ����Ҫ�õ��Ļ�������
    needGet=NumberOfMatchsticksNeedToYouMustNeedToGet(allnum,iMaxnum,iMinnum,heMaxnum,heMinnum);
    
    //����ÿ�ε�����Ҫ���ߵĳ�ȡ��Ŀ�ͶԷ�ѡ����г�ȡ����Ŀ
    int iNeedDo,heHasDone;
    //��һ���Ǽ������ȡ,�����׼����ȡ����һ������ִ�ĵ㴦
    iNeedDo=needGet[needGet.size()-1];
    if (iNeedDo==-1)
    {
        //��ʼ���Ϲ��״̬...
        cout<<"�ò��ֻ�δ���,�������û�а�����Ӯ�������µľ��߻�δ���.\n����������......";
    }
    else//��ʼ��ʤ
    {
        nownum=iNeedDo;//������Ѿ�����˵�һ��ץȡ
        cout<<"�����ѡ��ץȡ"<<iNeedDo<<"��,Ŀǰ�ܹ���ץȥ��"<<nownum<<"��,��ʣ"<<(allnum-nownum)<<
        "��,��ô.\n";
        for(int i=1;nownum<=allnum;i++)
        {
            cout<<"�����׼����ȡ����?"<<endl;
            //�˴���Ҫ�жϽƻ��������Ƿ������������ȷ������
            while(1)
            {
                cin>>heHasDone;
                if(heHasDone>=heMinnum&&heHasDone<=heMaxnum&&heHasDone<=(allnum-nownum))
                {
                    break;
                }
                else
                {
                    cout<<"+-+-+-+-+-+-+-����!������벻���Ϲ���!\n�鿴һ����������Ƿ����"<<heMinnum
                    <<"��"<<heMaxnum<<"֮��,�������ȡ����Ŀ�Ƿ�����˻�����ʣ����.\n";
                    cout<<"�����Ϻ�,����������:\n";
                }
            }
            nownum+=heHasDone;
            cout<<"��ѡ���ȡ��"<<heHasDone<<"��,�ܲ���,Ŀǰ�ܹ���ȡ����"<<nownum<<"��,��ʣ"
            <<(allnum-nownum)<<"��.\n";
            if(nownum==allnum)
            {
                cout<<"================================\n";
                cout<<"�޴������డ,��ʧ���� 0_0"<<endl;
                cout<<"================================\n";
                break;
            }
            //�������Ҫ���о���,��Ȼ��Ҫ������ȱ,׼ʱ������һ����λ
            iNeedDo=needGet[needGet.size()-1-i]-nownum;
            nownum+=iNeedDo;
            cout<<"----------.\n�����������ȡ"<<iNeedDo<<"��,Ŀǰ�ܹ�����ȥ��"<<nownum<<"��,��ʣ"
            <<(allnum-nownum)<<"��,��ô\n";
            if(nownum==allnum)
            {
                cout<<"��Ȼ����һ�е�������Զ���ᱻ��ӡ����"<<endl;
            }

        }
    }
    selfPause();
    return 0;
}

//�����㷨--�������Ҫ��ȡ�Ļ�������
vector<int> NumberOfMatchsticksNeedToYouMustNeedToGet(int numOfALlMatchSticks,int IMaxGet,
int IMinGet, int heMaxGet, int heMinGet)
{
    vector<int> numNeedGet;
    if(numOfALlMatchSticks<=IMaxGet)
    {
        //���һ�ξ���ץȡ���Ļ�,��ô��ץ��ֻʣ�����һ��,�����Է�
        numNeedGet.push_back(numOfALlMatchSticks-heMinGet);
        return numNeedGet;
    }
    else if(numOfALlMatchSticks==(min(IMaxGet+heMinGet,IMinGet+heMaxGet)+IMinGet))
	{
		cout<<"������ܾ����б�����Ϸ���뻻���������֣���Ϸ������\n";
		numNeedGet.push_back(-1);
        return numNeedGet;
	 } 
    
    else
    {
        //һ�β���ץȡ��,��ʱ�Ҿͱ�����뷽�跨ץ����n-1��.һ����ץ����n-1����Ҫ����,���۶Է�ѡ��ץ���ٸ�(1->h),�Ҷ�
        //����ץ����n-1��.����Ҫ���������ҷ�����ץȡ�Ļ���������Χ����С�ڶԷ�����ץȡ�Ļ���������Χ.
        if((IMaxGet-IMinGet)>=(heMaxGet-heMinGet))
        {
            //���Ի�þ���ʤ��
            cout<<"\n------\n�������,��������Ի����ȫʤ��\n------\n"<<endl;
            //��ȡ����ʤ�������һϵ�е����ֵ
            int needGet=numOfALlMatchSticks-heMinGet;
            while(1)
            {
                numNeedGet.push_back(needGet);

                needGet=needGet-(min(IMaxGet+heMinGet,IMinGet+heMaxGet));

                if(needGet<=0)
                {
                    //��ý���
                    return numNeedGet;
                }
            }
        }
        else
        {
            //�����������������,����ֻ�ܸ���������ҵĴ���ˮƽ������ʤ��Ĳ�����.�����ֲ�����ʽ��,�������������
            //���ѱ�֤���յ�ʤ��
            cout<<"�������,�����Ϸ�Լ����û�а�����Ӯ,�������ʾ�����������Ϸ."<<endl;
            numNeedGet.push_back(-1);//-1�����޷���,��ʱ�������⴦��
            return numNeedGet;
        }
    }
}

void selfPause()
{
    string liangzi;
    cout<<"�������н���,���������ַ��˳���...... = =!\n";
    cin>>liangzi;
}
