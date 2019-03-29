import numpy as np

class DanChunXingFa:
    '''
    单纯形法更新:
    此次更新的内容包括:
    1) 将之前只能使用在标准型的要求去除,添加了二阶段法人工变量来解决更一般的线性规划问题;
    2) 添加了Bland准则以解决退化问题,但是该准则只会在出现了退化现象之后持续使用.目前测试良好,还没有发觉到死循环的例子.
    3) 改进了输入,以使其更加人性化,更加符合常理.通过对矩阵的类型特殊限定让计算的精度更加稳定.
    liangzid
    2019.3.19
    '''

    '''
    单纯形法的尝试编程。
    使用:将目标优化(max)函数的系数传入,将已经化为标准型的系数|增广矩阵传入,使用run()方法得到结果.
    比如:
    test=DanChunXingFa(target=np.array([2.,3.,0.,0.,0.]),constraint=np.array([[1.,2.,1.,0.,0.,8.],
                                                                              [4.,0.,0.,1.,0.,16.],
                                                                              [0.,4.,0.,0.,1.,12.]]))
    test.run()

    其结果为:
    /home/liangzi/anaconda3/bin/python3.7 /home/liangzi/code/mycode/OperationResearch/demo.py
    ====================================================================================================
    ====================================================================================================
    ====================================================================================================
    ====================================================================================================
    已经找到唯一的最优解
    最优的决策向量为:
    [[4.]
    [2.]
    [0.]
    [0.]
    [4.]],
    此时的目标函数最大值为:
    14.0

    liangzid
    2019.3.12

    '''

    def __init__(self,usingAimMax=True,target=np.zeros((0,0)),equ=np.zeros((0,20)),XiaoYu=np.zeros((0,20)),Dayu=np.zeros((0,20))):

        self.equ = equ
        '''
        if (equ.shape)[0] != 0:
        
            
            if (equ[:,-1] < 0).any:
                self.equ = 0. - equ
        else:
            self.equ = equ
            '''
        self.XiaoYu=XiaoYu
        self.DaYu=Dayu

        constraint,NumOfAdd,y=self._BiaoZhunHua(self.equ,self.DaYu,self.XiaoYu)

        self.M = np.max(constraint)*100
        print('在这个计算中大M的数值为{}'.format(self.M))
        # print(self.target.shape)
        self.target=np.zeros(y+NumOfAdd)
        self.usingAimMax = usingAimMax
        if self.usingAimMax:
            self.target[0:(y)] = target
        else:
            self.target[0:(y)] = 0. - target  # 将极小化问题转化为极大化问题

        #self.target=np.hstack((self.target,np.zeros((NumOfAdd,))))

        self.constraintA=constraint[:,:-1] # 约束条件的系数矩阵
        #print(self.constraintA)
        self.constraintb=constraint[:,-1] # 约束方程的b向量
        m,n=self.XiaoYu.shape
        # print(self.constraintb.shape)
        self.constraintb[m:-1]=(-1)*self.M
        self.x,self.y=self.constraintA.shape

        # 这些是在迭代时需要被频繁使用的变量
        self.bestSolution=np.zeros((self.y,1),dtype=np.double) # 存储历史最优解（由于贪心算法，所以历史最优解就是当前解）的向量

        self.basisVaries=np.zeros((self.x,1),dtype=np.double)  # 存储基变量的索引，从零开始
        self.basisNotVaries=np.zeros((self.y-self.x,1),dtype=np.double) # 存储非基变量的索引

        self.basisMatrix=np.eye(self.x,dtype=np.double) # 基向量形成的方阵
        self.bVector=np.zeros((self.x,1),dtype=np.double) # 每次迭代过程中的b向量
        self.basisNotMatrix=np.zeros((self.x,self.y-self.x),dtype=np.double) # 非基变量对应的系数的矩阵

        self.jueCeMax=np.zeros((1,self.y-self.x+1),dtype=np.double) # 使用非决策变量表示目标函数得到的系数，其中首位为常数项（即为目标最大值），
                                                    # 剩下的顺序按照非决策变量的顺序进行
        self.maxValueNow=self.jueCeMax[0]
        self.selectInput=-100 # 决策使用那个变量作为选入
        self.jueCeWhichOutputMatrix = np.zeros((self.x, self.y - self.x + 1),dtype=np.double)  # 存储考虑选出哪一个变量时进行的操作
        self.selectOutput=-100 # 决策使用那个变量作为选出

        self.yibusinong=1e-3 # 常量

    def _BiaoZhunHua(self,dengyu,dayu,xiaoyu):
        y=0
        dengyu_x,dengyu_y=dengyu.shape
        dengyu_y -= 1

        dayu_x,dayu_y=dayu.shape
        dayu_y -= 1

        xiaoyu_x,xiaoyu_y=xiaoyu.shape
        xiaoyu_y -= 1

        # if dengyu_y!=dayu_y|dengyu_y!=xiaoyu_y|dayu_y!=xiaoyu_y:
        #    raise Exception('错误!输入有误!\n您必须保证您输入的约束条件中决策变量的数目是相同的,缺少的量可以考虑用0代替.\n')
        #print(dengyu,dayu,xiaoyu,xiaoyu_y)
        if dengyu_x > 0:
            y = dengyu_y
        else:
            if dayu_x > 0:
                y = dayu_y
            else:
                if xiaoyu_x > 0:
                    # print('=========111111==================='+str(xiaoyu_y))
                    y = xiaoyu_y
                else:
                    # print(xiaoyu_y)
                    raise Exception("没有任何的有效输入,无法运行程序.\n")
        # print('=========111111==================='+str(y))
        # print(xiaoyu_y)


        if dengyu_y != y:
            dengyu = np.zeros((0, y + 1))
        if dayu_y!=y:
            dayu=np.zeros((0,y+1))
        if xiaoyu_y!=y:
            xiaoyu=np.zeros((0,y+1))

        # print(xiaoyu)
        # 查看一下是否存在b向量中小于0的项,如果存在的话,将其转化为大于0
        temp_vec=np.zeros((1,xiaoyu_y+1))
        i=0
        while 1:
            temp_rows,_=xiaoyu.shape
            if i==temp_rows:
                break

            if xiaoyu[i, -1] < 0:
                # print('是呀!')
                temp_vec=(-1)*xiaoyu[i,:]
                np.delete(xiaoyu,i,axis=0)
                np.insert(dayu,values=temp_vec,axis=0)
            else:
                i+=1


        i=0
        while 1:
            if i==dayu_x:
                break

            if dayu[i,-1]<0:
                # print('zhaodaonile')
                temp_vec=(-1)*dayu[i,:]
                np.delete(dayu,i,axis=0)
                np.insert(xiaoyu,values=temp_vec,axis=0)
            else:
                i+=1

        # 更新小于大于的参数
        engyu_x, dengyu_y = dengyu.shape
        dayu_x, dayu_y = dayu.shape
        xiaoyu_x, xiaoyu_y = xiaoyu.shape

        dengyu_y -= 1
        dayu_y -= 1
        xiaoyu_y -= 1


        # print(xiaoyu_x)
        constraint = np.zeros((dengyu_x + dayu_x + xiaoyu_x, y + dengyu_x + dayu_x * 2 + xiaoyu_x),
                              dtype=np.double)
        # print(constraint.shape)
        # print(xiaoyu_x)
        # 填补上这些数字
        constraint[0:xiaoyu_x, 0:(y)] = xiaoyu[:, 0:(xiaoyu_y)]
        constraint[xiaoyu_x:(xiaoyu_x + dayu_x), 0:(y)] = dayu[:, 0:-1]
        constraint[(xiaoyu_x + dayu_x):-1, 0:(y)] = dengyu[:, 0:-1]

        #对对大于的等式添加-1
        for i in range(dayu_x):
            constraint[xiaoyu_x+i,y+i]=-1
        #对小于的等式添加1
        for i in range(xiaoyu_x):
            constraint[0+i,y+dayu_x+i]=1

        #添加人工变量
        for i in range(dayu_x+dengyu_x):
            constraint[xiaoyu_x+i,y+dayu_x+xiaoyu_x+i]=1

        # print(xiaoyu.shape,dayu.shape,dengyu.shape)
        b=np.zeros((xiaoyu_x+dengyu_x+dayu_x,1))

        if dengyu_x!=0 :
            b=xiaoyu[:,-1]
            if dayu_x!=0:
                print('=====================')
                print(b.shape,dayu[:,-1])
                b=np.vstack((b,dayu[:,-1]))
                if dengyu_x!=0:
                    b=np.vstack((b,dengyu[:,-1]))
            else:
                if dengyu_x!=0:
                    b = np.vstack((b, dengyu[:, -1]))
        else:
            if dayu_x!=0:
                b=dayu[:,-1]
                if dengyu_x!=0:
                    b=np.vstack((b,dengyu[:,-1]))
            else:
                if dengyu_x!=0:
                    b =dengyu[:, -1]


        # b=np.zeros((dengyu_x+dayu_x+xiaoyu_x,1))
        # b=np.vstack((np.vstack((xiaoyu[:,-1],dayu[:,-1])),dengyu[:,-1]))

        NumOfAdd=dengyu_x+dayu_x*2+xiaoyu_x

        # print(b.shape)
        # print(constraint.shape)

        return np.hstack((constraint,b)),NumOfAdd,y


    def IterationInit(self):
        '''
        用来读取单位阵，从而获取初始解
        :return:
        '''
        if not self._equation((self.constraintA[:,(self.y-self.x)::]-np.eye(self.x)).all(),0): # 如果系数矩阵右边不是单位阵的话
            print('没有转化成标准形式，需要改进算法\n')
            return -1

        # 否则，定义出相应数量的状态变量，并解出其对应的数值
        self.basisVaries=np.arange(self.y)[(self.y-self.x)::]
        self.basisNotVaries=np.arange(self.y)[0:(self.y-self.x)]
        self.basisMatrix=np.eye(self.x)
        self.basisNotMatrix=self.constraintA[:,0:(self.y-self.x)]
        self.bVector=self.constraintb

        # 求出此时的解,即初始解:
        self.bestSolution=self._CalculateSolutionVector(self.basisVaries,self.basisNotVaries,
                                                       self.basisMatrix,self.basisNotMatrix,self.bVector)
        # 根据求得的解得到非决策变量表示决策变量的形式(非常简单)
        self.jueCeWhichOutputMatrix=self._CalculateBiaoShiMatrix(self.basisMatrix,self.basisNotMatrix,self.bVector)

        # 根据表示矩阵和优化目标得到目标函数用非决策变量表示的表达式的系数向量
        self.jueCeMax=self._CalculateMuBiaoFunction(self.jueCeWhichOutputMatrix,self.target)
        # print("line61:{}".format(self.jueCeMax))

        # 初始化完毕，退出
        return 0

    def run(self,MaxTime=999):
        '''
        调用此方法来使用单纯形法

        :return:
        '''
        if self.IterationInit()==-1: #用来对各个迭代变量初始化并获得一组初始解
            print('程序退出\n')
            return -1
        # 开始迭代
        ii=0
        while True:
            print('='*100)
            #print('#'*8+'第{}次迭代'.format(ii)+'#'*8)
            #print(self.jueCeMax[0,1:], self.basisNotVaries)
            # 找到最优解
            if (self.jueCeMax[0,1:]< 0).all():
                print("已经找到唯一的最优解")
                print('最优的决策向量为:\n{0},\n此时的目标函数最大值为:\n{1}\n'.format(self.bestSolution,self.jueCeMax[0,0]))
                return self.bestSolution

            elif (self.jueCeMax[0,1:]<=0).all() & (self.jueCeMax[0,1:]>-self.yibusinong).any():
                print('存在无数个最优解\n')
                print('已知的一组为:\n{0},\n此时的目标函数最大值为:\n{1}\n'.format(self.bestSolution, self.jueCeMax[0]))
                return self.bestSolution
            else:# 没有找到最优解
                # 找寻将被换入的基向量

                self.selectInput=self._findInputInJueCeVector(self.jueCeMax,self.basisNotVaries)
                # 找寻将被换出的基变量

                self.selectOutput=self._findOuputInBiaoShiMatrix(self.jueCeWhichOutputMatrix,self.basisVaries,
                                                             self.selectInput,self.basisNotVaries)

                # 利用选入选出更新四个参数,这里完全可以不加返回值,我主要是怕自己弄混才这样写的,这样写不好.
                self.basisMatrix,self.basisNotMatrix,\
                self.basisVaries,self.basisNotVaries=self._OperateInputOutput(self.selectInput,self.selectOutput)

                # 进行初等航变换,求出此时的解,:
                self.bestSolution = self._CalculateSolutionVector(self.basisVaries, self.basisNotVaries,
                                                                  self.basisMatrix, self.basisNotMatrix, self.bVector)
                #print('当前迭代次数下的最优解为:\n{0},最优值为:\n{1}'.format(self.bestSolution,self.jueCeMax[0,0]))
                # 根据求得的解得到非决策变量表示决策变量的形式(非常简单)
                self.jueCeWhichOutputMatrix = self._CalculateBiaoShiMatrix(self.basisMatrix, self.basisNotMatrix,
                                                                           self.bVector)
                # print('当前的非基变量表示基变量的矩阵为:{}'.format(self.jueCeWhichOutputMatrix))
                # 根据表示矩阵和优化目标得到目标函数用非决策变量表示的表达式的系数向量
                self.jueCeMax = self._CalculateMuBiaoFunction(self.jueCeWhichOutputMatrix, self.target)

                # 迭代一次完成
                ii+=1
                if ii==MaxTime:
                    print('最大迭代次数已达到,程序退出\n')
                    print('目前最优的决策向量为:\n{0},\n此时的目标函数最大值为:\n{1}\n'.format(self.bestSolution,self.jueCeMax[0]))
                    return self.bestSolution

    def _max(self,a,b):
        if a>=b:
            return a,1
        else:
            return b,-1

    def _findInputInJueCeVector(self,jVector,bNVec):
        _,m=jVector.shape
        max=self.yibusinong
        InputIndex=-100

        for i in range(m):

            if i==0: # 第一项是目标函数最大值,直接跳过
                continue
            if jVector[:,i]<=0:

                continue
            else:
                max,tihuan=self._max(jVector[:,i],max)

                if tihuan==1:
                    InputIndex=i

        return bNVec[InputIndex-1]

    def _min(self,a,b):
        if a<=b:
            return a,1
        else:
            return b,-1

    def _findOuputInBiaoShiMatrix(self,BiaoShiMatrix,basisV,selectIp,basisNV):
        m,n=BiaoShiMatrix.shape
        # 需要得到选入变量对应的那一列的索引

        XuanRuIndex=self._getIndexOfAVector(basisNV,selectIp)
        # 存储最小的选出变量的变化的数值
        min=np.inf
        OutputIndex=-100
        print(m)
        for i in range(m): # 对每一个基变量进行计算,得到最应该选出的基变量
            print('============2222222==========================')

            if self._equation(BiaoShiMatrix[i,XuanRuIndex+1],self.yibusinong)==1:
                a=0
            else:
                a=BiaoShiMatrix[i,0]*(-1)/(BiaoShiMatrix[i,XuanRuIndex+1])
            #print(a)
            if a<=0: # 保证向量的恒正
                print('该问题无界,无解.\n')
                continue
            #   print(a)
            min,tihuan=self._min(a,min)
            #print(min,tihuan,OutputIndex)
            if tihuan==1:
                OutputIndex=i
            print('i:{0},a:{1},min:{2},OutputIndex:{3}.'.format(i,a,min,OutputIndex))

        return basisV[OutputIndex]

    def _OperateInputOutput(self,input,output):

        RuIndex=self._getIndexOfAVector(self.basisNotVaries, input)
        ChuIndex=self._getIndexOfAVector(self.basisVaries, output)
        RuVec=self.basisNotMatrix[:,RuIndex]
        ChuVec=self.basisMatrix[:,ChuIndex].copy()

        self.basisMatrix[:,ChuIndex]=RuVec

        self.basisNotMatrix[:,RuIndex]=ChuVec
        self.basisVaries[ChuIndex]=input
        self.basisNotVaries[RuIndex]=output

        # 这里完全不需要返回,加上为明显
        return self.basisMatrix,self.basisNotMatrix,self.basisVaries,self.basisNotVaries

    def _equation(self,x,y):
        if abs(x-y)<=self.yibusinong:
            return True
        else:
            return False

    def _getIndexOfAVector(self,vector,value):
        length,=vector.shape
        for i in range(length):
            if vector[i]==value:
                return i
        print('没有找到参数!\n')
        return None

    def _CalculateSolutionVector(self,basisV,basisNV,basisM,basisNM,bVec):
        m,n=basisM.shape

        if (np.abs(basisM-np.eye(m))<=self.yibusinong).all(): # 单位矩阵
            vector=bVec
        else: # 如果不是就需要使用初等行变换啊啊啊啊啊啊啊啊啊啊!麻烦死了!
            if self.selectInput==-100 & self.selectOutput==-100: # 也就是没有进行选入选出操作
                pass # 突然想到这是不会发生的
            else:
                # 进行初等行变换把左侧变为单位阵
                changeIndex=self._getIndexOfAVector(basisV,self.selectInput)

                if basisM[changeIndex,changeIndex]!=1:
                    x=basisM[changeIndex,changeIndex]
                    self.basisMatrix[changeIndex,:]/=x # 首先将对角元素进行归一化
                    self.bVector[changeIndex]/=x
                    self.basisNotMatrix[changeIndex,:]/=x

                # 使用那一行的数值乘以若干倍加在其他行上,使得整体变为单位阵
                for iii in range(m):
                    if iii==changeIndex:
                        continue # 已经处理过,考虑的是将这一行的数乘加到其他行
                    else:

                        if not self._equation(basisM[iii,changeIndex],0):
                            xx=self.basisMatrix[iii,changeIndex]
                            #print(xx)
                            self.basisMatrix[iii,:]-=self.basisMatrix[changeIndex,:]*xx
                            self.bVector[iii]-=self.bVector[changeIndex]*xx
                            self.basisNotMatrix[iii,:]-=self.basisNotMatrix[changeIndex]*xx
            vector=self.bVector

        # 通过添加0得到输出最优解的向量
        long,_=self.bestSolution.shape
        SolutionVector=self.bestSolution.copy()
        for i in range(long):
            if (i == self.basisVaries).any():
                index=self._getIndexOfAVector(self.basisVaries,i)

                # print('index is'+str(index))
                SolutionVector[i]=vector[index]
            else:
                SolutionVector[i]=0.0

        return SolutionVector
    
    def _CalculateBiaoShiMatrix(self,basisM,basisNM,bVec):
        m,n=basisM.shape
        if (np.abs(basisM-np.eye(m))<=self.yibusinong).all():
            AfterYiXiang=basisNM*(-1)
            matrix=np.zeros((self.x,1+self.y-self.x))
            matrix[:,0]=bVec
            matrix[:,1:]=AfterYiXiang
            return matrix
        else:
            print("出现错误.在用非基变量表示基变量时.你上一步没有处理好.\n")
            return None

    def _CalculateMuBiaoFunction(self,BiaoShiMatrix,MuBiaoVector):
        BiaoShiVector=np.zeros((1, self.y-self.x+1))

        for Variabelbasis in self.basisVaries:
            if MuBiaoVector[Variabelbasis]!=0:
                BiaoShiVector+=MuBiaoVector[Variabelbasis]*\
                               BiaoShiMatrix[self._getIndexOfAVector(self.basisVaries,Variabelbasis),:]

        for VariabelNotbasis in self.basisNotVaries:
            if MuBiaoVector[VariabelNotbasis]!=0:
                BiaoShiVector[0,self._getIndexOfAVector(self.basisNotVaries,VariabelNotbasis)+1]+=MuBiaoVector[VariabelNotbasis]
        return BiaoShiVector



'''
可以直接运行 python DanChunXingFa.py 运行下面的程序
'''


print('==========================示例1=================================')
test1=DanChunXingFa(usingAimMax=True,target=np.array([3,-1,-1]),XiaoYu=np.array([[1,-2,1,11]]),
                    Dayu=np.array([[-4,1,2,3]]),
                    equ=np.array([[-2,0,1,1]]) )
test1.run()

print('===========================示例二================================')
test2=DanChunXingFa(usingAimMax=True,target=np.array([[40],[45],[24]]),XiaoYu=np.array([[2,3,1,100],
                                                                                [3,3,2,120]]))
test2.run()

print('============================示例三================================')





import torch
import numpy as np


def liangzi(haha):
    liangziii=