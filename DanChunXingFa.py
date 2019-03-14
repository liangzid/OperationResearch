import numpy as np

class DanChunXingFa:
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

    def __init__(self,target,constraint):
        self.target=target # 决策目标，必须是最大化问题
        self.constraintA=constraint[:,:-1] # 约束条件的系数矩阵
        #print(self.constraintA)
        self.constraintb=constraint[:,-1] # 约束方程的b向量
        self.x,self.y=self.constraintA.shape

        # 这些是在迭代时需要被频繁使用的变量
        self.bestSolution=np.zeros((self.y,1)) # 存储历史最优解（由于贪心算法，所以历史最优解就是当前解）的向量

        self.basisVaries=np.zeros((self.x,1))  # 存储基变量的索引，从零开始
        self.basisNotVaries=np.zeros((self.y-self.x,1)) # 存储非基变量的索引

        self.basisMatrix=np.eye(self.x) # 基向量形成的方阵
        self.bVector=np.zeros((self.x,1)) # 每次迭代过程中的b向量
        self.basisNotMatrix=np.zeros((self.x,self.y-self.x)) # 非基变量对应的系数的矩阵

        self.jueCeMax=np.zeros((1,self.y-self.x+1)) # 使用非决策变量表示目标函数得到的系数，其中首位为常数项（即为目标最大值），
                                                    # 剩下的顺序按照非决策变量的顺序进行
        self.maxValueNow=self.jueCeMax[0]
        self.selectInput=-100 # 决策使用那个变量作为选入
        self.jueCeWhichOutputMatrix = np.zeros((self.x, self.y - self.x + 1))  # 存储考虑选出哪一个变量时进行的操作
        self.selectOutput=-100 # 决策使用那个变量作为选出

        self.yibusinong=1e-3 # 常量

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
        for i in range(m): # 对每一个基变量进行计算,得到最应该选出的基变量
            #print(i)
            if self._equation(BiaoShiMatrix[i,XuanRuIndex+1],self.yibusinong)==1:
                a=0
            else:
                a=BiaoShiMatrix[i,0]*(-1)/(BiaoShiMatrix[i,XuanRuIndex+1])
            if a<=0: # 保证向量的恒正
                continue
            #   print(a)
            min,tihuan=self._min(a,min)
            if tihuan==1:
                OutputIndex=i

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
        BiaoShiVector=np.zeros((1, 3))

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

#==================================例子============================
test=DanChunXingFa(target=np.array([2.,3.,0.,0.,0.]),constraint=np.array([[1.,2.,1.,0.,0.,8.],
                                                                          [4.,0.,0.,1.,0.,16.],
                                                                          [0.,4.,0.,0.,1.,12.]]))
test.run()






