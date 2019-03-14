import DanChunXingFa
import numpy as np

test=DanChunXingFa.DanChunXingFa(target=np.array([2.,3.,0.,0.,0.]),constraint=np.array([[1.,2.,1.,0.,0.,8.],
                                                                                        [4.,0.,0.,1.,0.,16.],
                                                                                        [0.,4.,0.,0.,1.,12.]]))
test.run()


