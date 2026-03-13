import numpy as np
import matplotlib.pyplot as plt
import csv
import measurements as msr

data1 = msr.ConfResults('thermalizetest1.h5',0,0)
data1.computeMag()

data2 = msr.ConfResults('thermalizetest2.h5',0,0)
data2.computeMag()

data3 = msr.ConfResults('thermalizetest3.h5',0,0)
data3.computeMag()

data4 = msr.ConfResults('thermalizetest4.h5',0,0)
data4.computeMag()

plt.plot(data1.phi[:,0])
plt.plot(data2.phi[:,0])
plt.plot(data3.phi[:,0])
plt.plot(data4.phi[:,0])
plt.show()
