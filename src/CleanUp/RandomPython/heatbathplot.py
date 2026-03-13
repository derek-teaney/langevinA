import numpy as np
import matplotlib.pyplot as plt
import csv
import measurements as msr

data= np.loadtxt('hcrit1_N032_m-04813_h010000_l10000_averages_ave.out')

print("The magnetization is %f +- %f with corrtime %f" %(data[0], data[1], data[2]*100.))

hbdata = msr.ConfResults('heatbathtest.h5',200,0)
hbdata.computeMag()
print("The new magnetization is %f +- %f\n", hbdata.mag, hbdata.magErr)

plt.plot(hbdata.phi[:,0], '.')
plt.axhline(y=0.2977)
plt.axhline(y=0.2980)
plt.axhline(y=0.2974)
plt.show()

