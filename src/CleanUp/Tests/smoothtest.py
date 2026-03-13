import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py

plt.ioff()

def ls(name) :
    print(name)

def printh5():
    with h5py.File('smoothtest_ic.h5','r') as f:
        f.visit(ls)

with h5py.File('smoothtest1_phi.h5','r') as f:
    p1 = f['phi_final'][16,:,:,1]
    p9 = f['phi_final'][16,16,:,4]
    p10 = f['phi_final'][16,16,:,9]

with h5py.File('smoothtest2_phi.h5','r') as f:
    p2 = f['phi_final'][16,:,:,1]
    p7 = f['phi_final'][16,16,:,4]
    p8 = f['phi_final'][16,16,:,9]

with h5py.File('smoothtest3_phi.h5','r') as f:
    p3 = f['phi_final'][16,:,:,1]

with h5py.File('smoothtest3_phi.h5','r') as f:
    p4 = f['phi_final'][16,16,:,4]
    p5 = f['phi_final'][16,16,:,9]
    p6 = np.linspace(0,31, 32)

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6), (ax7,ax8, ax9)) = plt.subplots(3, 3) 

ax1.matshow(p1) ;
ax2.matshow(p2) ;
ax3.matshow(p3) ;
ax4.matshow(p1-p1) ;
ax5.matshow(p2-p1) ;
ax6.matshow(p3-p1) ;
ax7.plot(p6,p4, p6, p5, p6, p7, p6, p8, p6,p9, p6, p10) 
plt.show()  

