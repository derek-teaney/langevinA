import numpy as np
import matplotlib.pyplot as plt
import h5py


def test1(ax):
    x = np.loadtxt('t1_zero.out') 

    print(x[:,0])

    # the histogram of the data
    n, bins, patches = ax.hist(x[:,0], 30)
    print(patches)

def test2(ax):
    with h5py.File('q_heatbathtest_phi.h5','r') as f:
        p1 = f['phi_ic'][16,:,:,5]
    ax.matshow(p1)

def test3(ax):
    with h5py.File('q_heatbathtest_phi.h5','r') as f:
        p1 = f['phi_final'][16,:,:,5]
    ax.matshow(p1)

def test4(ax):
    with h5py.File('q_heatbathtest_phi.h5','r') as f:
        p1 = f['phi_final'][16,16,:,5]
        p2 = f['phi_ic'][16,16,:,5]
        p3 = np.linspace(0,31, 32)
    ax.plot(p3,p1,p3,p2)

def test5(ax):
    with h5py.File('q_heatbathtest_phi.h5','r') as f:
        p1 = f['phi_ic_sol'][16,:,:,5]
        p2 = f['phi_ic'][16,:,:,5]
    ax.matshow(p1)

def test6(ax):
    with h5py.File('q_heatbathtest_phi.h5','r') as f:
        p1 = f['phi_final_sol'][16,16,:,5]
        p2 = f['phi_final'][16,16,:,5]
        p3 = f['phi_ic_sol'][16,16,:,5]
        p4 = f['phi_ic'][16,16,:,5]
    p0 = np.linspace(0,31,32) 
    ax.plot(p0, p1,  p0, p2, '.', p0, p3) 

#fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2, 3) 

#test1(ax1)
#test2(ax2)
#test3(ax3)
#test4(ax4)
test6(plt) 

plt.show()
