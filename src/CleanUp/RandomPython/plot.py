import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py


f = h5py.File('output.h5','r')  # open the file

solution = f["Timestepsolution/solution"]  #grab the solution
phi=f["Timestepsolution/phi"]
c00=f["Timestepsolution/c00"]

print(solution)  # printout what it is (some hdf5 type)
print(phi)       # printout what it is (some hdf5 type)

timesteps_p0 = solution[20,:,:,8,1]  # read the dataset from the file into a numpy matrix

plt.plot(phi[:,4])          # plot the magnetization
plt.plot(c00[20,:])          # plot the spatial correlation function
plt.contourf(timesteps_p0)  # make a contour plot
plt.show()



