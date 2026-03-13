import h5py 
import numpy as np

fh = h5py.File("x2ktest.h5", "w") 

size = 8
x = np.arange(size,dtype=np.float64)

print(x) 

data = np.zeros((2,3,size), dtype=np.float64) 

for k  in [0,1, 2]:
    data[0,k,:] = np.cos(2*np.pi*k*x/float(size))
    data[1,k,:] = np.sin(2*np.pi*k*x/float(size))
print(data)

fh.create_dataset("x2ktest", data = data )



