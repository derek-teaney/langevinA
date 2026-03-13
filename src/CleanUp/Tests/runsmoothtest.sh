

time mpiexec -n 8 ./smoothtest.exe input=smoothtest.in outputfiletag=smoothtest1 evolverType=1 
time mpiexec -n 8 ./smoothtest.exe input=smoothtest.in outputfiletag=smoothtest2 evolverType=2 
time mpiexec -n 8 ./smoothtest.exe input=smoothtest.in outputfiletag=smoothtest3 evolverType=3 
#time mpiexec -n 1 ./smoothtest.exe input=smoothtest.in outputfiletag=smoothtest4 evolverType=2 
