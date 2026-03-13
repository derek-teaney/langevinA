#mpiexec -n 8 ./ModelA-Beuler.exe input=thermalizetest.in evolverType=1 outputfiletag=thermalizetest1 deltat=0.1 finaltime=20
#mpiexec -n 8 ./ModelA-Beuler.exe input=thermalizetest.in evolverType=2 outputfiletag=thermalizetest2 deltat=0.01 finaltime=80
#mpiexec -n 8 ./ModelA-Beuler.exe input=thermalizetest.in evolverType=3 outputfiletag=thermalizetest3 deltat=0.05 finaltime=20
#mpiexec -n 8 ./ModelA-Beuler.exe input=thermalizetest.in evolverType=4 outputfiletag=thermalizetest4 detlat=0.04 finaltime=80

python thermalizetest.py

