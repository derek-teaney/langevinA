rm testrestart*_*.rng
rm testrestart*.h5
mpiexec -n 2 ./SuperPions.exe input=testrestart1.in
mpiexec -n 2 ./SuperPions.exe input=testrestart1.in
mpiexec -n 2 ./SuperPions.exe input=testrestart1.in
mpiexec -n 2 ./SuperPions.exe input=testrestart2.in
h5dump -d /phi testrestart1.h5  testrestart2.h5
h5dump -d /wallX_phi_5 testrestart1.h5  testrestart2.h5
