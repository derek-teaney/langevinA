make -f Makefile test_random.exe
mpiexec -n 4 ./test_random.exe -log_view
