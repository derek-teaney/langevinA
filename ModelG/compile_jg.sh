#!  /usr/bin/bash

source Makefiles/setupprlm.sh
make -f Makefiles/Makefile.prlm clean
make -f Makefiles/Makefile.prlm SuperPions.exe
mv SuperPions.exe ../../
