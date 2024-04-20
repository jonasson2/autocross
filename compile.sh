#!/bin/bash
EXE=run/redfit-x
OPTIM=-O3
gfortran -Jmodules -Imodules $OPTIM -o$EXE src/redfit-x.f90 src/modules.f90
echo wrote executable to $EXE
