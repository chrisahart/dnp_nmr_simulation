:: This windows batch script creates .pyd files from fortran code

@echo off

SET output=f2py_cross_effect

python -m numpy.f2py -c -m %output% expokit.o mataid.o functions.f90 interactions.f90 cross_effect_dynamics.f90 --fcompiler=gnu95 --compiler=mingw32 --f90flags="-fopenmp -O2" -lgomp -llapack -lblas

PAUSE