:: This windows batch script creates .pyd files from fortran code

@echo off

SET input=f2py_dynamics
SET output=f2py_dynamics

python -m numpy.f2py -c -m %output% expokit.o mataid.o functions.f90 interactions.f90 solid_effect_dynamics.f90 --fcompiler=gnu95 --compiler=mingw32 --f90flags="-fopenmp -O2" -lgomp -llapack -lblas

PAUSE