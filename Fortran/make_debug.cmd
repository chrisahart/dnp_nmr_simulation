gfortran -o solid_effect expokit.o mataid.o f2py_functions.f90 f2py_dynamics.f90 solid_effect_main.f90 -lgomp -llapack -lblas -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan


solid_effect.exe

::rm f2py_functions.mod
::rm f2py_dynamics.mod
::rm solid_effect.exe

PAUSE