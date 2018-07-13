gfortran -o solid_effect expokit.o mataid.o functions.f90 interactions.f90 solid_effect_dynamics.f90 solid_effect_main.f90 -lgomp -llapack -lblas  -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan


solid_effect.exe

rm functions.mod
rm interactions.mod
rm solid_effect_dynamics.mod
rm solid_effect.exe

PAUSE