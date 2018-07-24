:: Temporary makefile

gfortran -o solid_effect expokit.o mataid.o functions.f90 interactions.f90 solid_effect_dynamics.f90 solid_effect_main.f90 -fopenmp -lgomp -llapack -lblas -O2 -fcheck=all

solid_effect.exe

rm functions.mod
rm interactions.mod
rm solid_effect_dynamics.mod
rm solid_effect.exe

python solid_effect_plotting.py

PAUSE