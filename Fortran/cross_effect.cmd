:: Temporary makefile

gfortran -o cross_effect expokit.o mataid.o functions.f90 interactions.f90 cross_effect_dynamics.f90 cross_effect_main.f90 -lgomp -fopenmp -llapack -lblas -O2 -fcheck=all -Wall

cross_effect.exe

rm functions.mod
rm interactions.mod
rm cross_effect_dynamics.mod
rm cross_effect.exe

python cross_effect_plotting.py

PAUSE