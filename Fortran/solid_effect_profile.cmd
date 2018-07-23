:: Temporary makefile

gfortran -pg -o solid_effect expokit.o mataid.o functions.f90 interactions.f90 solid_effect_dynamics.f90 solid_effect_main.f90 -lgomp -llapack -lblas

solid_effect.exe

gprof solid_effect.exe > solid_effect_profile.out

:: rm functions.mod
:: rm interactions.mod
:: rm solid_effect_dynamics.mod
:: solid_effect.exe
:: gmon.out

PAUSE