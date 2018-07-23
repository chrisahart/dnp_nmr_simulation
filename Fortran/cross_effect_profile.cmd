:: Temporary makefile

gfortran -pg -o cross_effect expokit.o mataid.o functions.f90 interactions.f90 cross_effect_dynamics.f90 cross_effect_main.f90 -lgomp -llapack -lblas

cross_effect.exe

gprof cross_effect.exe > cross_effect_profile.out

rm functions.mod
rm interactions.mod
rm cross_effect_dynamics.mod
rm cross_effect.exe
rm gmon.out

PAUSE