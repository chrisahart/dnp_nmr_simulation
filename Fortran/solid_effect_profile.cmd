:: Temporary makefile

gfortran -g -pg -o solid_effect expokit.o mataid.o functions.f90 interactions.f90 solid_effect_dynamics.f90 solid_effect_main.f90 -O3 -lgomp -llapack -lblas

solid_effect.exe

gprof solid_effect.exe > solid_effect_profile.out

gprof solid_effect.exe  | python gprof2dot.py | dot -Tpng -o solid_effect_profile.png

rm functions.mod
rm interactions.mod
rm solid_effect_dynamics.mod
rm solid_effect.exe
rm gmon.out

PAUSE