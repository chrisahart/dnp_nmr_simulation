:: This windows batch script creates .pyd files from fortran code

@echo off

SET input=f2py_dynamics
SET output=f2py_dynamics

:: SET canopy=C:\Users\storm\AppData\Local\Enthought\Canopy\edm\envs\User\python
SET canopy=C:\Users\Chris\AppData\Local\Enthought\Canopy\edm\envs\User\python

%canopy% -m numpy.f2py -c -m %output% %input%.f90 --fcompiler=gnu95 --compiler=mingw32 --f90flags="-fopenmp -O2" -lgomp -llapack

PAUSE
