:: This windows batch script creates .pyd files from fortran code, for both conda (python 3.6) and canopy (python 3.5)
:: Edit input/output filenames, conda and canopy paths as required

@echo off

SET input=f2py_dynamics
SET output=f2py_dynamics

SET canopy=C:\Users\storm\AppData\Local\Enthought\Canopy\edm\envs\User\python

%canopy% -m numpy.f2py -c -m %output% Kronecker.f90 %input%.f90 --fcompiler=gnu95 --compiler=mingw32 --f90flags="-fopenmp -O2" -lgomp -llapack

PAUSE
