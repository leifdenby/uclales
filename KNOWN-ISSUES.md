- Cannot currently be compiled without MPI support: Dummy-interface for when
  MPI is disabled doesn't work with `src/particles.f90` and so UCLA-LES can't
  be compiled without MPI at the moment.

- Highest version of gfortran compiler that is verified to work is `9.3.0`.
  There are type missmatches for subroutine calls in `src/rfft.f90` which raise
  warnings for gfortran `7.5.0` and `9.3.0`, but are raised as error in `11.2.0`.
