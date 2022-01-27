# Changelog

## [Unreleased](https://github.com/leifdenby/uclales/tree/HEAD...v2.0.0)

**NOTE**: this changelog represents changes since commit
[370af9c9f1ffd244adf8d4ada460e2368538359d](https://github.com/leifdenby/uclales/commit/370af9c9f1ffd244adf8d4ada460e2368538359d)
which has been tagged on the [leifdenby/uclales
fork](https://github.com/leifdenby/uclales) as `v2.0.0`


[Full Changelog](https://github.com/leifdenby/uclales/compare/v2.0.0...HEAD)

*bugfixes*

- Ensure valid date for time-coordinate origin in the netCDF output
  [a48f430](https://github.com/leifdenby/uclales/commit/a48f4302aa33711c65e32c2ec5b481b1ac2d2ce8)

*improvements*

- Fixes for compilation on UK's ARCHER HPC system (which uses Cray Fortran) and
  cmake configuration for ARCHER
  [0e185ec](https://github.com/leifdenby/uclales/commit/0e185ecb975eba0425c407c23b3f77af5e9dda24)

- Fixes to support newer gfortran compilers (tested with `7.5.0` and `9.3.0`)
  [\#2](https://github.com/leifdenby/uclales/pull/2)

*maintenance*

- Clean-up CMakeLists.txt to use cmake to detect packages, include MPI wrapper
  and source-files in `src/` explicitly
  [\#4](https://github.com/leifdenby/uclales/pull/4)

- remove unused namelist variable `dtcross` (correct namelist variable for
  setting cross-section output frequency is `frqcross`)
  [2f85744](https://github.com/leifdenby/uclales/commit/2f85744cb92d333d2151202416811d88afa8c128)

- add makefiles for `reducets` and `reduceps` post-processing fortran programs
  [b9cb6d8](https://github.com/leifdenby/uclales/commit/b9cb6d81e445a451e473c965969f57eaa864ef13)


## [v2.0.0](https://github.com/convml/convml_tt/tree/v2.0.0)

*NOTE* the UCLA-LES source code states that the version in `2.0` in the [very
first
commit](https://github.com/uclales/uclales/blame/e690bbe2eb4c6737ccf7243ca02e6b3c375399c4/src/ncio.f90#L55).
There is no record of changes besides the commit history.
[370af9c9f1ffd244adf8d4ada460e2368538359d](https://github.com/leifdenby/uclales/commit/370af9c9f1ffd244adf8d4ada460e2368538359d)
was tagged as `v2.0.0` so that changes since that commit may be recorded. In
tagging this commit as `v2.0.0` this fork of `UCLA-LES` will be adopting
[semver](https://semver.org/) for versioning going forward.
