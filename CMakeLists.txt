############################
# MESOHH CMAKE FILE
# Supports (thusfar) GNU, Intel and IBM compilers
# To set the compiler to a specific one, set the FC environment variable
# NETCDF can be set in the config/default.cmake, or through command line options
# NETCDF_INCLUDE, NETCDF_LIB, FFTW_INCLUDE and FFTW_LIB, respectively.
# Command line options are fed to CMake using the -D switch
# Options include:
# * MPI (TRUE/FALSE) to toggle the use of MPI. True by default.
# * PROFILER (SCALASCA/MARMOT) to enable a specific profiler. Default is none.
# * CMAKE_BUILD_TYPE (RELEASE/DEBUG) to build in optimized or debug mode. Default is RELEASE
# All options are persistent in the sense that once set, they will last until changed by the user.
# For daily use "cmake .. && make" should suffice to build the code, but a command line could look like:
# export FC="gfortran" && cmake -D SYST=tornado -D MPI=TRUE -D PROFILER=SCALASCA -D CMAKE_BUILD_TYPE=DEBUG .. && make
############################
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/config)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

cmake_minimum_required (VERSION 3.3.0)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


#  make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE )
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
  set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
      "Choose the type of build, options are: None Debug Release."  FORCE)
endif ()

#Load system specific settings
if (NOT SYST)
  set (SYST default CACHE STRING
      "Choose the location: mpipc thunder fedora archlinux supermuc jugene" FORCE )
endif ()
include (${SYST} OPTIONAL)

# Start the project only after all the variables are set
project (uclales Fortran)
enable_language(Fortran)
MESSAGE(STATUS "Fortran Compiler " ${CMAKE_Fortran_COMPILER})

## netCDF
find_package(NetCDF COMPONENTS Fortran REQUIRED)
set(INCLUDE_DIRS ${netCDF_INCLUDE_DIR})


## MPI
# Set MPI on or off. By default it is on.
set (MPI TRUE CACHE BOOL "Choose to enable MPI or not.")
if(MPI)
 MESSAGE(STATUS "MPI is Enabled.")
 find_package(MPI REQUIRED)
 add_definitions(${MPI_Fortran_COMPILE_FLAGS})
 include_directories(${MPI_Fortran_INCLUDE_PATH})

 set(INCLUDE_DIRS ${INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_DIRS})
 set(MPI_INTERFACE_SOURCEFILE "src/mpi/mpi_interface.f90")
else()
 MESSAGE(STATUS "MPI is Disabled.")
 # load source-file which has dummy subroutines for mpi calls
 set(MPI_INTERFACE_SOURCEFILE "src/seq/mpi_interface.f90")
endif()

# Compiler flags
set(CMAKE_Fortran_FLAGS ${USER_Fortran_FLAGS} )
set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} )
set(CMAKE_Fortran_FLAGS_DEBUG ${USER_Fortran_FLAGS_DEBUG})
MESSAGE(STATUS "Build Type " ${CMAKE_BUILD_TYPE})
if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE})
else()
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG})
endif()


FILE(GLOB sourcefiles ${MPI_INTERFACE_SOURCEFILE} "src/*.?90")
add_executable(${PROJECT_NAME} ${sourcefiles})
include_directories(${INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE NetCDF::NetCDF_Fortran)
if(MPI)
  target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_Fortran)
endif()


# For printing out all cmake variables:
set(PRINT_CMAKE_VARS FALSE)
if(PRINT_CMAKE_VARS)
  get_cmake_property(_variableNames VARIABLES)
  list (SORT _variableNames)
  foreach (_variableName ${_variableNames})
      message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endif()
