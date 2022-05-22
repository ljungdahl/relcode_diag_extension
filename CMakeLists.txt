cmake_minimum_required(VERSION 3.16)

project(relcode_diag_extension)
enable_language(Fortran)


set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(MKL_DIR /home/ljungdahl/intel/mkl)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

include_directories(${MKL_DIR}/include/intel64/lp64)
include_directories(${MKL_DIR}/include)

set(CMAKE_C_COMPILER gcc-9)
set(CMAKE_Fortran_COMPILER gfortran-9)

set(CMAKE_Fortran_FLAGS
        "-Og -fopenmp -cpp -ffree-line-length-none -fbacktrace -fbounds-check")


set(LIBMKL ${MKL_DIR}/lib/intel64/libmkl_blas95_lp64.a
        ${MKL_DIR}/lib/intel64/libmkl_lapack95_lp64.a
        -Wl,--start-group
        ${MKL_DIR}/lib/intel64/libmkl_gf_lp64.a
        ${MKL_DIR}/lib/intel64/libmkl_gnu_thread.a
        ${MKL_DIR}/lib/intel64/libmkl_core.a
        -Wl,--end-group
        -lpthread
        -lm
        -ldl)

set(LFLAGS -fopenmp -m64)

link_libraries(${LIBMKL} ${LFLAGS})

# fixed format .f
set(FIXEDSRC 
        src/coulcc.f
        src/coul_gen_cc_rel.f)

# free format .f90
set(F90SRC
        src/main.f90)

# fortran module source files *_mod.f90
set(MODSRC
        src/perturbed_wavefunction_mod.f90
        src/photons_mod.f90
        src/in_data_mod.f90
        src/intrinsics_mod.f90
        src/types_mod.f90
        src/ionisation_mod.f90
        src/output_data_mod.f90
        src/utility_mod.f90)

# this is only for IDE purposes to parse the macros.
set(PREPROCESSOR_MACRO_SRC src/logger.h)


set(OUT_EXEC relcode_diag_extension.exe)
add_executable(${OUT_EXEC} ${MODSRC} ${FIXEDSRC} ${F90SRC} ${PREPROCESSOR_MACRO_SRC})
