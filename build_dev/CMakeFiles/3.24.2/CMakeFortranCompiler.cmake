set(CMAKE_Fortran_COMPILER "/opt/intel/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "19.1.0.20200306")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/opt/apps/gcc/8.3.0/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/opt/apps/gcc/8.3.0/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/intel/compilers_and_libraries_2020.1.217/linux/pstl/include;/opt/intel/compilers_and_libraries_2020.1.217/linux/daal/include;/opt/intel/compilers_and_libraries_2020.1.217/linux/tbb/include;/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/include;/opt/intel/compilers_and_libraries_2020.1.217/linux/ipp/include;/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/include/intel64;/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/include/icc;/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/include;/usr/local/include;/opt/apps/gcc/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include;/opt/apps/gcc/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include-fixed;/opt/apps/gcc/8.3.0/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/intel/debugger_2020/libipt/intel64/lib;/opt/intel/compilers_and_libraries_2020.1.217/linux/daal/lib/intel64_lin;/opt/intel/compilers_and_libraries_2020.1.217/linux/tbb/lib/intel64/gcc4.8;/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin;/opt/intel/compilers_and_libraries_2020.1.217/linux/ipp/lib/intel64;/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin;/opt/apps/gcc/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0;/opt/apps/gcc/8.3.0/lib64;/lib64;/usr/lib64;/opt/apps/gcc/8.3.0/x86_64-pc-linux-gnu/lib;/opt/apps/gcc/8.3.0/lib;/lib;/usr/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
