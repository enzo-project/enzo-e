#!/bin/bash

HDF5=-DH5_USE_16_API 
ARCH="-DLINUX -fPIC"
INT=int32
FLOAT=float64

if [ $INT == "int32" ]; then
   CINT="-DSMALL_INTS"
   FINT="-DSMALL_INTS"
fi
if [ $INT == "int64" ]; then
   CINT="-DLARGE_INTS"
   FINT="-fdefault-integer-8 -DLARGE_INTS"
fi

CLIB="$ARCH $CINT  -DCONFIG_BFLOAT_8  -O2 -shared -DSHARED_LIBRARY -I/usr/local/hdf5/1.8.2s/include     -I."
FLIB="-fno-second-underscore -ffixed-line-length-132 -fdefault-real-8 -fdefault-double-8 $FINT -O2 -shared $ARCH -DH5_USE_16_API   -DCONFIG_BFLOAT_8 "

	g++  -c -o calculate_cooling_time.o $CLIB calculate_cooling_time.C
g++  -c -o calculate_cooling_time_table.o $CLIB calculate_cooling_time_table.C
g++  -c -o calculate_gamma.o $CLIB calculate_gamma.C
g++  -c -o calculate_pressure.o $CLIB calculate_pressure.C
g++  -c -o calculate_pressure_table.o $CLIB calculate_pressure_table.C
g++  -c -o calculate_temperature.o $CLIB calculate_temperature.C
g++  -c -o calculate_temperature_table.o $CLIB calculate_temperature_table.C
gfortran  -c -o calc_rates_g.o $FLIB calc_rates_g.F
gfortran  -c -o calc_tdust_1d_g.o $FLIB calc_tdust_1d_g.F
gfortran  -c -o cie_thin_cooling_rate_g.o $FLIB cie_thin_cooling_rate_g.F
gfortran  -c -o colh2diss_g.o $FLIB colh2diss_g.F
gfortran  -c -o coll_rates_g.o $FLIB coll_rates_g.F
gfortran  -c -o cool1d_cloudy_g.o $FLIB cool1d_cloudy_g.F
gfortran  -c -o cool1d_multi_g.o $FLIB cool1d_multi_g.F
gfortran  -c -o cool_multi_time_g.o $FLIB cool_multi_time_g.F
g++  -c -o initialize_chemistry_data.o $CLIB initialize_chemistry_data.C
g++  -c -o initialize_cloudy_data.o $CLIB initialize_cloudy_data.C
g++  -c -o initialize_UVbackground_data.o $CLIB initialize_UVbackground_data.C
g++  -c -o set_default_chemistry_parameters.o $CLIB set_default_chemistry_parameters.C
g++  -c -o solve_chemistry.o $CLIB solve_chemistry.C
g++  -c -o solve_chemistry_table.o $CLIB solve_chemistry_table.C
gfortran  -c -o solve_rate_cool_g.o $FLIB solve_rate_cool_g.F
g++  -c -o update_UVbackground_rates.o $CLIB update_UVbackground_rates.C
g++   -O2 -shared -o libgrackle.so *.o -L/usr/local/hdf5/1.8.2s/lib -lhdf5 -lz  -lgfortran  
