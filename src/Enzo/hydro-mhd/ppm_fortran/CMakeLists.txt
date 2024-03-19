# See LICENSE_CELLO file for license and copyright information

# define the Enzo:ppm_fortran subcomponent
add_library(Enzo_ppm_fortran
  ppm_fortran.hpp # included for completeness (but probably not necessary)

  calcdiss.F
  calc_dt.F
  calc_eigen.F
  cicinterp.F
  cic_deposit.F
  cic_interp.F
  euler.F
  flux_hll.F
  flux_hllc.F
  flux_twoshock.F
  interp1d.F interp2d.F interp3d.F
  inteuler.F
  intpos.F
  intprim.F
  intvar.F
  int_lin3d.F
  ngpinterp.F
  pgas2d.F
  pgas2d_dual.F
  ppm_de.F
  tscint1d.F tscint2d.F tscint3d.F
  twoshock.F
  xeuler_sweep.F
  yeuler_sweep.F
  zeuler_sweep.F
)
add_library(Enzo::ppm_fortran ALIAS Enzo_ppm_fortran)

# the following links a LOT of unnecessary stuff!
target_link_libraries(Enzo_ppm_fortran PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_ppm_fortran PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_ppm_fortran PRIVATE ${Cello_TARGET_LINK_OPTIONS})
