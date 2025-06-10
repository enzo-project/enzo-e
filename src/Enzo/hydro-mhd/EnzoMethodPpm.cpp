// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpm class

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/hydro-mhd/hydro-mhd.hpp"

#include "Enzo/hydro-mhd/ppm_fortran/ppm_fortran.hpp" // FORTRAN_NAME(calc_dt)

// #define DEBUG_PPM
// #define COPY_FIELDS_TO_OUTPUT

#ifdef DEBUG_PPM
#  define TRACE_PPM(MESSAGE)						\
  CkPrintf ("%s:%d TRACE_PPM %s %s\n",					\
	    __FILE__,__LINE__,block->name().c_str(),MESSAGE);		\
  fflush (stdout);
#else
#  define TRACE_PPM(MESSAGE) /* ... */
#endif

#define COPY_FIELD(BLOCK,FIELD,FIELD_COPY)                              \
  {                                                                     \
    Field field = BLOCK->data()->field();                               \
    enzo_float * f = (enzo_float *) field.values(FIELD);                \
    enzo_float * f_copy = (enzo_float *) field.values(FIELD_COPY);      \
    if (f_copy) {                                                       \
      int mx,my,mz;                                                     \
      field.dimensions(0,&mx,&my,&mz);                                  \
      for (int i=0; i<mx*my*mz; i++) f_copy[i]=f[i];                    \
    }                                                                   \
  }

//----------------------------------------------------------------------

EnzoMethodPpm::EnzoMethodPpm (bool store_fluxes_for_corrections,
                              ParameterGroup p)
  : Method(),
    comoving_coordinates_(enzo::cosmology() != nullptr),
    store_fluxes_for_corrections_(store_fluxes_for_corrections),
    diffusion_(p.value_logical("diffusion", false)),
    flattening_(p.value_integer("flattening", 3)),
    pressure_free_(p.value_logical("pressure_free", false)),
    steepening_(p.value_logical("steepening", false)),
    use_minimum_pressure_support_(p.value_logical
                                  ("use_minimum_pressure_support",false)),
    minimum_pressure_support_parameter_(p.value_integer
                                        ("minimum_pressure_support_parameter",
                                         100))
{
  this->set_courant(p.value_float("courant",1.0));

  // check compatability with EnzoPhysicsFluidProps
  EnzoPhysicsFluidProps* fluid_props = enzo::fluid_props();
  const EnzoDualEnergyConfig& de_config = fluid_props->dual_energy_config();
  ASSERT("EnzoMethodPpm::EnzoMethodPpm",
         "selected formulation of dual energy formalism is incompatible",
         de_config.is_disabled() || de_config.bryan95_formulation());

  // technically, the PPM solver doesn't directly use any of the functionality
  // implemented by EnzoEOSIdeal, (it implements the required functionality in
  // fortran), but this check is here to ensure that the EOS is handled
  // consistently by different Methods in a given Enzo-E simulation
  ASSERT("EnzoMethodPpm::EnzoMethodPpm",
         "PPM solver is currently incompatible with a non-ideal EOS",
         fluid_props->eos_variant().holds_alternative<EnzoEOSIdeal>());

  const int rank = cello::rank();

  cello::define_field("density");
  cello::define_field("total_energy");
  cello::define_field("internal_energy");
  cello::define_field("pressure");
  if (rank >= 1) {
      cello::define_field("velocity_x");
      cello::define_field("acceleration_x");
  }
  if (rank >= 2) {
      cello::define_field("velocity_y");
      cello::define_field("acceleration_y");
  }
  if (rank >= 3) {
      cello::define_field("velocity_z");
      cello::define_field("acceleration_z");
  }

  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("density");
  refresh->add_field("velocity_x");
  refresh->add_field("velocity_y");
  refresh->add_field("velocity_z");
  refresh->add_field("total_energy");
  refresh->add_field("internal_energy");
  refresh->add_field("pressure");
  refresh->add_field("acceleration_x");
  refresh->add_field("acceleration_y");
  refresh->add_field("acceleration_z");

  // add all color fields to refresh
  refresh->add_all_fields("color");

   // PPM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPpm::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | comoving_coordinates_;
  p | store_fluxes_for_corrections_;
  p | diffusion_;
  p | flattening_;
  p | pressure_free_;
  p | steepening_;
  p | use_minimum_pressure_support_;
  p | minimum_pressure_support_parameter_;
}

//----------------------------------------------------------------------

void EnzoMethodPpm::compute ( Block * block) throw()
{
  TRACE_PPM("BEGIN compute()");
#ifdef COPY_FIELDS_TO_OUTPUT
  const int rank = cello::rank();
  COPY_FIELD(block,"density","density_in");
  COPY_FIELD(block,"velocity_x","velocity_x_in");
  COPY_FIELD(block,"velocity_y","velocity_y_in");
  if (rank >= 3) COPY_FIELD(block,"velocity_z","velocity_z_in");
  COPY_FIELD(block,"total_energy","total_energy_in");
  COPY_FIELD(block,"internal_energy","internal_energy_in");
  COPY_FIELD(block,"pressure","pressure_in");
  COPY_FIELD(block,"acceleration_x","acceleration_x_in");
  COPY_FIELD(block,"acceleration_y","acceleration_y_in");
  if (rank >= 3) COPY_FIELD(block,"acceleration_z","acceleration_z_in");
#endif

  bool single_flux_array = true;
  if (store_fluxes_for_corrections_){
    Field field = block->data()->field();

    auto field_names = field.groups()->group_list("conserved");
    const int nf = field_names.size();
    std::vector<int> field_list;
    field_list.resize(nf);
    for (int i=0; i<nf; i++) {
      field_list[i] = field.field_id(field_names[i]);
    }

    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    block->data()->flux_data()->allocate(nx,ny,nz,field_list,single_flux_array);
  }

  if (block->is_leaf()) {

    EnzoBlock * enzo_block = enzo::block(block);

    // (this should go in interpolation / restriction not here)
    //
    // // restore energy consistency if dual energy formalism used
    //
    // if (enzo::fluid_props()->dual_energy_config().bryan95_formulation()) {
    //   int mx,my,mz;
    //   field.dimensions(0,&mx,&my,&mz);
    //   const int m = mx*my*mz;
    //   enzo_float * ie = (enzo_float*) field.values("internal_energy");
    //   enzo_float * te = (enzo_float*) field.values("total_energy");
    //   enzo_float * vxa = (enzo_float*) field.values("velocity_x");
    //   enzo_float * vya = (enzo_float*) field.values("velocity_y");
    //   enzo_float * vza = (enzo_float*) field.values("velocity_z");
    //   const int rank = cello::rank();
    //   if (rank == 1) {
    //     for (int i=0; i<m; i++) {
    //       const enzo_float vx = vxa[i]*vxa[i];
    //       te[i] = ie[i] + 0.5*(vx*vx);
    //     }
    //   } else if (rank == 2) {
    //     for (int i=0; i<m; i++) {
    //       const enzo_float vx = vxa[i]*vxa[i];
    //       const enzo_float vy = vya[i]*vya[i];
    //       te[i] = ie[i] + 0.5*(vx*vx+vy*vy);
    //     }
    //   } else if (rank == 3) {
    //     for (int i=0; i<m; i++) {
    //       const enzo_float vx = vxa[i]*vxa[i];
    //       const enzo_float vy = vya[i]*vya[i];
    //       const enzo_float vz = vza[i]*vza[i];
    //       te[i] = ie[i] + 0.5*(vx*vx+vy*vy+vz*vz);
    //     }
    //   }
    // }

    TRACE_PPM ("BEGIN SolveHydroEquations");

    const double time = block->state()->time();
    const double dt   = block->state()->dt();
    EnzoMethodPpm::SolveHydroEquations 
      ( *enzo_block, time, dt, comoving_coordinates_,
        single_flux_array, diffusion_, flattening_, pressure_free_,
        steepening_,
        use_minimum_pressure_support_,
        minimum_pressure_support_parameter_);

    TRACE_PPM ("END SolveHydroEquations");

  }

#ifdef COPY_FIELDS_TO_OUTPUT
  COPY_FIELD(block,"density","density_out");
  COPY_FIELD(block,"velocity_x","velocity_x_out");
  COPY_FIELD(block,"velocity_y","velocity_y_out");
  if (rank >= 3) COPY_FIELD(block,"velocity_z","velocity_z_out");
  COPY_FIELD(block,"total_energy","total_energy_out");
  COPY_FIELD(block,"internal_energy","internal_energy_out");
  COPY_FIELD(block,"pressure","pressure_out");
  COPY_FIELD(block,"acceleration_x","acceleration_x_out");
  COPY_FIELD(block,"acceleration_y","acceleration_y_out");
  if (rank >= 3) COPY_FIELD(block,"acceleration_z","acceleration_z_out");
#endif
  TRACE_PPM("END compute()");

  block->compute_done();

}

//----------------------------------------------------------------------

double EnzoMethodPpm::timestep ( Block * block ) throw()
{

  TRACE_PPM("timestep()");

  EnzoBlock * enzo_block = enzo::block(block);

  enzo_float cosmo_a = 1.0, cosmo_dadt=0.0;

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  ASSERT ("EnzoMethodPpm::timestep()",
	  "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
	  ! (comoving_coordinates_ && (cosmology == NULL)) );

  if (comoving_coordinates_) {

    enzo_float time = (enzo_float)enzo_block->state()->time();
    cosmology->compute_expansion_factor (&cosmo_a, &cosmo_dadt, time);

  }

  enzo_float dtBaryons = ENZO_HUGE_VAL;

  /* Compute the pressure. */

  enzo_float gamma = enzo::fluid_props()->gamma();

  EnzoComputePressure compute_pressure(gamma, comoving_coordinates_);
  compute_pressure.compute(enzo_block);

  Field field = enzo_block->data()->field();

  int rank = cello::rank();

  enzo_float * density    = (enzo_float *)field.values("density");
  enzo_float * velocity_x = (rank >= 1) ?
    (enzo_float *)field.values("velocity_x") : NULL;
  enzo_float * velocity_y = (rank >= 2) ?
    (enzo_float *)field.values("velocity_y") : NULL;
  enzo_float * velocity_z = (rank >= 3) ?
    (enzo_float *)field.values("velocity_z") : NULL;
  enzo_float * pressure = (enzo_float *) field.values("pressure");

  /* calculate minimum timestep */

  int pressure_free_int = pressure_free_;

  FORTRAN_NAME(calc_dt)(&rank,
			enzo_block->GridDimension,
			enzo_block->GridDimension+1,
			enzo_block->GridDimension+2,
			enzo_block->GridStartIndex,
			enzo_block->GridEndIndex,
			enzo_block->GridStartIndex+1,
			enzo_block->GridEndIndex+1,
			enzo_block->GridStartIndex+2,
			enzo_block->GridEndIndex+2,
			&enzo_block->CellWidth[0],
			&enzo_block->CellWidth[1],
			&enzo_block->CellWidth[2],
			&gamma, &pressure_free_int, &cosmo_a,
			density, pressure,
			velocity_x,
			velocity_y,
			velocity_z,
			&dtBaryons);

  TRACE1 ("dtBaryons: %f",dtBaryons);

  dtBaryons *= courant_;

  double dt = dtBaryons;

  return dt;
}
