// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpm class

#include "cello.hpp"
#include "enzo.hpp"

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
    Field field = BLOCK->data()->field();				\
    enzo_float * f = (enzo_float *) field.values(FIELD);            \
    enzo_float * f_copy = (enzo_float *) field.values(FIELD_COPY);  \
    int mx,my,mz;                                                       \
    field.dimensions(0,&mx,&my,&mz);                                    \
    for (int i=0; i<mx*my*mz; i++) f_copy[i]=f[i];              \
  }

#define DEBUG_BLOCK "B00:1_00:0_00:1"
#define DEBUG_CYCLE 101

//----------------------------------------------------------------------

EnzoMethodPpm::EnzoMethodPpm ()
  : Method(),
    comoving_coordinates_(enzo::config()->physics_cosmology)
{
  // Initialize default Refresh object

  cello::simulation()->new_refresh_set_name(ir_post_,name());
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

  FieldDescr * field_descr = cello::field_descr();

  // add all colour fields to refresh
  for (int i = 0; i < field_descr->field_count(); i++){
    std::string name = field_descr->field_name(i);
    if (field_descr->groups()->is_in(name,"colour")){
      refresh->add_field(name);
    }
  }

   // PPM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPpm::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | comoving_coordinates_;
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
  if (block->is_leaf()) {
    EnzoBlock * enzo_block = enzo::block(block);
    TRACE_PPM ("BEGIN SolveHydroEquations");
    enzo_block->SolveHydroEquations 
      ( block->time(), block->dt(), comoving_coordinates_ );
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

double EnzoMethodPpm::timestep ( Block * block ) const throw()
{

  TRACE_PPM("timestep()");

  EnzoBlock * enzo_block = enzo::block(block);

  enzo_float cosmo_a = 1.0, cosmo_dadt=0.0;

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  ASSERT ("EnzoMethodPpm::timestep()",
	  "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
	  ! (comoving_coordinates_ && (cosmology == NULL)) );

  if (comoving_coordinates_) {

    cosmology->compute_expansion_factor
      (&cosmo_a, &cosmo_dadt,(enzo_float)enzo_block->time());
    
  }

  enzo_float dtBaryons = ENZO_HUGE_VAL;

  /* Compute the pressure. */

  const int in = cello::index_static();

  EnzoComputePressure compute_pressure
    (EnzoBlock::Gamma[in],comoving_coordinates_);
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
			&EnzoBlock::Gamma[in], &EnzoBlock::PressureFree[in], &cosmo_a,
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
