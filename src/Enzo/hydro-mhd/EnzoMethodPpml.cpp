// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpml class

//----------------------------------------------------------------------

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/hydro-mhd/hydro-mhd.hpp"

#include "Enzo/hydro-mhd/ppml_fortran/ppml_fortran.hpp" // FORTRAN_NAME(calc_dt_ppml)

//----------------------------------------------------------------------

EnzoMethodPpml::EnzoMethodPpml(ParameterGroup p)
  : Method(),
    comoving_coordinates_(enzo::cosmology() != nullptr)
{
  // Initialize the default Refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());

  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  // check compatability with EnzoPhysicsFluidProps
  EnzoPhysicsFluidProps* fluid_props = enzo::fluid_props();
  const EnzoDualEnergyConfig& de_config = fluid_props->dual_energy_config();
  ASSERT("EnzoMethodPpml::EnzoMethodPpml",
         "incompatible with dual energy formalism", de_config.is_disabled());

  // technically, the PPML solver doesn't directly use any of the functionality
  // implemented by EnzoEOSIsothermal, (it implements the required
  // functionality in fortran), but this check is here to ensure that the EOS
  // is handled consistently by different Methods in a given Enzo-E simulation
  ASSERT("EnzoMethodPpml::EnzoMethodPpml",
         "PPML solver is currently incompatible with a non-isothermal EOS",
         fluid_props->eos_variant().holds_alternative<EnzoEOSIsothermal>());

  this->set_courant(p.value_float("courant",1.0));
}

//----------------------------------------------------------------------

void EnzoMethodPpml::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
  p | comoving_coordinates_;
}

//----------------------------------------------------------------------

void EnzoMethodPpml::compute ( Block * block ) throw()
{

  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = enzo::block(block);

  EnzoMethodPpml::SolveMHDEquations ( *enzo_block, block->state()->dt() );

  enzo_block->compute_done();

}

//----------------------------------------------------------------------

double EnzoMethodPpml::timestep (Block * block) throw()
{
 
  EnzoBlock * enzo_block = enzo::block(block);

  /* initialize */
 
  enzo_float dt;
  // enzo_float dtTemp;
  enzo_float dtBaryons      = ENZO_HUGE_VAL;
  // enzo_float dtViscous      = ENZO_HUGE_VAL;
  // enzo_float dtParticles    = ENZO_HUGE_VAL;
  // enzo_float dtExpansion    = ENZO_HUGE_VAL;
  // enzo_float dtAcceleration = ENZO_HUGE_VAL;
  // int dim, i, result;
 
  /* Compute the field size. */
 
  // int size = 1;
  // for (dim = 0; dim < GridRank; dim++)
  //   size *= GridDimension[dim];
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  enzo_float cosmo_a = 1.0, cosmo_dadt = 0.0;
  
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  ASSERT ("EnzoMethodPpml::timestep()",
	  "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
	  ! (comoving_coordinates_ && (cosmology == NULL)) );

  if (cosmology) {
    const auto time = enzo_block->state()->time();
    cosmology->compute_expansion_factor (&cosmo_a, &cosmo_dadt, time);
  }
  //  float afloat = float(a);
 
  /* 1) Compute Courant condition for baryons. */

  Field field = enzo_block->data()->field();
  if (field.num_permanent() > 0) { // TODO: revisit if-clause. This could be
                                   // improved. (plus we probably want to
                                   // report an error when false)
 
    /* Find fields: density, total energy, velocity1-3. */
 
    // int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    // float *pressure_field;

    // if (HydroMethod != PPML_Isothermal3D) {
    //   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
    // 					 Vel3Num, TENum) == FAIL) {
    // 	fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
    // 	exit(FAIL);
    //   }
 
    // /* Compute the pressure. */
 
    //   pressure_field = new float[size];
    //   if (DualEnergyFormalism)
    // 	result = this->ComputePressureDualEnergyFormalism(Time, pressure_field);
    //   else
    // 	result = this->ComputePressure(Time, pressure_field);
 
    //   if (result == FAIL) {
    // 	fprintf(stderr, "Error in grid->ComputePressure.\n");
    // 	exit(EXIT_FAILURE);
    //   }
    // }

    /* Call fortran routine to do calculation. */

    enzo_float * d  = (enzo_float *) field.values("density");
    enzo_float * vx = (enzo_float *) field.values("velox");
    enzo_float * vy = (enzo_float *) field.values("veloy");
    enzo_float * vz = (enzo_float *) field.values("veloz");
    enzo_float * bx = (enzo_float *) field.values("bfieldx");
    enzo_float * by = (enzo_float *) field.values("bfieldy");
    enzo_float * bz = (enzo_float *) field.values("bfieldz");

    FORTRAN_NAME(calc_dt_ppml)
      (enzo_block->GridDimension, 
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
       d,
       vx, vy, vz,
       bx, by, bz,
       &dtBaryons);
    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
 
    dtBaryons *= courant_;
    
  }
 
  /* 5) calculate minimum timestep */

  dt = std::numeric_limits<enzo_float>::max();

  dt = MIN(dt, dtBaryons);


  return dt;
}

//----------------------------------------------------------------------
