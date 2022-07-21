// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpmlIG.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @author   Alexei Kritsuk (akritsuk@ucsd.edu) 
/// @date     Tue Sep 18 14:16:01 PDT 2018
/// @brief    Implements the EnzoMethodPpmlIG class

//----------------------------------------------------------------------

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

// #define DEBUG_FIELDS

#ifdef DEBUG_FIELDS
#   define CHECK_FIELD(VALUES,NAME)             \
  ASSERT1("CHECK_FIELD",                        \
          "Field %s must be defined",           \
          NAME,                                 \
          (VALUES != nullptr));

#   define FIELD_STATS(NAME,VALUES,mx,my,mz,gx,gy,gz)           \
  {                                                             \
   double avg=0.0, max=-1.0, min=1e9;                           \
   int count=0;                                                 \
   for (int iz=gz; iz<mz-gz; iz++) {                            \
     for (int iy=gy; iy<my-gy; iy++) {                          \
       for (int ix=gx; ix<mx-gx; ix++) {                        \
         const int i=ix+mx*(iy+my*iz);                          \
         avg += VALUES[i];                                      \
         max = std::max(max,VALUES[i]);                         \
         min = std::min(min,VALUES[i]);                         \
         count++;                                               \
       }                                                        \
     }                                                          \
   }                                                            \
   avg /= count;                                                \
   CkPrintf ("FIELD_STATS %s  %g %g %g\n",NAME,min,avg,max);    \
  }
#else
#   define CHECK_FIELD(VALUES,NAME) /* ... */
#   define FIELD_STATS(NAME,VALUES,mx,my,mz,gx,gy,gz) /* ... */
#endif

//----------------------------------------------------------------------

EnzoMethodPpmlIG::EnzoMethodPpmlIG () 
  : Method(),
    comoving_coordinates_(enzo::config()->physics_cosmology),
    gamma_               (enzo::config()->field_gamma),
    b0_()
{
  b0_[0] = enzo::config()->method_ppml_b0[0];
  b0_[1] = enzo::config()->method_ppml_b0[1];
  b0_[2] = enzo::config()->method_ppml_b0[2];
  // Initialize the default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());

  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();
}

//----------------------------------------------------------------------

void EnzoMethodPpmlIG::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
  p | comoving_coordinates_;
  p | gamma_;
  PUParray(p,b0_,3);
}

//----------------------------------------------------------------------

void EnzoMethodPpmlIG::compute ( Block * block ) throw()
{

  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  enzo_block->SolveMHDEquationsIG ( block->dt(), gamma_, b0_ );

  enzo_block->compute_done();

}

//----------------------------------------------------------------------

double EnzoMethodPpmlIG::timestep (Block * block) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

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

  ASSERT ("EnzoMethodPpmlIG::timestep()",
	  "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
	  ! (comoving_coordinates_ && (cosmology == NULL)) );

  if (cosmology) {
    cosmology->compute_expansion_factor (&cosmo_a, &cosmo_dadt,enzo_block->time());
  }
  //  float afloat = float(a);
 
  /* 1) Compute Courant condition for baryons. */
 
  const int in = cello::index_static();
  if (EnzoBlock::NumberOfBaryonFields[in] > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    // int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    // float *pressure_field;

    // if (HydroMethod != PPML_Isothermal3D) {
    //   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
    // 					 Vel3Num, TENum) == FAIL) {
    // 	fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
    // 	exit(FAIL);

 
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
 
    Field field = enzo_block->data()->field();

    enzo_float * d  = (enzo_float *) field.values("density");
    enzo_float * vx = (enzo_float *) field.values("velox");
    enzo_float * vy = (enzo_float *) field.values("veloy");
    enzo_float * vz = (enzo_float *) field.values("veloz");
    enzo_float * bx = (enzo_float *) field.values("bfieldx");
    enzo_float * by = (enzo_float *) field.values("bfieldy");
    enzo_float * bz = (enzo_float *) field.values("bfieldz");
    enzo_float * pr = (enzo_float *) field.values("pressure");
    CHECK_FIELD(d,"density");
    CHECK_FIELD(vx,"velox");
    CHECK_FIELD(vy,"veloy");
    CHECK_FIELD(vz,"veloz");
    CHECK_FIELD(bx,"bfieldx");
    CHECK_FIELD(by,"bfieldy");
    CHECK_FIELD(bz,"bfieldz");
    CHECK_FIELD(pr,"pressure");

    /// Compute pressure
    // const int in = cello::index_static();
    // EnzoComputePressure compute_pressure (EnzoBlock::Gamma[in], false);
    // compute_pressure.compute(enzo_block, pr);

    FORTRAN_NAME(calc_dt_ppml_ig)
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
       bx, by, bz, pr,
       b0_, &gamma_,
       &dtBaryons);

    dtBaryons *= courant_;
    
  }
 
  /* 5) calculate minimum timestep */

  dt = std::numeric_limits<enzo_float>::max();

  dt = MIN(dt, dtBaryons);

  return dt;
}

//----------------------------------------------------------------------
