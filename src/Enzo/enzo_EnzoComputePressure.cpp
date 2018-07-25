// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputePressure.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoComputePressure class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputePressure::EnzoComputePressure (double gamma,
					  bool comoving_coordinates)
  : Compute(),
    gamma_(gamma),
    comoving_coordinates_(comoving_coordinates)
{
}

//----------------------------------------------------------------------

void EnzoComputePressure::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);

  p | gamma_;
  p | comoving_coordinates_;

}

//----------------------------------------------------------------------

void EnzoComputePressure::compute ( Block * block) throw()
{
  compute_(block);
}

//----------------------------------------------------------------------

void EnzoComputePressure::compute_(Block * block)
{

  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = enzo_block->data()->field();

  enzo_float * p = (enzo_float*) field.values("pressure");
  enzo_float * d = (enzo_float*) field.values("density");

  const int rank = enzo_block->rank();

  enzo_float * v3[3] =
    { (enzo_float*) (              field.values("velocity_x")),
      (enzo_float*) ((rank >= 2) ? field.values("velocity_y") : NULL),
      (enzo_float*) ((rank >= 3) ? field.values("velocity_z") : NULL) };

  enzo_float * te = (enzo_float*) field.values("total_energy");

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  int m = (nx+2*gx) * (ny+2*gy) * (nz+2*gz);
  enzo_float gm1 = gamma_ - 1.0;
  for (int i=0; i<m; i++) {
    enzo_float e= te[i];
    e -= 0.5*v3[0][i]*v3[0][i];
    if (rank >= 2) e -= 0.5*v3[1][i]*v3[1][i];
    if (rank >= 3) e -= 0.5*v3[2][i]*v3[2][i];
    p[i] = gm1 * d[i] * e;
  }


#ifdef CONFIG_USE_GRACKLE
  /* Correct pressure for significant amounts of H_2 */
  if( grackle_data->primordial_chemistry > 1){
    const EnzoUnits * enzo_units = (const EnzoUnits*)
         block->simulation()->problem()->units();

    enzo_float * HI_density    = (enzo_float*) field.values("HI_density");
    enzo_float * HII_density   = (enzo_float*) field.values("HII_density");
    enzo_float * HeI_density   = (enzo_float*) field.values("HeI_density");
    enzo_float * HeII_density  = (enzo_float*) field.values("HeII_density");
    enzo_float * HeIII_density = (enzo_float*) field.values("HeIII_density");
    enzo_float * H2I_density   = (enzo_float*) field.values("H2I_density");
    enzo_float * H2II_density  = (enzo_float*) field.values("H2II_density");
    enzo_float * e_density     = (enzo_float*) field.values("e_density");

    enzo_float temperature_units = cello::mass_hydrogen*
                             (enzo_units->velocity() * enzo_units->velocity())
                            / cello::kboltz;

    for(int i=0; i<m; i++){
      enzo_float number_density
        = 0.25*(HeI_density[i] + HeII_density[i] + HeIII_density[i]) +
           HI_density[i] + HII_density[i] + e_density[i];

      enzo_float nH2 = 0.5*(H2I_density[i] + H2II_density[i]);

      enzo_float temperature = std::max(temperature_units * p[i] / (number_density+nH2),1.0);

      enzo_float GammaH2Inverse = 0.5 * 5.0;
			enzo_float GammaInverse   = 1.0 / ( grackle_data->Gamma - 1.0);
      if (nH2 / GammaH2Inverse > 1.0E-3){
        enzo_float x = 6100.0 / temperature;
        if (x < 10.0){
          GammaH2Inverse = 0.5*(5.0 + 2.0 *x*x *exp(x)/( (exp(x)-1)*(exp(x)-1)));
        }
      }
      enzo_float Gamma1 = 1.0 + (nH2 + number_density) /
                                (nH2 + GammaH2Inverse + number_density *GammaInverse);

      // correct pressure
      p[i] *= (Gamma1 - 1.0) / (grackle_data->Gamma - 1.0);
    }
  }
#endif
}
