/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoMethodAccretionRemoveGas.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date
/// @brief  Implements a Method class for removing accreted gas from the gas
///         density field
///

#include "cello.hpp"
#include "enzo.hpp"


EnzoMethodAccretionRemoveGas::EnzoMethodAccretionRemoveGas()
  : Method()
{

  // Refresh operation accumulates values of density_accreted from neighbouring
  // blocks into density_accreted accumulate. 
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  ParticleDescr * particle_descr = cello::particle_descr();
  refresh->set_accumulate(true);
  refresh->add_field_src_dst("density_accreted","density_accreted_accumulate");
}

void EnzoMethodAccretionRemoveGas::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
 
  return;
}


void EnzoMethodAccretionRemoveGas::compute ( Block *block) throw()
{
  if (enzo::simulation()->cycle() == enzo::config()->initial_cycle){
    
    // Check if accretion_compute method precedes accretion_remove_gas
    // method
    ASSERT("EnzoMethodAccretionRemoveGas",
	   "accretion_compute must precede accretion_remove_gas",
	   enzo::problem()->method_precedes("accretion_compute",
					    "accretion_remove_gas"));
  }

  if (block->is_leaf()){
    this->compute_(block);
  }
  block->compute_done();

  return;
}

// Required
double EnzoMethodAccretionRemoveGas::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodAccretionRemoveGas::compute_(Block * block)
{

  // Get pointers to field data
  Field field = block->data()->field();

  const int id   = field.field_id("density");
  const int ida  = field.field_id("density_accreted");
  const int idaa = field.field_id("density_accreted_accumulate");

  int nx,ny,nz;
  int gx,gy,gz;
  field.dimensions (0,&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  int m = mx * my * mz;

  enzo_float * density = (enzo_float*) field.values(id);
  enzo_float * density_accreted = (enzo_float*) field.values(ida);
  enzo_float * density_accreted_accumulate = (enzo_float*) field.values(idaa);

  // Iterate over active cells. density_accreted and density_accreted_accumulate
  // have negative values, and are added to density to remove gas which has been
  // accreted onto a star particle

  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){

	int i = INDEX(ix,iy,iz,mx,my);
	density[i] += density_accreted[i];
	density[i] += density_accreted_accumulate[i];
      }
    }
  }
	 
  // Now set density_accreted and density_accreted_accumulate to zero in all cells,
  // including in ghost zones

  for (int i = 0; i < m; i++){
    density_accreted[i] = 0.0;
    density_accreted_accumulate[i] = 0.0;
  }
  return;
}
