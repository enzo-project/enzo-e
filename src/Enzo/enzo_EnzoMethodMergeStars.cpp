/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoMethodMergeStars.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   18 January 2022
/// @brief  Implementation of EnzoMethodMergeStars class, a method for merging
///         star particles separated by a distance less than the merging
///         radius.
///
///         After copying particles from
///         neighbouring blocks during the refresh, this method then uses an
///         iterative procedure to merge particles together, so that no
///         two particles are within a 'merging radius' of each other.

#include "cello.hpp"
#include "enzo.hpp"
#include "FofLib.hpp"
#include <time.h>

//#define DEBUG_MERGESTARS


EnzoMethodMergeStars::EnzoMethodMergeStars(double merging_radius_cells)
  : Method(),
    merging_radius_cells_(merging_radius_cells)
{
  // This method requires three dimensions.
  ASSERT("EnzoMethodMergeStars::EnzoMethodMergeStars()",
	 "EnzoMethodMergeStars requires that we run a 3D problem (Domain: rank = 3)",
	 cello::rank());
  
  const EnzoConfig * enzo_config = enzo::config();
  ASSERT("EnzoMethodMergeStars::EnzoMethodMergeStars()",
	 "EnzoMethodMergeStars requires unigrid mode (Adapt : max_level = 0). "
	 "In future, we may put in a refinement condition that blocks containing "
	 "star particles or neighbouring such a block is at highest refinement "
	 "level", enzo_config->mesh_max_level == 0);

  // Refresh copies all star particles from neighbouring blocks
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  ParticleDescr * particle_descr = cello::particle_descr();
  refresh->add_particle(particle_descr->type_index("star"));
  refresh->set_particles_are_copied(true);
}

void EnzoMethodMergeStars::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | merging_radius_cells_;
 
  return;
}


void EnzoMethodMergeStars::compute ( Block *block) throw()
{

  
  if (block->is_leaf()){
    this->compute_(block);
  }
  block->compute_done();

  return;
}

// Required
double EnzoMethodMergeStars::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodMergeStars::compute_(Block * block)
{
  Hierarchy * hierarchy = cello::hierarchy();  
  EnzoBlock * enzo_block = enzo::block(block);
  
  // Get width of block in x-direction, used to mark particles for deletion
  double block_xm, block_ym, block_zm, block_xp, block_yp, block_zp;
  enzo_block->lower(&block_xm,&block_ym,&block_zm);
  enzo_block->upper(&block_xp,&block_yp,&block_zp);
  const double block_width_x = block_xp - block_xm;
  
  Particle particle = enzo_block->data()->particle();
  int it = particle.type_index("star");
  int num_particles = particle.num_particles(it);
  
#ifdef DEBUG_MERGESTARS
  CkPrintf("In total, there are %d particles on Block %s \n",num_particles,
	   block->name().c_str());
#endif

  // Declare pointers to particle attributes
  enzo_float *px, *py, *pz, *pvx, *pvy, *pvz;
  enzo_float *plifetime, *pcreation, *pmass, *pmetal;
  int64_t *pid;
  
    // Get attribute indices
  const int ia_m   = particle.attribute_index (it, "mass");
  const int ia_x   = particle.attribute_index (it, "x");
  const int ia_y   = particle.attribute_index (it, "y");
  const int ia_z   = particle.attribute_index (it, "z");
  const int ia_vx  = particle.attribute_index (it, "vx");
  const int ia_vy  = particle.attribute_index (it, "vy");
  const int ia_vz  = particle.attribute_index (it, "vz");
  const int ia_l   = particle.attribute_index (it, "lifetime");
  const int ia_c   = particle.attribute_index (it, "creation_time");
  const int ia_mf  = particle.attribute_index (it, "metal_fraction");
  const int ia_id  = particle.attribute_index(it,"id");
  
  // Attribrute stride lengths
  const int dm   = particle.stride(it, ia_m);
  const int dp   = particle.stride(it, ia_x);
  const int dv   = particle.stride(it, ia_vx);
  const int dl   = particle.stride(it, ia_l);
  const int dc   = particle.stride(it, ia_c);
  const int dmf  = particle.stride(it, ia_mf);
  const int did  = particle.stride(it, ia_id);
      
  // Array giving the FoF group number of each particle
  int * group_index = new int[num_particles];
  
  // group_sizes will be an array giving the number of particles in each
  // group
  int *  group_size;
  
  // group_lists will be an 'array of arrays'. Each element will be an
  // array containing the indices of particles belonging to a particular
  // group
  int ** group_lists;
  
  // Array containing particle positions in 'block units'
  enzo_float * particle_coordinates = new enzo_float[3 * num_particles];
      
  // Fill in particle coordinates array.
  // This handles periodic boundary conditions, by taking the periodic
  // image of particle positions if necessary.
  
  get_particle_coordinates_(enzo_block, it, particle_coordinates);

  // Get the max cell width (across three dimensional axes), used to calculate
  // the merging radius
  
  double cell_width_x, cell_width_y,cell_width_z;
  enzo_block->cell_width(&cell_width_x,&cell_width_y,&cell_width_z);
  const double max_cell_width = std::max(
				std::max(cell_width_x,cell_width_y),
				cell_width_z);
  const enzo_float merging_radius = merging_radius_cells_ * max_cell_width;

  // FofList runs the Friends-of-Friends algorithm on particle positions
  // (given by the particle_coordinates_block_units array), with
  // merging_radius_block_units as the linking length. This function
  // fills in the (already allocated) group_index array, and allocates
  // and fills in the group_size and group_lists arrays.

  int ngroups = FofList(num_particles, particle_coordinates,merging_radius, 
			group_index, &group_size, &group_lists);
  
#ifdef DEBUG_MERGESTARS
  CkPrintf("The %d particles on Block %s are in %d FoF groups \n",num_particles,
	   block->name().c_str(),ngroups);
#endif 
  
  for (int i = 0; i < ngroups; i++){
    
#ifdef DEBUG_MERGESTARS
    CkPrintf("Group %d out of %d on block %s: Group size = %d \n",i+1, ngroups,
	     block->name().c_str(),group_size[i]);
#endif

    // Only need to merge particles if there are two or more particles in the
    // group 
    if (group_size[i] > 1){
      
      ASSERT("EnzoMethodMergeStars::compute_()",
	     "There is a FoF group containing a pair of star particles "
	     "in non-neighbouring blocks. Since this cannot be properly "
	     "dealt with we exit the program here. This has likely "
	     "happened because the merging radius is too large in "
	     "comparison to the block size.",
	     particles_in_neighbouring_blocks_(enzo_block,particle_coordinates,
					       group_lists,group_size,i));
      
      // ib1 and ip1 index the first particle in this group
      int ib1, ip1;
      particle.index(group_lists[i][0],&ib1,&ip1);
      
      // We get the attribrutes of this particle and store them in a new
      // variable
      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib1);
      px = (enzo_float *) particle.attribute_array(it, ia_x, ib1);
      py = (enzo_float *) particle.attribute_array(it, ia_y, ib1);
      pz = (enzo_float *) particle.attribute_array(it, ia_z, ib1);
      pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib1);
      pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib1);
      pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib1);
      plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib1);
      pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib1);
      pmetal    = (enzo_float *) particle.attribute_array(it, ia_mf, ib1);
      pid = (int64_t *) particle.attribute_array(it, ia_id, ib1);
      
      enzo_float pmass1 = pmass[ip1*dm];
      enzo_float px1 = px[ip1*dp];
      enzo_float py1 = py[ip1*dp];
      enzo_float pz1 = pz[ip1*dp];
      double     pos1[3] = {px1,py1,pz1};
      enzo_float pvx1 = pvx[ip1*dv];
      enzo_float pvy1 = pvy[ip1*dv];
      enzo_float pvz1 = pvz[ip1*dv];
      enzo_float plifetime1 = plifetime[ip1*dl];
      enzo_float pcreation1 = pcreation[ip1*dc];
      enzo_float pmetal1 = pmetal[ip1*dmf];
      int64_t    pid1    = pid[ip1*did];
      
      // now loop over the rest of the particles in this group, and merge
      // them in to the first particle
      
      for (int j = 1; j < group_size[i]; j++){

	// ib2 and ip2 are used to index the other particles in this group 
	int ib2, ip2;
	particle.index(group_lists[i][j],&ib2,&ip2);
	
	// get attributes of this particle
	pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib2);
	px = (enzo_float *) particle.attribute_array(it, ia_x, ib2);
	py = (enzo_float *) particle.attribute_array(it, ia_y, ib2);
	pz = (enzo_float *) particle.attribute_array(it, ia_z, ib2);
	pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib2);
	pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib2);
	pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib2);
	plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib2);
	pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib2);
	pmetal    = (enzo_float *) particle.attribute_array(it, ia_mf, ib2);
	pid    = (int64_t *) particle.attribute_array(it, ia_id, ib2);
		
	enzo_float pmass2 = pmass[ip2*dm];
	enzo_float px2 = px[ip2*dp];
	enzo_float py2 = py[ip2*dp];
	enzo_float pz2 = pz[ip2*dp];
	enzo_float pvx2 = pvx[ip2*dv];
	enzo_float pvy2 = pvy[ip2*dv];
	enzo_float pvz2 = pvz[ip2*dv];
	double     pos2[3] = {px2,py2,pz2};
	enzo_float plifetime2 = plifetime[ip2*dl];
	enzo_float pcreation2 = pcreation[ip2*dc];
	enzo_float pmetal2 = pmetal[ip2*dmf];
	int64_t    pid2 = pid[ip2*did];
	
	enzo_float f1 = pmass1 / (pmass1 + pmass2);
	enzo_float f2 = 1.0 - f1;
	
#ifdef DEBUG_MERGESTARS
	CkPrintf("Merger in Group %d: \n Particle 1: Mass = %g, "
		 "Position = (%g,%g,%g), ID = %ld \n"
		 "Particle 2: Mass = %g, "
		 "Position = (%g,%g,%g), ID = %ld \n", i+1, pmass1,
		 px1,py1,pz1,pid1,pmass2,px2,py2,pz2,pid2);
	
#endif
	// Get the nearest periodic image of particle 2 to particle 1
	double npi[3];
	hierarchy->get_nearest_periodic_image(pos2,pos1,npi);
	
	// Compute new properties of 'particle 1'
	px1 = f1 * px1 + f2 * npi[0];
	py1 = f1 * py1 + f2 * npi[1];
	pz1 = f1 * pz1 + f2 * npi[2];
	pvx1 = f1 * pvx1 + f2 * pvx2;
	pvy1 = f1 * pvy1 + f2 * pvy2;
	pvz1 = f1 * pvz1 + f2 * pvz2;
	plifetime1 = std::min(plifetime1, plifetime2);
	pcreation1 = std::min(pcreation1, pcreation2);
	pmetal1 = f1 * pmetal1 + f2 * pmetal2;
	pid1 = std::min(pid1,pid2);
	pmass1 += pmass2;
	
#ifdef DEBUG_MERGESTARS
	CkPrintf("Particle number %d in group %d out of %d on Block %s "
		 "is merged into 0th particle. New properties of 0th particle: "
		 "Mass = %g. Position = (%g,%g,%g)\n",
		 j,i,ngroups, block->name().c_str(),pmass1,px1,py1,pz1);
#endif


	// To mark this particle for deletion after it has been merged into
	// another particle, we shift its position out of the block, so that
	// when we call particle->delete_non_local_particles_, this particle
	// is deleted. Yes, I know its a bit of a hack

	px[ip2*dp] += 2.0 * block_width_x;
	
      } // Loop over particles in group
      
      // Set new properties for the first particle in the group, into which
      // the other particles have been merged.
      
      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib1);
      px = (enzo_float *) particle.attribute_array(it, ia_x, ib1);
      py = (enzo_float *) particle.attribute_array(it, ia_y, ib1);
      pz = (enzo_float *) particle.attribute_array(it, ia_z, ib1);
      pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib1);
      pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib1);
      pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib1);
      plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib1);
      pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib1);
      pmetal    = (enzo_float *) particle.attribute_array(it, ia_mf, ib1);
      pid    = (int64_t *) particle.attribute_array(it, ia_id, ib1);
      
      pmass[ip1*dm] = pmass1;
      
      double folded_pos[3];
      pos1[0] = px1;
      pos1[1] = py1;
      pos1[2] = pz1;
      // Fold position within the domain if periodic boundary positions
      hierarchy->get_folded_position(pos1,folded_pos);
      px[ip1*dp] = folded_pos[0];
      py[ip1*dp] = folded_pos[1];
      pz[ip1*dp] = folded_pos[2];      
      pvx[ip1*dv] = pvx1;
      pvy[ip1*dv] = pvy1;
      pvz[ip1*dv] = pvz1;
      plifetime[ip1*dl] = plifetime1;
      pcreation[ip1*dc] = pcreation1;
      pmetal[ip1*dmf] = pmetal1;
      pid[ip1*did] = pid1;
      
    }// if (group_size[i] > 1)
      
  }// Loop over Fof groups
  
  // Delete the dynamically allocated arrays
  
  delete [] group_index;
  free(group_size);
  free(group_lists);
  delete [] particle_coordinates;
  
#ifdef DEBUG_MERGESTARS
  CkPrintf("Block %s: After merging, num_particles = %d \n",
	   block->name().c_str(),particle.num_particles(it));
#endif
  
  // Now we delete non-local particles, i.e. particles with positions out-of-bounds
  // of the block.
  int delete_count = enzo_block->delete_non_local_particles_(it);
  cello::simulation()->data_delete_particles(delete_count);
      
#ifdef DEBUG_MERGESTARS
  CkPrintf("Block %s: After deletion, num_particles = %d \n",
	   block->name().c_str(),particle.num_particles(it));
#endif
  
  return;
      
}

// This fills a 1D array, which must have already been allocated with length
// 3 * num_particles, with x, y, z coordinates of star particles in the block.
// This also takes care of periodic boundary conditions by taking the nearest
// periodic image of coordinates if necessary.

void EnzoMethodMergeStars::get_particle_coordinates_
  (EnzoBlock * enzo_block, int it,
   enzo_float * particle_coordinates)
{
  Hierarchy * hierarchy = cello::hierarchy();
  
  // Get coordinates of the centre of the block
  double block_xm, block_ym, block_zm, block_xp, block_yp, block_zp;
  enzo_block->lower(&block_xm,&block_ym,&block_zm);
  enzo_block->upper(&block_xp,&block_yp,&block_zp);
  const double block_centre_x = 0.5 * (block_xm + block_xp);
  const double block_centre_y = 0.5 * (block_ym + block_yp);
  const double block_centre_z = 0.5 * (block_zm + block_zp);
  double block_centre[3] = {block_centre_x,block_centre_y,block_centre_z};
  
  Particle particle = enzo_block->data()->particle();
  const int ia_x = particle.attribute_index (it, "x");
  const int ia_y = particle.attribute_index (it, "y");
  const int ia_z = particle.attribute_index (it, "z");
  const int dp = particle.stride(it, ia_x);
  enzo_float *px, *py, *pz;

  // ip_block is particle index within the block
  // ib is batch index
  // ip_batch is particle index within a batch
  int ip_block, ib, ip_batch;

  const int num_particles = particle.num_particles(it);

  // Loop over all particles in block
  for (ip_block = 0; ip_block < num_particles; ip_block++){

    // Get the particle's batch index and its index within the batch
    particle.index(ip_block,&ib,&ip_batch);
    
    // Get pointers to the attribute arrays
    px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

    // Get the nearest periodic image to the block centre. If
    // boundary conditions are non-periodic, this just returns the particle
    // coordinates
    double npi[3];
    double pos[3] = {px[ip_batch*dp],py[ip_batch*dp],pz[ip_batch*dp]};
    hierarchy->get_nearest_periodic_image(pos,block_centre,npi);

    // Now we can set particle coordinates

    particle_coordinates[3*ip_block]     = npi[0];
    particle_coordinates[3*ip_block + 1] = npi[1];
    particle_coordinates[3*ip_block + 2] = npi[2];
    
  } // ip_block loop
  
  return;
}

// Checks if all the particles within a group (specified by group_index)
// are in neighbouring blocks
bool EnzoMethodMergeStars::particles_in_neighbouring_blocks_
(EnzoBlock * enzo_block,
 enzo_float * particle_coordinates,
 int ** group_lists, int * group_size,
 int group_index)
{
  bool return_val = 1;
  
  // Get block widths and block centre
  double block_xm, block_ym, block_zm, block_xp, block_yp, block_zp;
  enzo_block->lower(&block_xm,&block_ym,&block_zm);
  enzo_block->upper(&block_xp,&block_yp,&block_zp);
  const double block_width_x = block_xp - block_xm;
  const double block_width_y = block_yp - block_ym;
  const double block_width_z = block_zp - block_zm;
  
  // Loop over all particles, getting their positions in a block-centred
  // frame-of-reference, where the 'left' and 'right' faces of the block, in all
  // 3 dimensions, have coordinates 0 and 1 respectively. Checking if a particle
  // is in the block is equivalent to its x,y,z coordinates in this
  // frame-of-reference being between 0 and 1.
  for (int j = 0; j < group_size[group_index]; j++){
    const int ind_1 = group_lists[group_index][j];
    
    const enzo_float px1 =
      (particle_coordinates[3*ind_1]     - block_xm) / block_width_x;
    const enzo_float py1 =
      (particle_coordinates[3*ind_1 + 1] - block_ym) / block_width_y;
    const enzo_float pz1 =
      (particle_coordinates[3*ind_1 + 2] - block_zm) / block_width_z;

    // if particle is in bounds, then there is no problem, don't
    // need to check all the pairs containing this particle
    if (px1 > 0.0 && px1 < 1.0 &&
	py1 > 0.0 && py1 < 1.0 &&
	pz1 > 0.0 && pz1 < 1.0) continue;
    
    // Otherwise need to loop over all particles which have not already
    // been considered, checking if the pair (j,k) are on non-neighbouring
    // blocks.
    for (int k = j; k < group_size[group_index]; k++){
      const int ind_2 = group_lists[group_index][k];
      const enzo_float px2 =
	(particle_coordinates[3*ind_2]     - block_xm) / block_width_x;
      const enzo_float py2 =
	(particle_coordinates[3*ind_2 + 1] - block_ym) / block_width_y;
      const enzo_float pz2 =
	(particle_coordinates[3*ind_2 + 2] - block_zm) / block_width_z;

      // In each dimension, check if the coordinate of one of pair is less than 0
      // and the other greater than 1.
      if (px1 < 0.0 && px2 > 1.0 ||
	  px2 < 0.0 && px1 > 1.0 ||
	  py1 < 0.0 && py2 > 1.0 ||
	  py2 < 0.0 && py1 > 1.0 ||
	  pz1 < 0.0 && pz2 > 1.0 ||
	  pz2 < 0.0 && pz1 > 1.0) {
	return_val = 0;
	break; // break out of the k loop
      }
      
    } // k loop
    
    if (return_val == 0) break; // break out of the j loop

  } // j loop
  
  return return_val;
}
