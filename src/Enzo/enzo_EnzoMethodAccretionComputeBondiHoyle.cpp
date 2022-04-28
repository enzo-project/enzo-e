/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretionComputeBondiHoyle.cpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @author     John Regan (john.regan@mu.ie)
/// @date
/// @brief      Computes accretion rates according to Bondi-Hoyle model.
///             See Krumholz+ 2004, ApJ, 611, 399 for details.
///

#include "cello.hpp"
#include "enzo.hpp"

//------------------------------------------------------------------

EnzoMethodAccretionComputeBondiHoyle::EnzoMethodAccretionComputeBondiHoyle
(double accretion_radius_cells,
 double density_threshold,
 double max_mass_fraction,
 bool conserve_angular_momentum)
  : EnzoMethodAccretionCompute(accretion_radius_cells,
			       density_threshold,
			       max_mass_fraction,
			       conserve_angular_momentum)
{

}

//-------------------------------------------------------------------

void EnzoMethodAccretionComputeBondiHoyle::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodAccretionCompute::pup(p); // call parent class pup

  return;
}

//--------------------------------------------------------------------

void EnzoMethodAccretionComputeBondiHoyle::compute (Block * block) throw()
{
  if (enzo::simulation()->cycle() == enzo::config()->initial_cycle)
    do_checks_();

  // Only call compute_ if block is at highest refinement level.
  // Currently this method can only be used if refinement is turned off
  // (unigrid mode), but in future, we will have a refinement condition
  // which forces blocks containing accreting sink particles to be
  // at the highest refinement level, and using this method will
  // require this refinement condition to be activated.
  //if (block->level() == enzo::config()->mesh_max_level) {
  // this->compute_(block);
  //}

  block->compute_done();
  return;
}

//-------------------------------------------------------------------------

// void EnzoMethodAccretionComputeBondiHoyle::compute_(Block * block)

// {
//   // Get pointers to field data
//   Field field = block->data()->field();
//   const int id   = field.field_id("density");
//   const int ida  = field.field_id("density_accreted");
//   const int idvx  = field.field_id("velocity_x");
//   const int idvy  = field.field_id("velocity_y");
//   const int idvz  = field.field_id("velocity_z");
//   const int idie  = field.field_id("internal_energy");

//   enzo_float * density = (enzo_float*) field.values(id);
//   enzo_float * density_accreted = (enzo_float*) field.values(ida);
//   enzo_float * vx = (enzo_float*) field.values(idvx);
//   enzo_float * vy = (enzo_float*) field.values(idvy);
//   enzo_float * vz = (enzo_float*) field.values(idvz);
//   enzo_float * ie = (enzo_float*) field.values(idie);


//   // Get field dimensions
//   int mx,my,mz;
//   field.dimensions (0,&mx,&my,&mz);

//   // Get cell widths and cell volume
//   double hx, hy, hz;
//   block->cell_width(&hx, &hy, &hz);
//   const double cell_volume = hx * hy * hz;
//   // Also need these quantities in cgs units
//   const double hx_cgs = hx * enzo::units()->length();
//   const double hy_cgs = hy * enzo::units()->length();
//   const double hz_cgs = hz * enzo::units()->length();
//   const double cell_volume_cgs = hx_cgs * hy_cgs * hz_cgs;

//   Particle particle = block->data()->particle();
//   int it = particle.type_index("sink");
//   int num_particles = particle.num_particles(it);

//   if (num_particles > 0) {

//     // Declare pointers to particle attributes
//     enzo_float *pmass, *px, *py, *pz, *pvx, *pvy, *pvz;

//     // Get attribute indices
//     const int ia_m   = particle.attribute_index (it, "mass");
//     const int ia_x   = particle.attribute_index (it, "x");
//     const int ia_y   = particle.attribute_index (it, "y");
//     const int ia_z   = particle.attribute_index (it, "z");
//     const int ia_vx  = particle.attribute_index (it, "vx");
//     const int ia_vy  = particle.attribute_index (it, "vy");
//     const int ia_vz  = particle.attribute_index (it, "vz");
//     const int ia_acc = particle.attribute_index (it, "accretion_rate");

//     // Attribrute stride lengths
//     const int dm   = particle.stride(it, ia_m);
//     const int dp   = particle.stride(it, ia_x);
//     const int dv   = particle.stride(it, ia_vx);
//     const int dacc = particle.stride(it, ia_acc);

//     // Loop over batches
//     const int nb = particle.num_batches(it);
//     for (int ib=0; ib < nb; ib++){

//       pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
//       px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
//       py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
//       pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
//       pvx   = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
//       pvy   = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
//       pvz   = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
//       pacc  = (enzo_float *) particle.attribute_array(it, ia_acc, ib);

//       // Loop over particles in this batch
//       const int np = particle.num_particles(it,ib);
//       for (int ip=0; ip < np; ip++){

//	// Might want to have a check here that this is an accreting particle

// 	// Get the indices of the cell containing the particle
// 	const std::vector<int> host_cell_indices =
// 	  get_host_cell_indices_(block, px[ip*dp], py[ip*dp], pz[ip*dp]);

// 	// Get the indices of cells in the accretion zone, as well as the
// 	// squares of the distances from the particle
// 	std::pair<const std::vector<int>,const std::vector<double>> acc_zone =
// 	  get_accretion_zone_(block, px[ip*dp], py[ip*dp], pz[ip*dp]);

// 	const std::vector<int> acc_zone_indices = acc_zone.first;
// 	const std::vector<int> acc_zone_r2      = acc_zone.second;

// 	// Get the Bondi-Hoyle radius
// 	const double v_inf_2 =
// 	  compute_v_inf_2_(block, pvx[ip*dv], pvy[ip*dv], pvz[ip*dv], host_cell_indices);

// 	const double c_s_inf_2 = compute_c_s_inf_2_(block, host_cell_indices);

// 	const double r_BH_cgs =
// 	  compute_bondi_hoyle_radius_(pmass[ip*dm], v_inf_2, c_s_inf_2);

// 	// Get normalised weights for the cells in accretion zone
// 	const std::vector<double> weights =
// 	  compute_weights_(block, acc_zone_r2, r_BH_cgs);

// 	// For each cell in accretion zone, set new value for density_accreted field
// 	// in each cell, and return the amount of mass removed from each cell
// 	const std::vector<enzo_float> mass_removed =
// 	  compute_mass_removed_(block, weights, acc_zone_indices,
// 				 host_cell_indices, r_BH_cgs,
// 				 pmass[ip*dm] * enzo::units()->mass());

// 	// Also compute the total momentum removed (this "vector" is an actual vector!)
// 	const std::vector<enzo_float> momentum_removed =
// 	  compute_momentum_removed_(block, mass_removed);

// 	// Add the total removed mass and momentum to the particle
// 	enzo_float total_mass_removed = 0.0;
// 	for (auto m : mass_removed) total_mass_removed += m;

// 	const double old_particle_momentum_x = pmass[ip*dm] * pvx[ip*dv];
// 	const double old_particle_momentum_y = pmass[ip*dm] * pvy[ip*dv];
// 	const double old_particle_momentum_z = pmass[ip*dm] * pvz[ip*dv];

// 	const double new_particle_momentum_x =
// 	  old_particle_momentum_x + momentum_removed[0];
// 	const double new_particle_momentum_y =
// 	  old_particle_momentum_y + momentum_removed[1];
// 	const double new_particle_momentum_z =
// 	  old_particle_momentum_z + momentum_removed[2];

// 	pmass[ip*dm] += total_mass_removed;
// 	pvx[ip*dv] = new_particle_momentum_x / pmass[ip*dm];
// 	pvy[ip*dv] = new_particle_momentum_y / pmass[ip*dm];
// 	pvz[ip*dv] = new_particle_momentum_z / pmass[ip*dm];



//       } // Loop over particles in this batch
//     } // Loop over batches
//   } // if (num_particles > 0)
//   return;
// }

// // ------------------------------------------------------------------------------------

// const std::vector<int>
// EnzoMethodAccretionComputeBondiHoyle::get_host_cell_indices_
// (const Block * block, enzo_float x, enzo_float y, enzo_float z) throw()

// {
//   // Get field dimensions
//   int mx,my,mz;
//   block->data()->field().dimensions (0,&mx,&my,&mz);

//   // Get cell widths
//   double hx, hy, hz;
//   block->cell_width(&hx, &hy, &hz);

//   // Get the coordinates of the "front-lower-left" corner of the block
//   // (including ghost zones)
//   double xm, ym, zm;
//   block->data()->lower(&xm,&ym,&zm);
//   int gx, gy, gz;
//   field.ghost_depth(0,&gx,&gy,&gz);
//   double min_x = xm - gx * hx;
//   double min_y = ym - gy * hy;
//   double min_z = zm - gz * hz;

//   int ind_x = floor((x - min_x) / hx);
//   int ind_y = floor((y - min_y) / hy);
//   int ind_z = floor((z - min_z) / hz);

//   const std::vector<int> indices{ind_x,ind_y,ind_z};
//   return indices;
// }

// // -------------------------------------------------------------------------------------------

// std::pair<const std::vector<int>,const std::vector<double>>
// EnzoMethodAccretionComputeBondiHoyle::get_accretion_zone_
// (const Block * block, enzo_float x, enzo_float y, enzo_float z) throw()

// {
//   // Get field dimensions
//   int mx,my,mz;
//   block->data()->field().dimensions (0,&mx,&my,&mz);

//   // Get cell widths
//   double hx, hy, hz;
//   block->cell_width(&hx, &hy, &hz);

//   // Set accretion radius to be accretion_radius_cells_ multipled by minimum cell width
//   const double min_cell_width = std::min(hx,std::min(hy,hz));
//   const double accretion_radius = accretion_radius_cells_ * min_cell_width;

//   // Get the coordinates of the center of the cell with index 0 (outermost ghost zone cell)
//   double xm, ym, zm;
//   block->data()->lower(&xm,&ym,&zm);
//   int gx, gy, gz;
//   field.ghost_depth(0,&gx,&gy,&gz);
//   double first_cell_center_x = (0.5 - gx) * hx + xm;
//   double first_cell_center_y = (0.5 - gy) * hy + ym;
//   double first_cell_center_z = (0.5 - gz) * hz + zm;

//   // Find indices which specify a region which bounds the accretion zone
//   const int min_ind_x =
//     ceil((x - first_cell_center_x - accretion_radius) / hx);
//   const int min_ind_y =
//     ceil((y - first_cell_center_y - accretion_radius) / hy);
//   const int min_ind_z =
//     ceil((z - first_cell_center_z - accretion_radius) / hz);
//   const int max_ind_x =
//     floor((x - first_cell_center_x + accretion_radius) / hx);
//   const int max_ind_y =
//     floor((y - first_cell_center_y + accretion_radius) / hy);
//   const int max_ind_z =
//     floor((z - first_cell_center_z + accretion_radius) / hz);

//   std::vector<int> indices;
//   std::vector<double> r2_vec;


//   // Loop over all cells in this region
//   for (int k = min_ind_z; k <= max_ind_z; k++){
//     for (int j = min_ind_y; j <= max_ind_y; j++){
//       for (int i = min_ind_x; i <= max_ind_x; i++){

// 	// Check if center of cell is within accretion radius of particle
// 	const double cell_center_x = first_cell_center_x + i * hx;
// 	const double cell_center_y = first_cell_center_y + j * hy;
// 	const double cell_center_z = first_cell_center_z + k * hz;

// 	const double r2 =
// 	  (px[ip*dp] - cell_center_x) * (px[ip*dp] - cell_center_x)
// 	  + (py[ip*dp] - cell_center_y) * (py[ip*dp] - cell_center_y)
// 	  + (pz[ip*dp] - cell_center_z) * (pz[ip*dp] - cell_center_z);

// 	if (r2 < accretion_radius * accretion_radius){
// 	  const int index = INDEX(i,j,k,mx,my);
// 	  indices.push_back(index);
// 	  r2_vec.push_back(r2);
// 	} // if (r2 < accretion_radius * accretion_radius)      
//       }
//     }
//   } // Loop over cells in region bounding accretion zone
//   return indices;
// }

// // ----------------------------------------------------------------------------------

// const double compute_bondi_hoyle_radius_
// (const double pmass, const double v_inf_2, const double c_s_2) throw()
// {

//   // Convert quantities into cgs units
//   const double v_inf_cgs_2 = v_inf_2 *
//     enzo::units()->velocity() * enzo::units()->velocity();

//   const double c_s_inf_cgs_2 = c_s_inf_2 *
//     enzo::units()->velocity() * enzo::units()->velocity();

//   const double pmass_cgs = pmass * enzo::units()->mass();

//   // Compute the Bondi-Hoyle radius in cgs units
//   // (Equation 10 in Krumholz paper)
//   const double r_BH_cgs = cello::grav_constant * pmass_cgs /
//     (v_inf_cgs_2 + c_s_inf_cgs_2);

//   // Return the Bondi-Hoyle radius in code units
//   return r_BH_cgs / enzo::units()->length();
// }

// // ----------------------------------------------------------------------------------

// const double compute_v_inf_2_
// (const Block * block, const double pvx, const double pvy,
//  const double pvz, const std::vector<int> cell_indices) throw() {

//   Field field = block->data()->field();
//   enzo_float * vx = (enzo_float*) field.values("velocity_x");
//   enzo_float * vy = (enzo_float*) field.values("velocity_y");
//   enzo_float * vz = (enzo_float*) field.values("velocity_z");

//   int mx,my,mz;
//   field.dimensions (0,&mx,&my,&mz);
//   const int cell_index =
//     INDEX(cell_indices[0],cell_indices[1],cell_indices[2],mx,my);

//   return (vx[cell_index] - pvx) * (vx[cell_index] - pvx) +
//     (vy[cell_index] - pvy) * (vy[cell_index] - pvy) +
//     (vz[cell_index] - pvz) * (vz[cell_index] - pvz);

// }

// // ----------------------------------------------------------------------------------

// const double compute_c_s_inf_2_
// (const Block * block, const std::vector<int> cell_indices) throw() {

//   Field field = block->data()->field();
//   enzo_float * ie = (enzo_float*) field.values("internal_energy");

//   int mx,my,mz;
//   field.dimensions (0,&mx,&my,&mz);
//   const int cell_index =
//     INDEX(cell_indices[0],cell_indices[1],cell_indices[2],mx,my);

//   // Now compute the square of the sound speed
//   const double gamma = enzo::config()->field_gamma;
//   return  gamma * (gamma - 1.0) * internal_energy[cell_index];

// }

// // ----------------------------------------------------------------------------------

// const std::vector<double> EnzoMethodAccretionComputeBondiHoyle::compute_weights_
// (const Block * block,
//  const std::vector<double>& r2_vec,
//  const double r_BH) throw()
// {
//   // Get min cell width and accretion radius
//   double hx, hy, hz;
//   block->cell_width(&hx, &hy, &hz);
//   const double min_cell_width = std::min(hx,std::min(hy,hz));
//   const double accretion_radius = accretion_radius_cells_ * min_cell_width;

//   // Set the kernel radius (Equation 13 of Krumholz paper)
//   const double kernel_radius =
//     (r_BH < min_cell_width / 4.0) ? min_cell_width / 4.0 :
//     ((r_BH <= accretion_radius / 2.0) ? r_BH : accretion_radius / 2.0);

//   const double kernel_radius_2 = kernel_radius * kernel_radius;

//   // Compute the weights using the kernel radius
//   // (Equation 14 of Krumholz paper)
//   std::vector<double> weights;
//   double sum_of_weights = 0.0;

//   for (auto r2 : r2_vec){
//     const double w = exp(-r2 / kernel_radius_2);
//     weights.push_back(w);
//     sum_of_weights += w;
//   }

//   // Normalize the weights
//   for (auto &w : weights) w /= sum_of_weights;

//   return weights;
// }

// // ----------------------------------------------------------------------------------------

// const std::vector<enzo_float>
// EnzoMethodAccretionComputeBondiHoyle::compute_mass_removed_
// (Block * block,
//  const std::vector<enzo_float>& weights,
//  const std::vector<int>& acc_zone_indices,
//  const std::vector<int>& host_cell_indices,
//  const double r_BH,
//  const double pmass,
//  const double v_inf_2,
//  const double c_s_inf_2) throw()
// {

//   // Compute the weighted mean density
//   double weighted_mean_density = 0.0;
//   for (decltype(weights.size()) i = 0; i < weights.size(); ++i)
//     weighted_mean_density += weights[i] * density[acc_zone_indices[i]];

//   // Get min cell width
//   double hx, hy, hz;
//   block->cell_width(&hx, &hy, &hz);
//   const double min_cell_width = std::min(hx,std::min(hy,hz));

//   // Get Bondi-Hoyle radius in code units
//   const double r_BH = r_BH_cgs / enzo::units()->length();

//   // Compute rho_inf using weighted mean density and "alpha function"
//   // (Equation 12 of Krumholz paper)
//   const double rho_inf = weighted_mean_density /
//     alpha_(1.2 * min_cell_width / r_BH);

//   // Now compute accretion rate from Equation 11 of same paper
//   const double lambda_ = 0.25 * exp(1.5);
//   const double acc_rate =
//     4.0 * cello::pi * rho_inf * r_BH * r_BH *
//     sqrt(lambda_ * lambda_ * c_s_inf_2 + v_inf_2); 

//   // Not sure if timestep should be multiplied by some cosmological factor
//   const double accreted_mass = acc_rate * block->dt();
// }

// // -----------------------------------------------------------------------------------------

// const std::vector<enzo_float>
// EnzoMethodAccretionComputeBondiHoyle::compute_momentum_removed_
// (const Block * block, const std::vector<enzo_float>& mass_removed) throw();

// // ----------------------------------------------------------------------------------------


// /* Ported (copy and paste) direct from enzo-dev */
// double alpha_(double x) {

// #define XMIN 0.01
// #define XMAX 2.0
// #define NTABLE 51

//   double lambda_c, xtable, xtablep1, alpha_exp;
//   int idx;

//   /* This is a precomputed table of alpha values.  These correspond to x values
//      that run from 0.01 to 2.0 with uniform logarithmic spacing.  The reason for
//      this choice of range is that the asymptotic expressions are accurate to
//      better than 2% outside this range */

//   double alphatable[NTABLE] = {820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427,
// 			      277.295, 237.845, 204.1, 175.23, 150.524, 129.377, 111.27, 95.7613,
// 			      82.4745, 71.0869, 61.3237, 52.9498, 45.7644, 39.5963, 34.2989,
// 			      29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, 13.0714,
// 			      11.4717, 10.0903, 8.89675, 7.86467, 6.97159, 6.19825, 5.52812,
// 			      4.94699, 4.44279, 4.00497, 3.6246, 3.29395, 3.00637, 2.75612,
// 			      2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, 1.70804,
// 			      1.62439};

//   // A constant that appears in the following formulae.  This hardcoded value is
//   // valid for an isothermal gas.
//   lambda_c = 0.25*exp(1.5);

//   // deal with the off-the-table cases
//   if (x < XMIN)
//     return lambda_c / sqrt(2.*x*x);
//   else if (x >= XMAX)
//     return exp(1./x);
//   else {
//     // we are on the table

//     idx = floor((NTABLE-1)*log(x/XMIN)/log(XMAX/XMIN));
//     xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1));
//     xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1));
//     alpha_exp = log(x/xtable) / log(xtablep1/xtable);

//     return alphatable[idx] * pow(alphatable[idx+1]/alphatable[idx],alpha_exp);
//   }

// #undef NTABLE
// #undef XMIN
// #undef XMAX

// }
