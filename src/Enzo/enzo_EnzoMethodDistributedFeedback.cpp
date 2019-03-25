/// See LICENSE_CELLO file for license and copyright information

///
///
///
///

#include "cello.hpp"
#include "enzo.hpp"

EnzoMethodDistributedFeedback::EnzoMethodDistributedFeedback
()
  : Method()
{
  FieldDescr * field_descr = cello::field_descr();
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  // Initialize default refresh object
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
                             enzo_sync_id_method_feedback);
  refresh(ir)->add_all_fields();

  dual_energy_         = enzo_config->ppm_dual_energy;

  total_ejecta_mass_   = enzo_config->method_feedback_ejecta_mass * cello::mass_solar /
                         enzo_units->mass();

  total_ejecta_energy_ = enzo_config->method_feedback_supernova_energy * 1.0E51 /
                         enzo_units->mass() / enzo_units->velocity() /
                         enzo_units->velocity();

  kinetic_fraction_    = enzo_config->method_feedback_ke_fraction;

  ejecta_metal_fraction_ = enzo_config->method_feedback_ejecta_metal_fraction;

  stencil_                   = enzo_config->method_feedback_stencil;
  stencil_rad_               = ( (int) (stencil_ - 1 / 2.0));
  number_of_feedback_cells_  = stencil_ * stencil_ * stencil_;
  shift_cell_center_         = enzo_config->method_feedback_shift_cell_center;

  // Do error checking here to make sure all required
  // fields exist..
  // NOTE: Good idea - make ALL method objects have a 'required fields' list
  //       defined in the header file and then loop through this in a more general
  //       fashion at initialization of problem instead of cluttering every
  //       method object init with this

  return;
}

void EnzoMethodDistributedFeedback::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | total_ejecta_mass_;
  p | total_ejecta_energy_;
  p | ejecta_metal_fraction_;
  p | kinetic_fraction_;

  p | dual_energy_;
  p | stencil_;
  p | stencil_rad_;
  p | number_of_feedback_cells_;
  //p | mass_per_cell;
  //p | energy_per_cell;
  p | shift_cell_center_;

  return;
}

void EnzoMethodDistributedFeedback::compute (Block * block) throw()
{

  if (block->is_leaf()){
    this->compute_(block);
  }

  block->compute_done();

  return;
}

void EnzoMethodDistributedFeedback::compute_ (Block * block)
{

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  Particle particle = enzo_block->data()->particle();

  EnzoUnits * enzo_units = enzo::units();

  Field field = block->data()->field();
  // Obtain grid sizes and ghost sizes

  enzo_float * d           = (enzo_float *) field.values("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");

  enzo_float * v3[3]       = { (enzo_float *) field.values("velocity_x"),
                               (enzo_float *) field.values("velocity_y"),
                               (enzo_float *) field.values("velocity_z")   };

  enzo_float * metal = (enzo_float *) field.values("metal_density");

  int mx, my, mz, gx, gy, gz, nx, ny, nz;
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;

  // Block sizes (exlcuding ghost zones)
  field.size(&nx,&ny,&nz);

  field.ghost_depth(0,&gx,&gy,&gz);
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;

  // We will probably never be in the situation of constant acceleration
  // and cosmology, but just in case.....
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  enzo_float cosmo_a = 1.0;

  const int rank = cello::rank();
>>>>>>> 5899fa6868f2e4610229f4614ad5e1d2e93aa24e

  double current_time  = block->time();

  // apply feedback depending on particle type
  // for now, just do this for all star particles

  int it = particle.type_index("star");
  int count = 0;

  if (particle.num_particles(it) > 0 ){

    const int ia_m = particle.attribute_index (it, "mass");

    const int ia_x  = particle.attribute_index (it, "x");
    const int ia_y  = particle.attribute_index (it, "y");
    const int ia_z  = particle.attribute_index (it, "z");
    const int ia_vx = particle.attribute_index (it, "vx");
    const int ia_vy = particle.attribute_index (it, "vy");
    const int ia_vz = particle.attribute_index (it, "vz");

    const int ia_l = particle.attribute_index (it, "lifetime");
    const int ia_c = particle.attribute_index (it, "creation_time");

    const int dm = particle.stride(it, ia_m);
    const int dp = particle.stride(it, ia_x);
    const int dv = particle.stride(it, ia_vx);
    const int dl = particle.stride(it, ia_l);
    const int dc = particle.stride(it, ia_c);

    const int nb = particle.num_batches(it);

    for (int ib=0; ib<nb; ib++){
      enzo_float *px=0, *py=0, *pz=0, *pvx=0, *pvy=0, *pvz=0;
      enzo_float *plifetime=0, *pcreation=0, *pmass=0;

      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);

      px  = (enzo_float *) particle.attribute_array(it, ia_x, ib);
      py  = (enzo_float *) particle.attribute_array(it, ia_y, ib);
      pz  = (enzo_float *) particle.attribute_array(it, ia_z, ib);
      pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
      pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
      pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

      plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
      pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib);

      int np = particle.num_particles(it,ib);

      for (int ip=0; ip<np; ip++){
        // AE: Check and see if these differ....
        int ipdp = ip*dp;
        int ipdm = ip*dm;
        int ipdv = ip*dv;
        int ipdl = ip*dl;
        int ipdc = ip*dc;

        // negative lifetime are particles that have alreahy gone SN
        // creation time must be > 0
        // only go SN if age >= lifetime
        if ( (plifetime[ipdl] <= 0.0) || (pcreation[ipdc] <= -100.0) ||
             (current_time - pcreation[ipdc]) < plifetime[ipdl]) continue;
        count++;

        // Update particle properties here if needed
        plifetime[ipdl] *= -1.0;                  // set to negative - flag for alreahy gone SNe
        pmass[ipdm] = pmass[ipdm] - total_ejecta_mass_;

        // get corresponding grid position of particle
        // and shift it if it is too close to the grid boundaries
        double xpos = px[ipdp];
        double ypos = py[ipdp];
        double zpos = pz[ipdp];

        this->inject_feedback(block, xpos, ypos, zpos,
                              enzo_config->method_feedback_ejecta_mass,
                              enzo_config->method_feedback_supernova_energy,
                              enzo_config->method_feedback_ke_fraction);

        count++;
      } // end loop over particles
    } // end loop over batches

    if (count > 0){
      std::cout << "Number of feedback particles:   " << count << "\n";
    }

  } // end particle check

  return;
}

void EnzoMethodDistributedFeedback::inject_feedback(
                                          Block * block,
                                          double xpos, double ypos, double zpos,
                                          double m_eject, double E_51,
                                          double ke_fraction,
                                          enzo_float pvx,    //default -999
                                          enzo_float pvy,    //default -999
                                          enzo_float pvz){   //default -999
  /*

    Inject feedback with ejecta mass m_eject (in Msun) and total energy
    (in units of 10^51 erg) split between thermal and kinetic energy as
    determined by ke_fraction. Location of feedback is centered at
    xpos,ypos,zpos, but is shifted away from grid boundaries if
    shift_cell_center is ON.

    pvx,pvy, and pvz are the local velocity of the injection source. These are only
    needed if ke_fraction > 0. By default these are set to the local
    gas flow

  */

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  Field field = block->data()->field();
  // Obtain grid sizes and ghost sizes

  enzo_float * d           = (enzo_float *) field.values("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");

  enzo_float * v3[3]       = { (enzo_float *) field.values("velocity_x"),
                               (enzo_float *) field.values("velocity_y"),
                               (enzo_float *) field.values("velocity_z")   };

  enzo_float * metal = (enzo_float *) field.values("metal_density");

  int mx, my, mz, gx, gy, gz, nx, ny, nz;
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  nx = mx + 2*gx;
  ny = my + 2*gy;
  nz = mz + 2*gz;

  // We will probably never be in the situation of constant acceleration
  // and cosmology, but just in case.....
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  enzo_float cosmo_a = 1.0;

  const int rank = cello::rank();

  double current_time  = block->time();
  if (cosmology) {
    enzo_float cosmo_dadt = 0.0;
    double dt    = block->dt();
    cosmology->compute_expansion_factor(&cosmo_a,&cosmo_dadt,current_time+0.5*dt);
    if (rank >= 1) hx *= cosmo_a;
    if (rank >= 2) hy *= cosmo_a;
    if (rank >= 3) hz *= cosmo_a;
  }

  double inv_vol = 1.0 / (hx*hy*hz);


  //
  //
  // Set explosion properties
  //
  double energy_per_cell = E_51 * 1.0E51 / enzo_units->mass() /
                            (enzo_units->velocity() * enzo_units->velocity());
  double mass_per_cell   = m_eject * cello::mass_solar / enzo_units->mass();
  mass_per_cell /= ((double) number_of_feedback_cells_);
  energy_per_cell /= ((double) number_of_feedback_cells_);


  //  stencil_rad_ is integer separation from cell center
  //  and edge of injection region (i.e. 1 for 3x3 injection grid)
  if (  shift_cell_center_ &&
       ( ((xpos - (stencil_rad_+1)*hx) < xm) ||
         ((xpos + (stencil_rad_+1)*hx) > xp) ||
         ((ypos - (stencil_rad_+1)*hy) < ym) ||
         ((ypos + (stencil_rad_+1)*hy) > yp) ||
         ((zpos - (stencil_rad_+1)*hz) < zm) ||
         ((zpos + (stencil_rad_+1)*hz) > zp)    )  ) {

      xpos = std::min(  std::max(xpos, xm + (stencil_rad_ + 1 + 0.5)*hx),
                       xp - (stencil_rad_ + 1 + 0.5)*hx);
      ypos = std::min(  std::max(ypos, ym + (stencil_rad_ + 1 + 0.5)*hy),
                       yp - (stencil_rad_ + 1 + 0.5)*hy);
      zpos = std::min(  std::max(zpos, zm + (stencil_rad_ + 1 + 0.5)*hz),
                       zp - (stencil_rad_ + 1 + 0.5)*hz);
  }

  // compute coordinates of central feedback cell
  double xcell = (xpos - xm) / hx + gx - 0.5;
  double ycell = (ypos - ym) / hy + gy - 0.5;
  double zcell = (zpos - zm) / hz + gz - 0.5;

  // I believe -1's are needed (added from Fortran code due to index start differences )
  int ix       = ((int) floor(xcell + 0.5)) + gx - 1;
  int iy       = ((int) floor(ycell + 0.5)) + gy - 1;
  int iz       = ((int) floor(zcell + 0.5)) + gz - 1;

  double dxc   = ix + 0.5 - xcell;
  double dyc   = iy + 0.5 - ycell;
  double dzc   = iz + 0.5 - zcell;

  //
  // Set source velocity to local gas flow if all are left at
  // the default values
  //
  if (((pvx == pvy) && (pvy == pvz)) && (pvx == -9999.0)){
    int index = INDEX(ix,iy,iz,mx,my);
    pvx = v3[0][index]; pvy = v3[1][index]; pvz = v3[2][index];
  }

  double *u_local=0, *v_local=0, *w_local=0, *d_local=0, *ge_local=0, *te_local=0;
  double *ke_before=0, *metal_local=0;

  // number of feedback cells + 1 cell in each dimension
  int num_loc = (stencil_ + 1) * (stencil_ + 1) * (stencil_ + 1);
  if (!u_local){
    u_local     = new double[num_loc];
    v_local     = new double[num_loc];
    w_local     = new double[num_loc];
    d_local     = new double[num_loc];
    ge_local    = new double[num_loc];
    te_local    = new double[num_loc];
    metal_local = new double[num_loc];
    ke_before   = new double[num_loc];
  }

  // assign initial values to these
  for (int i = 0; i < num_loc; i++){
    u_local[i] = 0.0; v_local[i] = 0.0; w_local[i] = 0.0;
    d_local[i] = 0.0; te_local[i] = 0.0; ge_local[i] = 0.0;
    ke_before[i] = 0.0; metal_local[i] = 0.0;
  }

  double ke_f = 0.0;

  //
  if (ke_fraction < 0){

    // calculate variable kinetic energy fraction

    double avg_z = 0.0, avg_n = 0.0, avg_d = 0.0;

    for (int k = iz - stencil_rad_; k<= iz + stencil_rad_; k++){
      for (int j = iy - stencil_rad_; j <= iy + stencil_rad_; j++){
        for (int i = ix - stencil_rad_; i <= ix + stencil_rad_; i++){

//      AE TO DO: Actually calculate number density here with species

          int index = INDEX(i,j,k,mx,my);

          // NOTE: The presence of these statements throughout this routine
          //       is to generalize for situations where we don't have to
          //       kick particles away from grid edges once non-local
          //       blocks / processors know about particles that deposit
          //       feedback on their grids (this allows the loops to be
          //       simple - otherwise will have to continually recalc
          //       the min / max bounds of the loops to avoid edges )
          if ( (index < 0) || (index > nx*ny*nz)) continue;

          double mu_cell  = enzo_config->ppm_mol_weight;

          avg_z += metal[index]; // need to divide by d_tot below
          avg_n += d[index] / mu_cell;
          avg_d += d[index];
        }
      }
    } // end loop over local cells

    double inv_ncell = 1.0 / ((double) number_of_feedback_cells_);
    const double z_solar = 0.02;   // as assumed for these equations
    avg_z =  (avg_z / avg_d) / z_solar; // mass-weighted metallicity
    avg_n *= inv_ncell * enzo_units->density() / cello::mass_hydrogen; // in cgs
    avg_d *= inv_ncell * enzo_units->density(); // in cgs


    // Compute the time and radius of transtion to PDS phase for gas
    // with the computed properties.
    // t_PDS is in units of kyr  -    R_PDS is in units of pc

    double t_PDS = 0.0, R_PDS = 0.0;
    // For metal poor gas
    if (avg_z < 0.01){
      t_PDS = 306.0 * pow(E_51,0.125) * pow(avg_n,-0.75);
      R_PDS =  49.3 * pow(E_51,0.250) * pow(avg_n,-0.50);
    } else {
      t_PDS = 26.50 * pow(E_51,3.0/14.0) * pow(avg_z,-5.0/14.0) * pow(avg_n,-4.0/7.0);
      R_PDS = 18.50 * pow(E_51,2.0/7.0 ) * pow(avg_z,-1.0/7.0 ) * pow(avg_n,-3.0/7.0);
    } // end metallicity check

    double     R_resolve = hx*enzo_units->length() / cello::pc_cm;
    const double  n_resolve = 4.5; // number of cells needed to resolve R

    if (R_PDS > n_resolve * R_resolve){
      ke_f = 0.0;

    } else {

      ke_f = 3.97133E-6 * (avg_d / cello::mass_hydrogen) *
                   (1.0 / (R_resolve * R_resolve)) *
                   pow(R_PDS,7) * ( 1.0 / (t_PDS*t_PDS)) *
                   (1.0 / E_51);

    } // end resolved check

  } else {
    ke_f = kinetic_fraction_;

  }// end kinetic energy fraction check

  // apply kinetic energy fraction floor
  ke_f = ke_f < 1.0E-10 ? 0.0 : ke_f;

  double E_therm = (1.0 - ke_f) * energy_per_cell;

  // compute kinetic energy in the localized region on the grid before
  // the explosion
  int loc_index = 0;
  for (int k = iz - stencil_rad_ ; k <= iz + stencil_rad_ + 1; k++){
    for (int j = iy - stencil_rad_ ; j <= iy + stencil_rad_ + 1; j++){
      for (int i = ix - stencil_rad_ ; i <= ix + stencil_rad_ + 1; i++){

        int index = INDEX(i,j,k,mx,my);

        if ( (index < 0) || (index > nx*ny*nz)) continue;

        ke_before[loc_index] = 0.5 * d[index] * ( v3[0][index] * v3[0][index] +
                                                  v3[1][index] * v3[1][index] +
                                                  v3[2][index] * v3[2][index] );
        loc_index++;
      }
    }
  }

  // Now convert velocities in affected region to momentum in the
  // particle's reference frame
  this->convert_momentum(v3[0], v3[1], v3[2], d,
                         pvx, pvy, pvz,
                         nx, ny, nz,
                         ix, iy, iz, 1);

  // compute the total mass and energy in the cells before the explosion
  double sum_mass_init, sum_energy_init, sum_ke_init;
  this->sum_mass_energy(v3[0], v3[1], v3[2], d, ge, te,
                        nx, ny, nz, ix, iy, iz,
                        sum_mass_init, sum_energy_init, sum_ke_init);

  // compute the mass and momentum properties of adding feedback to the
  // psuedo-grid centered on the particle position before doing so for
  // the real grid. this is used to compute coefficient properties below
  // and prep for the CIC deposition
  this->add_feedback_to_grid(u_local, v_local, w_local, d_local,
                             ge_local, te_local, metal_local,
                             stencil_+1, stencil_+1, stencil_+1,          // nx,ny,nz for local grid
                             stencil_rad_, stencil_rad_, stencil_rad_, // local grid cell center  - should be 1 for 3x3x3 stencil (0,1,2) - 2 for 5x5x5 stencil (0,1, 2, 3,4)
                             dxc, dyc, dzc,
                             mass_per_cell, 1.0, 0.0);

  // momenum injection - compute coefficients for injection
  double mom_per_cell = 0.0;
  if (ke_f > 0){
    double A, B, C;
    this->compute_coefficients( v3[0], v3[1], v3[2], d, ge,
                                u_local, v_local, w_local, d_local,
                                nx, ny, nz, ix, iy, iz, A, B, C);

    A            = A - (sum_ke_init + ke_f*energy_per_cell);
    mom_per_cell = (-B + std::sqrt(B*B - 4.0 * A * C)) / (2.0 * C);
  } // add switch here? to do just momentum injection but no KE

  // Need to deposit metal(s) here

  this->add_feedback_to_grid(v3[0], v3[1], v3[2], d, ge, te, metal,
                             nx, ny, nz, ix, iy, iz, dxc, dyc, dzc,
                             mass_per_cell, mom_per_cell, E_therm);

  double sum_mass_final, sum_energy_final, sum_ke_final;

  this->sum_mass_energy(v3[0], v3[1], v3[2], d, ge, te,
                        nx, ny, nz, ix, iy, iz,
                        sum_mass_final, sum_energy_final, sum_ke_final);

  // AE NOTE:
  //    should do some error checking here
  //    with kinetic energy and momentum

  // Now convert momentum back to velocity
  this->convert_momentum(v3[0], v3[1], v3[2], d,
                         pvx, pvy, pvz,
                         nx, ny, nz, ix, iy, iz, 0);

  // Adjust total energy

  double ke_injected = 0.0;
  double delta_ke    = 0.0;
  double ke_after    = 0.0;
  loc_index = 0;
  for (int k = iz - stencil_rad_; k <= iz + stencil_rad_ + 1; k++){
    for (int j = iy - stencil_rad_; j <= iy + stencil_rad_ + 1; j++){
      for (int i = ix - stencil_rad_; i <= ix + stencil_rad_ + 1; i++){

        int index = (i + (j + k*ny)*nx);

        if (index < 0 || index >= nx*ny*nz) continue;

        ke_after = 0.5 * d[index] * ( v3[0][index] * v3[0][index] +
                                      v3[1][index] * v3[1][index] +
                                      v3[2][index] * v3[2][index]);

        delta_ke = ke_after - ke_before[loc_index];

        te[index] += delta_ke / d[index];

        ke_injected += delta_ke;

        loc_index++;
      }
    }
  }


  delete [] u_local;
  delete [] v_local;
  delete [] w_local;
  delete [] d_local;
  delete [] ge_local;
  delete [] te_local;
  delete [] metal_local;
  delete [] ke_before;

  return;
}


void EnzoMethodDistributedFeedback::convert_momentum(
                                    enzo_float * vx, enzo_float * vy, enzo_float * vz, enzo_float * d,
                                    const enzo_float & up, const enzo_float & vp, const enzo_float & wp,
                                    const int &mx, const int &my, const int &mz,
                                    const int &ix, const int & iy, const int & iz,
                                    int idir){

  int xo, yo, zo;
  xo = 1;
  yo = mx;
  zo = (mx * my);

  for (int k = -stencil_rad_; k <= stencil_rad_ + 1; k++){
    for (int j = -stencil_rad_; j <= stencil_rad_ + 1; j++){
      for (int i = -stencil_rad_; i <= stencil_rad_ + 1; i++){

        int index = (ix + i) + ( (iy + j) + (iz + k)*my)*mx;

        if ( (index < 0) || (index > mx*my*mz)) continue;

        if (idir >= 1){ // velocity -> momentum

          vx[index] = (vx[index] - up) * d[index];
          vy[index] = (vy[index] - vp) * d[index];
          vz[index] = (vz[index] - wp) * d[index];

        } else {

          vx[index] = vx[index] / d[index] + up;
          vy[index] = vy[index] / d[index] + vp;
          vz[index] = vz[index] / d[index] + wp;

        } // end dir check

      }
    }
  }


  return;
}

void EnzoMethodDistributedFeedback::sum_mass_energy(
                               enzo_float * px, enzo_float * py, enzo_float * pz, enzo_float * d,
                               enzo_float * ge, enzo_float * te,
                               const int & mx, const int & my, const int & mz,
                               const int & ix, const int & iy, const int & iz,
                               double & sum_mass, double & sum_energy, double & sum_ke){

  sum_mass = 0; sum_energy = 0; sum_ke = 0;

  for (int k = -stencil_rad_;  k <= stencil_rad_ + 1; k++){
    for (int j = -stencil_rad_; j <= stencil_rad_ + 1; j++){
      for (int i = -stencil_rad_; i <= stencil_rad_ + 1; i++){

        int index = (ix + i) + ( (iy + j) + (iz + k)*my)*mx;

        if ( (index < 0) || (index > mx*my*mz)) continue;

        double mass_term  = d[index];
        double mom_term   = px[index]*px[index] +
                            py[index]*py[index] +
                            pz[index]*pz[index];
        double ke         = mom_term / (2.0 * mass_term);

        sum_mass += mass_term;
        sum_ke   += ke;

        double energy = 0.0;
        if (dual_energy_){
          energy = ge[index]*d[index];
        } else {
          energy = te[index]*d[index] - ke;
        }

        sum_energy += ke + energy;

      }
    }
  } // end loop

  return;
}

void EnzoMethodDistributedFeedback::add_feedback_to_grid(
                                enzo_float * px, enzo_float * py, enzo_float * pz,
                                enzo_float * d, enzo_float * ge, enzo_float * te,
                                enzo_float * metal,
                                const int & mx, const int & my, const int &mz,
                                const int & ix, const int & iy, const int &iz,
                                const double & dxc, const double & dyc, const double & dzc,
                                const double & mass_per_cell, const double & mom_per_cell,
                                const double & therm_per_cell){





  for (int k = -stencil_rad_; k <= stencil_rad_; k++){
    for (int j = -stencil_rad_; j <= stencil_rad_; j++){
      for (int i = -stencil_rad_; i <= stencil_rad_; i++){

        for (int i1 = i; i1 <= i + 1; i1++){
          double dxc1 = dxc;
          if (i1 == (i + 1)) dxc1 = 1.0 - dxc;

          for (int j1 = j; j1 <= j + 1; j1++){
            double dyc1 = dyc;
            if (j1 == (j + 1)) dyc1 = 1.0 - dyc;

            for (int k1 = k; k1 <= k + 1; k1++){
              double dzc1 = dzc;
              if (k1 == (k + 1)) dzc1 = 1.0 - dzc;

              double delta_mass =    mass_per_cell * dxc1 * dyc1 * dzc1;
              // use sign of i,j,k to assign direction. No momentum in cener
              // cell
              double delta_pu = ( (i > 0) ? 1 : (i < 0 ? -1 : 0)) * mom_per_cell * dxc1 * dyc1 * dzc1;
              double delta_pv = ( (j > 0) ? 1 : (j < 0 ? -1 : 0)) * mom_per_cell * dxc1 * dyc1 * dzc1;
              double delta_pw = ( (k > 0) ? 1 : (k < 0 ? -1 : 0)) * mom_per_cell * dxc1 * dyc1 * dzc1;

              double delta_therm = therm_per_cell * dxc1 * dyc1 * dzc1;

              int index = ( ix + i1 ) + ( (iy + j1) + (iz + k1)*my)*mx;

              if ( (index < 0) || (index > mx*my*mz)) continue;

              double inv_dens = 1.0 / (d[index] + delta_mass);

              // scale factor to account for the fact that delta_p's may be
              // zero for cardinal directions along cardinal axes, but we
              // still want total momentum change (|delta_p|) to be the same
              // for all cells
              int mom_norm = 3;
              mom_norm -= (i == 0) ? 1 : 0;
              mom_norm -= (j == 0) ? 1 : 0;
              mom_norm -= (k == 0) ? 1 : 0;

              // but set to 1 if in center cell since delta's above are zero
              // already and we don't want to divide by zero below
              mom_norm = (mom_norm == 0) ? 1 : mom_norm;
              double mom_scale = 1.0 / ( sqrt((double) mom_norm));

              px[index] +=  delta_pu * mom_scale;
              py[index] +=  delta_pv * mom_scale;
              pz[index] +=  delta_pw * mom_scale;

              te[index]  = (te[index]*d[index] + delta_therm) * inv_dens;

              if (dual_energy_)
                  ge[index] = (ge[index]*d[index] + delta_therm)*inv_dens;

              d[index] = d[index] + delta_mass;

              if(metal) metal[index] = (metal[index] +
                                       delta_mass * ejecta_metal_fraction_);

              // account for multi-species (H,He,etc.) here, along with additional
              // metal species fields if they are present

            } // end k1 loop
          } // end j1 loop
        } // end i1 loop


      } // end i loop
    } // end j loop
  } // end k loop


  return;
}

void EnzoMethodDistributedFeedback::compute_coefficients(
                           enzo_float *px, enzo_float *py, enzo_float *pz, enzo_float *d,
                           enzo_float *ge, enzo_float *px_l, enzo_float *py_l, enzo_float *pz_l,
                           enzo_float *d_l,
                           const int & mx, const int &my, const int &mz,
                           const int &ix, const int &iy, const int &iz,
                           double &A, double &B, double &C){

  A = 0.0; B = 0.0; C = 0.0;

  double mass_term = 0.0, mom_term = 0.0, b_term = 0.0, c_term = 0.0;

  int loc_index = 0;

  for (int k = -stencil_rad_; k <= stencil_rad_+ 1; k++){
    for (int j = -stencil_rad_; j <= stencil_rad_ + 1; j++){
      for (int i = -stencil_rad_; i <= stencil_rad_ + 1; i++){

        int index = (ix + i) + ( (iy + j) + (iz + k)*my)*mx;

        if ( (index < 0) || (index >= mx*my*mz)){
          loc_index++;
          continue;
        }

        double mass_term = d[index];
        double mom_term  = px[index]*px[index] + py[index]*py[index] +
                           pz[index]*pz[index];

        mass_term += d_l[loc_index];

        b_term     = px[index]*px_l[loc_index] + py[index]*py_l[loc_index] +
                     pz[index]*pz_l[loc_index];

        c_term     = px_l[loc_index] * px_l[loc_index] +
                     py_l[loc_index] * py_l[loc_index] +
                     pz_l[loc_index] * pz_l[loc_index];

        double inv_mass = 1.0 / mass_term;
        A         += 0.5 * mom_term * inv_mass;
        B         += b_term * inv_mass;
        C         += 0.5 * c_term * inv_mass;

        loc_index++;
      }
    }
  } // end loop


  return;
}

double EnzoMethodDistributedFeedback::timestep (Block * block) const throw()
{
  // In general this is not needed, but could imagine putting timestep
  // limiters in situations where, for example, one would want
  // dt < star_lifetime (or something like that), especially if
  // important things happen throughout the star's lifetime.
  EnzoUnits * enzo_units = enzo::units();

//  return 1000.0 * cello::yr_s / enzo_units->time();
  return std::numeric_limits<double>::max();
}
