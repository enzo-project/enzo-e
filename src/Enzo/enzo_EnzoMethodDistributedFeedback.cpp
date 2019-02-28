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

  total_ejecta_mass_   = enzo_config->method_feedback_ejecta_mass * cello::mass_solar /
                         enzo_units->mass();

  total_ejecta_energy_ = enzo_config->method_feedback_supernova_energy * 1.0E51 /
                         enzo_units->mass() / enzo_units->velocity() /
                         enzo_units->velocity();

  ejecta_metal_fraction_ = enzo_config->method_feedback_ejecta_metal_fraction;

  stencil_size_             = enzo_config->method_feedback_stencil;
  istencil_                 = ( (int) (stencil_size_ - 1 / 2.0));
  number_of_feedback_cells_ = stencil_size_ * stencil_size_ * stencil_size_;


  m_eject_    = total_ejecta_mass_ / ((double) number_of_feedback_cells_);

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

  p | stencil_size_;
  p | istencil_;
  p | number_of_feedback_cells_;
  p | m_eject_;

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

void EnzoMethodDistributedFeedback::compute_ (Block * block) throw()
{

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  Field field = block->data()->field();
  // Obtain grid sizes and ghost sizes

  enzo_float * d           = (enzo_float *) field.values ("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");

  enzo_float * v3[3]       = { (enzo_float *) field.values("velocity_x"),
                               (enzo_float *) field.values("velocity_y"),
                               (enzo_float *) field.values("velocity_z")   };

  enzo_float * metal = field.is_field("metal_density") ?
                       (enzo_float *) field.values("metal_density") : NULL;

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

  Particle particle = enzo_block->data()->particle();

  // apply feedback depending on particle type
  // for now, just do this for all star particles

  int it = particle.type_index("star");
  int count = 0;

  if (particle.num_particles(it) > 0 ){

    const int ia_m = particle.attribute_index (it, "mass");

    const int ia_x = (rank >= 1) ? particle.attribute_index (it, "x") : -1;
    const int ia_y = (rank >= 2) ? particle.attribute_index (it, "y") : -1;
    const int ia_z = (rank >= 3) ? particle.attribute_index (it, "z") : -1;

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
        if ( (plifetime[ipdl] <= 0.0) || (pcreation[ipdc] <= 0.0) ||
             (current_time - pcreation[ipdc]) < plifetime[ipdl]) continue;

        // Update particle properties here if needed
        plifetime[ipdl] *= -1.0;                  // set to negative - flag for alreahy gone SNe
        pmass[ipdm] = pmass[ipdm] - ejecta_mass_; // subtract mass from particle


        // get corresponding grid position of particle
        // and shift it if it is too close to the grid boundaries
        double xpos = px[ipdp];
        double ypos = py[ipdp];
        double zpos = pz[ipdp];

        //  istencil is integer separation from cell center
        //  and edge of injection region (i.e. 1 for 3x3 injection grid)
        if (  enzo_config->method_feedback_shift_cell_center &&
             ( ((xfc - (istencil_+1)*hx) < xm) ||
               ((xfc + (istencil_+1)*hx) > xp) ||
               ((yfc - (istencil_+1)*hy) < ym) ||
               ((yfc + (istencil_+1)*hy) > yp) ||
               ((zfc - (istencil_+1)*hz) < zm) ||
               ((zfc + (istencil_+1)*hz) > zp)    )  ) {

            xpos = std::min(  std::max(xpos, xm + (istencil_ + 1 + 0.5)*hx),
                             xp - (istencil_ + 1 + 0.5)*hx));
            ypos = std::min(  std::max(ypos, ym + (istencil_ + 1 + 0.5)*hy),
                             yp - (istencil_ + 1 + 0.5)*hxy));
            zpos = std::min(  std::max(zpos, zm + (istencil_ + 1 + 0.5)*hz),
                             zp - (istencil_ + 1 + 0.5)*hz));

            return;
        }

        double xface = (xpos - xm) / hx + gx;
        double yface = (ypos - ym) / hy + gy;
        double zface = (zpos - zm) / hz + gz;

        int ixface    = ((int) floor(xface + 0.5)) + gx;
        int iyface    = ((int) floor(yface + 0.5)) + gy;
        int izface    = ((int) floor(zface + 0.5)) + gz;

        double dxf   = ixface + 0.5 - xface;
        double dyf   = iyface + 0.5 - yface;
        double dzf   = izface + 0.5 - zface;

        // compute coordinates of central feedback cell
        double xcell = (xpos - xm) / hx + gx;
        double ycell = (ypos - ym) / hy + gy;
        double zcell = (zpos - zm) / hz + gz;

        int ix       = ((int) floor(xpos + 0.5)) + gx;
        int iy       = ((int) floor(ypos + 0.5)) + gy;
        int iz       = ((int) floor(zpos + 0.5)) + gz;

        double dxc   = ix + 0.5 - xpos;
        double dyc   = iy + 0.5 - ypos;
        double dzc   = iz + 0.5 - zpos;

        // allocate local fields to keep track of energy
        // and momentum for conservation
        double *u_local, *v_local, *w_local, *d_local, *ge_local, *te_local;
        double *ke_before;

        // number of feedback cells + extra in case
        // feedback grid (which is centered on particle position) is
        // misaligned with the actual grid
        int num_loc = number_of_feedback_cells_ + 2*2*2;
        u_local   = new double[num_loc];
        v_local   = new double[num_loc];
        w_local   = new double[num_loc];
        d_local   = new double[num_loc];
        ge_local  = new double[num_loc];
        te_local  = new double[num_loc];
        ke_before = new double[num_loc];

        // assign initial values to these
        for (int i = 0; i < num_loc; i++){
          u_local[i] = 0.0; v_local[i] = 0.0; w_local[i] = 0.0;
          d_local[i] = 0.0; te_local[i] = 0.0; ge_local[i] = 0.0;
          ke_before[i] = 0.0;
        }

        //

        // LEFT OFF HERE
        //
        // DO COMPUTATION HERE TO COMPUTE E_thermal, E_kin, and p_mom
        //
        //

        // compute kinetic energy in the localized rergion
        int loc_index = 0;
        for (int k = iz - istencil_ - 1; k <= iz + istencil_ + 1; k++){
          for (int j = iy - istencil_ - 1; j <= iy + istencil_ + 1; j++){
            for (int i = ix - istencil_ - 1; i <= ix + istencil_ + 1; i++){

              int index = INDEX(i,j,k,mx,my);

              ke_before[loc_index] = 0.5 * d[index] * ( v3[0][index] * v3[0][index] +
                                                        v3[1][index] * v3[1][index] +
                                                        v3[2][index] * v3[2][index] );
              loc_index++
            }
          }
        }

        // Now convert velocities in affected region to momentum in the
        // particle's reference frame
        this->convert_momentum(v3[0], v3[1], v3[2], d,
                               pvx[ipdv], pvy[ipdv], pvz[ipdv],
                               nx, ny, nz,
                               ix, iy, iz, 1);

        // compute the total mass and energy in the cells before the explosion
        double sum_mass_init, sum_energy_init, sum_ke_init;

        this->sum_mass_energy(v3[0], v3[1], v3[2], d, ge, te,
                              nx, ny, nz, ix, iy, iz,
                              &sum_mass_init, &sum_energy_init, &sum_ke_init);

        this->add_feedback_to_grid(u_local, v_local, w_local, d_local,
                                   ge_local, te_local,
                                   istencil_+2, istencil_+2, istencil_+2,
                                   istencil_+1, istencil_+1, istencil_+1,
                                   dxc, dyc, dzc, m_eject, 1.0, 0.0);

        // momenum injection
        double mom_per_cell = 0.0;
        if (E_kin > 0){
          double A, B, C;
          this->compute_coefficients( v3[0], v3[1], v3[2], d, ge,
                                      u_local, v_local, w_local, d_local,
                                      nx, ny, nz, ix, iy, iz, A, B, C);

          A            = A - (sum_ke_init + E_kin);
          mom_per_cell = (-B + std::sqrt(B*B - 4.0 * A * C)) / (2.0 * C);

        } // add switch here? to do just momentum injection but no KE

        // Need to deposit metal(s) here

        this->add_feedback_to_grid(v3[0], v3[1], v3[2], d, ge, te,
                                   nx, ny, nz, ix, iy, iz, dxc, dyc, dzc,
                                   m_eject, mom_per_cell, E_thermal);

        double sum_mass_final, sum_energy_final, sum_ke_final;

        this->sum_mass_energy(v3[0], v3[1], v3[2], d, ge, te,
                              nx, ny, nz, ix, iy, iz,
                              &sum_mass_final, &sum_energy_final, &sum_ke_final);

        // Now convert momentum back to velocity
        this->convert_momentum(v3[0], v3[1], v3[2], d,
                               pvx[ipdv], pvy[ipdv], pvz[ipdv],
                               nx, ny, nz, ix, iy, iz, 0);

        // Adjust total energy

        double ke_injected = 0.0;
        double delta_ke    = 0.0;
        double ke_after    = 0.0;
        loc_index = 0;
        for (int k = iz - istencil_ - 1; k <= iz + istencil_ + 1; k++){
          for (int j = iy - istencil_ - 1; j <= iy + istencil_ + 1; j++){
            for (int i = ix - istencil_ - 1; i <= ix + istencil_ + 1; i++){

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

        count++;
      } // end loop over particles
    } // end loop over batches


  }

  delete [] u_local;
  delete [] v_local;
  delete [] w_local;
  delete [] d_local;
  delete [] ge_local;
  delete [] te_local;
  delete [] ke_before;


  if (count > 0){
      std::cout << "Number of feedback particles:   " << count << "\n";
  }


  return;
}

void EnzoMethodDistributedFeedback::convert_momentum(
                                    double * vx, double * vy, double * vz, double * d,
                                    const double & up, const double & vp, const double & wp,
                                    const int &nx, const int &ny, const int &nz,
                                    const int &ix, const int & iy, const int & iz,
                                    int idir){

  int xo, yo, zo;
  xo = 1;
  yo = nx;
  zo = (nx * ny);

  for (int k = -istencil_ - 1; k <= istencil_ + 1; k++){
    for (int j = -istencil_ - 1; j <= istencil_ + 1; j++){
      for (int i = -istencil_ -1; i <= istencil_ + 1; i++){

        int index = (ix + 1) + ( (iy + j) + (iz + k)*ny)*nx;
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
                               double * px, double * py, double * pz, double * d,
                               double * ge, double * te,
                               const int & nx, const int & ny, const int & nz,
                               const int & ix, const int & iy, const int & iz,
                               double & sum_mass, double & sum_energy, double & sum_ke){

  *sum_mass = 0; *sum_energy = 0; *sum_ke = 0;

  for (int k = -istencil_ -1; k <= istencil_ + 1; k++){
    for (int j = -istencil_ - 1; j <= istencil_ + 1; j++){
      for (int i = -istencil_ - 1; i < = istencil_ + 1; i++){

        int index = (ix + i) + ( (iy + j) + (iz + k)*ny)*nx;

        if ( (index < 0) || (index > nx*ny*nz)) continue;

        double mass_term  = d[index];
        double mom_term   = px[index]*px[index] +
                            py[index]*py[index] +
                            pz[index]*pz[index];
        double ke         = mom_term / (2.0 * mass_term);

        *sum_mass += mass_term;
        *sum_ke   += ke;

        double energy = 0.0;
        if (enzo_config->ppm_dual_energy){
          energy = ge[index]*d[index];
        } else {
          energy = te[index]*d[index] - ke;
        }

        *sum_energy += ke + energy;

      }
    }
  } // end loop

  return;
}

void EnzoMethodDistributedFeedback::add_feedback_to_grid(
                                double * px, double * py, double * pz,
                                double * d, double * ge, double * te,
                                const int & nx, const int & ny, const int &nz,
                                const int & ix, const int & iy, const int &iz,
                                const double & dxc, const double & dyc, const double & dzc,
                                const double & mass_per_cell, const double & mom_per_cell,
                                const double & therm_per_cell){

  for (int k = -istencil_ - 1; k <= istencil_ + 1; k++){
    for (int j = -istencil_ - 1; j <= istencil_ + 1; j++){
      for (int i = -istencil_ - 1; i <= istencil_ + 1; i++){

        for (int i_loc )
      }
    }
  }


  return;
}

void EnzoMethodDistributedFeedback::compute_coefficients(
                           double *px, double *py, double *pz, double *d,
                           double *ge, double *px_l, double *py_l, double *pz_l,
                           double *d_l,
                           const int & nx, const int &ny, const int &nz,
                           const int &ix, const int &iy, const int &iz,
                           double &A, double &B, double &C){

  A = 0.0; B = 0.0; C = 0.0;

  double mass_term = 0.0, mom_term = 0.0, b_term = 0.0, c_term = 0.0;

  int loc_index = 0;

  for (int k = -istencil_ - 1; k <= istencil_ + 1; k++){
    for (int j = -istencil_ - 1; j <= istencil_ + 1; j++){
      for (int i = -istencil_ - 1; i <= istencil_ + 1; i++){

        int index = (ix + i) + ( (iy + j) + (iz + k)*ny)*nx;

        if ( (index < 0) || index >= nx*ny*nz)) continue;

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

  return std::numeric_limits<double>::max();
}
