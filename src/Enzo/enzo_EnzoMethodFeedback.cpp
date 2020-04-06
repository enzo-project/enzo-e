/// See LICENSE_CELLO file for license and copyright information

///
///
///
///

#include "cello.hpp"
#include "enzo.hpp"

EnzoMethodFeedback::EnzoMethodFeedback
()
  : Method()
{
  FieldDescr * field_descr = cello::field_descr();
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  // AJE: This was the old way this was done
  // Initialize default refresh object
  // const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
  //                           enzo_sync_id_method_feedback);
  // refresh(ir)->add_all_fields();

  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  ejecta_mass_   = enzo_config->method_feedback_ejecta_mass * cello::mass_solar /
                      enzo_units->mass();

  ejecta_energy_ = enzo_config->method_feedback_supernova_energy * 1.0E51 /
                   enzo_units->mass() / enzo_units->velocity() /
                   enzo_units->velocity();

  ejecta_metal_fraction_ = enzo_config->method_feedback_ejecta_metal_fraction;

  return;
}

void EnzoMethodFeedback::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | ejecta_mass_;
  p | ejecta_energy_;
  p | ejecta_metal_fraction_;

  return;
}

void EnzoMethodFeedback::compute (Block * block) throw()
{

  if (block->is_leaf()){
    this->compute_(block);
  }

  block->compute_done();

  return;
}

void EnzoMethodFeedback::compute_ (Block * block) throw()
{

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  Field field = block->data()->field();
  // Obtain grid sizes and ghost sizes

  enzo_float * d           = (enzo_float *) field.values ("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");

  enzo_float * metal = field.is_field("metal_density") ?
                       (enzo_float *) field.values("metal_density") : NULL;

  int mx, my, mz, gx, gy, gz;
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

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
    const int dl = particle.stride(it, ia_l);
    const int dc = particle.stride(it, ia_c);

    const int nb = particle.num_batches(it);

    for (int ib=0; ib<nb; ib++){
      enzo_float *px=0, *py=0, *pz=0;
      enzo_float *plifetime=0, *pcreation=0, *pmass=0;

      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);

      px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
      py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
      pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

      plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
      pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib);

      int np = particle.num_particles(it,ib);

      for (int ip=0; ip<np; ip++){
        // AE: Check and see if these differ....
        int ipdp = ip*dp;
        int ipdm = ip*dm;
        int ipdl = ip*dl;
        int ipdc = ip*dc;

        // negative lifetime are particles that have already gone SN
        // creation time must be > 0
        // only go SN if age >= lifetime
        if ( (plifetime[ipdl] <= 0.0) || (pcreation[ipdc] <= 0.0) ||
             (current_time - pcreation[ipdc]) < plifetime[ipdl]) continue;

        plifetime[ipdl] *= -1.0;                  // set to negative - flag for already gone SNe
        pmass[ipdm] = pmass[ipdm] - ejecta_mass_; // subtract mass from particle

        // get corresponding grid position
        double xp = (px[ipdp] - xm) / hx;
        double yp = (py[ipdp] - ym) / hy;
        double zp = (pz[ipdp] - zm) / hz;

        // get 3D grid index for particle - account for ghost zones!!
        int ix = ((int) std::floor(xp))  + gx;
        int iy = ((int) std::floor(yp))  + gy;
        int iz = ((int) std::floor(zp))  + gz;

        // now deposit feedback in this cell
        int i = INDEX(ix,iy,iz,mx,my);

        // rescale tracer fields to maintain constant mass fraction
        // with the corresponding new density
        double d_old = d[i];
        d[i] += ejecta_mass_ * inv_vol;
        // catch here for multi-species
        if (metal) metal[i] += ejecta_metal_fraction_ *
                               ejecta_mass_ * inv_vol;

        //double scale = d[i] / d_old;
        //EnzoMethodStarMaker::rescale_densities(enzo_block, i, scale);

        // inject energy
        te[i] += ejecta_energy_ * inv_vol;
        if (enzo_config->ppm_dual_energy) ge[i] += ejecta_energy_ * inv_vol;

        count++;
      } // end loop over particles
    } // end loop over batches


  }

  if (count > 0){
      CkPrintf("Number of feedback particles:  %i \n",count);
  }


  return;
}

double EnzoMethodFeedback::timestep (Block * block) const throw()
{
  // In general this is not needed, but could imagine putting timestep
  // limiters in situations where, for example, one would want
  // dt < star_lifetime (or something like that), especially if
  // important things happen throughout the star's lifetime.

  return std::numeric_limits<double>::max();
}
