
/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodFeedbackSTARSS.cpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
//              Will Hicks (whicks@ucsd.edu)
/// @date
/// @brief  Implements the STARSS model for stellar feedback
///         as implemented into Enzo as the MechStars model
///         by Azton Wells. This is intended to be a direct copy-over of the
///         MechStars methods.
///         TODO: List differences / changes here
///


#include "cello.hpp"
#include "enzo.hpp"

#include <time.h>

//#ifdef NOTDEFINED // for now... since not done coding

#define DEBUG_FEEDBACK_STARSS

// =============================================================================
// splice these off to a different file (later)

int EnzoMethodFeedbackSTARSS::determineSN(double age, int* nSNII, int* nSNIA,
                double mass_Msun, double tunit, float dt){

    const EnzoConfig * enzo_config = enzo::config();
    std::mt19937 mt(std::time(nullptr));

    if (NEvents > 0){
        *nSNII = 1;
        NEvents -= 1;
        return 1;
    }
    /* else, calculate SN rate, probability and determine number of events */
    *nSNII = 0;
    *nSNIA = 0;
    double RII=0, RIA=0, PII=0, PIA=0, random = 0;
    if (enzo_config->method_feedback_single_sn && NEvents < 0)
    {
        // printf("Calculating rates\n");
        /* age-dependent rates */
        if (age < 3.401)
        {
            RII = 0.0;
            RIA = 0.0;
        }
        if (3.401 <= age && age< 10.37)
        {
                RII = 5.408e-4;
                RIA = 0.0;
        }
        if (10.37 <= age && age < 37.53)
        {
                RII = 2.516e-4;
                RIA = 0.0;
        }
        if (37.53 <= age)
        {
                RII = 0.0;
                RIA = 5.3e-8+1.6e-5*exp(-0.5*pow((age-50.0)/10.0, 2));
        }
	//        fprintf(stdout, "Rates: %f %f %f\n", age, RII, RIA);
        /* rates -> probabilities */
        if (RII > 0){
        // printf("Zcpl = %e", zCouple);
            PII = RII * mass_Msun / cello::Myr_s *tunit*dt;
            // printf("PII =%f\n %f %e %f\n", PII, RII, massmass_solar, age);
            double random = double(mt())/double(mt.max());
            if (PII > 1.0 && enzo_config->method_feedback_unrestricted_sn){
                int round = (int)PII;
                *nSNII = round;
                PII -= round;
            }
            if (PII > 1.0 && !enzo_config->method_feedback_unrestricted_sn){
                ERROR("MethodFeedbackSTARSS::determineSN()", "PII too large!");
            }
            int psn = *nSNII;
            if (random < PII){
                *nSNII = psn+1;
            }
        }
        #ifdef DEBUG_FEEDBACK_STARSS
          CkPrintf("MethodFeedbackSTARSS::determineSN() -- mass_Msun = %f; age = %f; RII = %f; RIA = %f\n",
                    mass_Msun, age, RII, RIA);
        #endif
        // printf("RANDOM = %f\n", random);
        // printf("N SNII=%d\n",*nSNII);

        if (RIA > 0){
            PIA = RIA*mass_Msun / cello::Myr_s *tunit*dt;
            float random = float(rand())/float(RAND_MAX);

            if (PIA > 1.0 && enzo_config->method_feedback_unrestricted_sn)
            {
                int round = int(PIA);
                *nSNIA = round;
                PIA -= round;
            }
            int psn = *nSNIA;

            if (random < PIA)
                *nSNIA = psn+1;
            // if (*nSNIA > 0)
            //     fprintf(stdout, "PIA = %f\n", PIA);
        }
        #ifdef DEBUG_FEEDBACK_STARSS
          CkPrintf("MethodFeedbackSTARSS::determineSN() -- PII = %f; PIA = %f\n", PII, PIA);
        #endif
    }
 
        return 1;
}

int EnzoMethodFeedbackSTARSS::determineWinds(double age, double * eWinds, double * mWinds, double * zWinds,
                      double mass_Msun, double metallicity_Zsun, double tunit, double dt) 
    {
    // age in Myr

    const EnzoConfig * enzo_config = enzo::config();
    bool oldEnough = (age < 0.0001)?(false):(true);
    double windE = 0,  windM = 0, windZ = 0.0;
    double wind_factor = 0.0;
    double e_factor = 0.0;

    if (mass_Msun > 11 && oldEnough){

        if (0.001 < age && age < 1.0){
            wind_factor =4.763 * std::min((0.01 + metallicity_Zsun), 1.0) ;
        }
        if (1 <= age && age < 3.5){
            wind_factor = 4.763*std::min(0.01+metallicity_Zsun, 1.0)* 
                pow(age, 1.45+0.08*std::min(log(metallicity_Zsun), 1.0));
        }
        if (3.5 <= age && age < 100){
            wind_factor = 29.4*pow(age/3.5, -3.25)+0.0042;
        
        }
        if (age < 100){
            float d = powl(1+age/2.5, 1.4);
            float a50 = powl(double(age)/10.0, 5.0);
            e_factor = 5.94e4 / d + a50 +4.83;
            
        }
        if (100 <= age){
            e_factor = 4.83;
            wind_factor = 0.42*pow(age/1000, -1.1)/(19.81/log(age));
        }
        windM = mass_Msun * wind_factor; //Msun/Gyr
        windM = windM*dt*tunit/(1e3 * cello::Myr_s); //Msun
        // if (debug) printf("First winds mass = %e\nFrom wf = %f, dt=%f Z = %e\n", windM, wind_factor, dtFixed, zZsun);
        //printf("eFactor = %f age = %f\n", e_factor, age);
        if (windM > mass_Msun){
            CkPrintf("Winds too large Mw = %e, Mp = %e age=%f, Z = %e\n",
                windM, mass_Msun, age, metallicity_Zsun);
            windM = 0.125*mass_Msun; // limit loss to huge if necessary.
        }
        windZ = std::max(cello::metallicity_solar, 0.016+0.0041*std::max(metallicity_Zsun, 1.65)+0.0118)*windM;
        windE = e_factor * 1e12 * windM;


    *mWinds = windM;
    *zWinds = windZ;
    *eWinds = windE;

}

    return 1;

}

void EnzoMethodFeedbackSTARSS::transformComovingWithStar(enzo_float * density, 
                                  enzo_float * velocity_x, enzo_float * velocity_y, enzo_float * velocity_z,
                                  const enzo_float up, const enzo_float vp, const enzo_float wp,
                                  const int mx, const int my, const int mz, int direction) const throw()
{
  int size = mx*my*mz;
  if (direction > 0)
  {
    // to comoving with star
    // NOTE: This transforms the velocity field into a momentum
    //       field for the sake of depositing momentum easily 
    for (int ind = 0; ind<size; ind++) {
      double mult = density[ind];
      velocity_x[ind] = (velocity_x[ind]-up)*mult;
      velocity_y[ind] = (velocity_y[ind]-vp)*mult;
      velocity_z[ind] = (velocity_z[ind]-wp)*mult;
    }
  }

  else if (direction < 0)
  {
    // back to "lab" frame
    for (int ind = 0; ind<size; ind++) {
      double mult = 1/density[ind];
      velocity_x[ind] = velocity_x[ind]*mult + up;
      velocity_y[ind] = velocity_y[ind]*mult + vp;
      velocity_z[ind] = velocity_z[ind]*mult + wp;
    }
  }

}

void EnzoMethodFeedbackSTARSS::createCouplingParticles(EnzoBlock * enzo_block, const int nCouple,
                                 double coupledEnergy, double coupledGasEnergy,
                                 double coupledMass, double coupledMetals, double coupledMomenta,
                                 enzo_float * xcp, enzo_float * ycp, enzo_float * zcp,
                                 enzo_float * unitx, enzo_float * unity, enzo_float * unitz) const throw()
{
  Particle particle = enzo_block->data()->particle();

  const int it = particle.type_index("starss_coupling");
  int ib = 0; // batch counter
  int ipp = 0; // counter

  //TODO: rename mass attribute to density???
  const int ia_m  = particle.attribute_index(it,"mass");
  const int ia_x  = particle.attribute_index(it,"x");
  const int ia_y  = particle.attribute_index(it,"y");
  const int ia_z  = particle.attribute_index(it,"z");
  const int ia_momentum_x  = particle.attribute_index(it,"px");
  const int ia_momentum_y  = particle.attribute_index(it,"py");
  const int ia_momentum_z  = particle.attribute_index(it,"pz"); 
  const int ia_e  = particle.attribute_index(it,"energy");
  const int ia_ge = particle.attribute_index(it,"gas_energy");
  const int ia_mz = particle.attribute_index(it,"metal_mass");

  enzo_float * pmass = 0;
  enzo_float * px = 0;
  enzo_float * py = 0;
  enzo_float * pz = 0;
  enzo_float * p_px = 0;
  enzo_float * p_py = 0;
  enzo_float * p_pz = 0;
  enzo_float * pe = 0;
  enzo_float * pge = 0;
  enzo_float * pmz = 0;

  double mult = 1.0/nCouple;

  for (int ip = 0; ip < nCouple; ip++) {
    int my_particle = particle.insert_particles(it,1);
    particle.index(my_particle, &ib, &ipp);
    int io = ipp;

    pmass = (enzo_float *) particle.attribute_array(it, ia_m , ib);
    px    = (enzo_float *) particle.attribute_array(it, ia_x , ib);
    py    = (enzo_float *) particle.attribute_array(it, ia_y , ib);
    pz    = (enzo_float *) particle.attribute_array(it, ia_z , ib);
    p_px    = (enzo_float *) particle.attribute_array(it, ia_momentum_x , ib);
    p_py    = (enzo_float *) particle.attribute_array(it, ia_momentum_y , ib);
    p_pz    = (enzo_float *) particle.attribute_array(it, ia_momentum_z , ib);
    pe    = (enzo_float *) particle.attribute_array(it, ia_e , ib);
    pge   = (enzo_float *) particle.attribute_array(it, ia_ge, ib);
    pmz   = (enzo_float *) particle.attribute_array(it, ia_mz, ib);

    pmass[io] = coupledMass * mult;
    px[io] = xcp[ip];
    py[io] = ycp[ip];
    pz[io] = zcp[ip];
    p_px[io] = coupledMomenta * unitx[ip];
    p_py[io] = coupledMomenta * unity[ip];
    p_pz[io] = coupledMomenta * unitz[ip];
    pe[io] = coupledEnergy * mult;
    pge[io] = coupledGasEnergy * mult;
    pmz[io] = coupledMetals * mult;
  }
}

 void EnzoMethodFeedbackSTARSS::deleteCouplingParticles(EnzoBlock * enzo_block) const throw()
{
  Particle particle = enzo_block->data()->particle();

  const int it = particle.type_index("starss_coupling");

  const int nb = particle.num_batches(it);
  int count = enzo_block->delete_particle_copies_(it);
  cello::simulation()->data_delete_particles(count);
} 

extern "C" void FORTRAN_NAME(cic_deposit)
  ( enzo_float (*px)[26], enzo_float (*py)[26], enzo_float (*pz)[26], const int * rank,
    const int * nCouple, enzo_float (*mass)[26], enzo_float * field, enzo_float (*left_edge)[3],
    int * mx, int * my, int * mz, double * hx, const double * cloudsize );
// =============================================================================
// =============================================================================
// =============================================================================

EnzoMethodFeedbackSTARSS::EnzoMethodFeedbackSTARSS
()
  : Method()
  , ir_feedback_(-1)
{
  FieldDescr * field_descr = cello::field_descr();
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  // required fields
  cello::define_field("density");
  cello::define_field("pressure");
  cello::define_field("internal_energy");
  cello::define_field("total_energy");
  cello::define_field("velocity_x");
  cello::define_field("velocity_y");
  cello::define_field("velocity_z");
  cello::define_field("metal_density");

  
  cello::define_field_in_group("metal_density","color");

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  sf_minimum_level_ = enzo_config->method_feedback_min_level;
  single_sn_        = enzo_config->method_feedback_single_sn;

  // Initialize refresh object
  cello::define_field("density_deposit");
  cello::define_field("velocity_x_deposit");
  cello::define_field("velocity_y_deposit");
  cello::define_field("velocity_z_deposit");
  cello::define_field("total_energy_deposit");
  cello::define_field("internal_energy_deposit");
  cello::define_field("metal_density_deposit"); 

  ir_feedback_ = add_refresh_();
  cello::simulation()->refresh_set_name(ir_feedback_,name()+":add");
  Refresh * refresh_fb = cello::refresh(ir_feedback_); 

  refresh_fb->set_accumulate(true);


  refresh_fb->add_field_src_dst
    ("density_deposit","density_deposit_copy");
  refresh_fb->add_field_src_dst
    ("velocity_x_deposit", "velocity_x_deposit_copy");
  refresh_fb->add_field_src_dst
    ("velocity_y_deposit", "velocity_y_deposit_copy");
  refresh_fb->add_field_src_dst
    ("velocity_z_deposit","velocity_z_deposit_copy");
  refresh_fb->add_field_src_dst
    ("total_energy_deposit","total_energy_deposit_copy");
  refresh_fb->add_field_src_dst
    ("internal_energy_deposit","internal_energy_deposit_copu");
  refresh_fb->add_field_src_dst
    ("metal_density_deposit","metal_density_deposit_copy");
 
  refresh_fb->set_callback(CkIndex_EnzoBlock::p_method_feedback_starss_end());
 
  std::mt19937 mt(std::time(nullptr)); // initialize random function
  
  // initialize NEvents parameter (mainly for testing). Sets off 'NEvents' supernovae,
  // with at most one supernova per star particle per cycle.
  this->NEvents = enzo_config->method_feedback_NEvents; 
  return;
}

void EnzoMethodFeedbackSTARSS::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | sf_minimum_level_;
  p | single_sn_;
  p | NEvents;
  p | ir_feedback_;

  return;
}

void EnzoMethodFeedbackSTARSS::compute (Block * block) throw()
{

  if (block->is_leaf()){
    this->compute_(block);
  }

  block->compute_done();

  return;
}

void EnzoMethodFeedbackSTARSS::add_accumulate_fields(EnzoBlock * enzo_block) throw()
{
  int mx, my, mz, gx, gy, gz, nx, ny, nz;

  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  Field field = enzo_block->data()->field();
  field.size(&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);
  enzo_block->data()->lower(&xm,&ym,&zm);
  enzo_block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;

  // add accumulated values over and reset them to zero

  enzo_float * d  = (enzo_float *) field.values("density");
  enzo_float * te = (enzo_float *) field.values("total_energy");
  enzo_float * ge = (enzo_float *) field.values("internal_energy");
  enzo_float * mf = (enzo_float *) field.values("metal_density");
  enzo_float * vx = (enzo_float *) field.values("velocity_x");
  enzo_float * vy = (enzo_float *) field.values("velocity_y");
  enzo_float * vz = (enzo_float *) field.values("velocity_z");
  
  enzo_float * d_dep  = (enzo_float *) field.values("density_deposit");
  enzo_float * te_dep = (enzo_float *) field.values("total_energy_deposit");
  enzo_float * ge_dep = (enzo_float *) field.values("internal_energy_deposit");
  enzo_float * mf_dep = (enzo_float *) field.values("metal_density_deposit");
  enzo_float * vx_dep = (enzo_float *) field.values("velocity_x_deposit");
  enzo_float * vy_dep = (enzo_float *) field.values("velocity_y_deposit");
  enzo_float * vz_dep = (enzo_float *) field.values("velocity_z_deposit");

  enzo_float * d_dep_c  = (enzo_float *) field.values("density_deposit_copy");
  enzo_float * te_dep_c = (enzo_float *) field.values("total_energy_deposit_copy");
  enzo_float * ge_dep_c = (enzo_float *) field.values("internal_energy_deposit_copy");
  enzo_float * mf_dep_c = (enzo_float *) field.values("metal_density_deposit_copy");
  enzo_float * vx_dep_c = (enzo_float *) field.values("velocity_x_deposit_copy");
  enzo_float * vy_dep_c = (enzo_float *) field.values("velocity_y_deposit_copy");
  enzo_float * vz_dep_c = (enzo_float *) field.values("velocity_z_deposit_copy");

//  for (int iz=gz; iz<nz+gz; iz++){
//    for (int iy=gy; iy<ny+gy; iy++){
//      for (int ix=gx; ix<nx+gx; ix++){
//        int i = INDEX(ix,iy,iz,mx,my);
      for (int i=0; i<mx*my*mz; i++){
       
        if (te_dep_c[i] <= 0) continue; 

          d [i] +=  d_dep_c[i];
          mf[i] += mf_dep_c[i];

          double cell_mass = d[i]*hx*hy*hz;
          double M_scale = d_dep_c[i]/d[i];
          te[i] += te_dep_c[i] / cell_mass; 
          ge[i] += ge_dep_c[i] / cell_mass;
          vx[i] += vx_dep_c[i] * M_scale;
          vy[i] += vy_dep_c[i] * M_scale;
          vz[i] += vz_dep_c[i] * M_scale;

      
           d_dep[i] = 0;
          mf_dep[i] = 0;
          te_dep[i] = 0;
          ge_dep[i] = 0;
          vx_dep[i] = 0;
          vy_dep[i] = 0;
          vz_dep[i] = 0;

           d_dep_c[i] = 0;
          mf_dep_c[i] = 0;
          te_dep_c[i] = 0;
          ge_dep_c[i] = 0;
          vx_dep_c[i] = 0;
          vy_dep_c[i] = 0;
          vz_dep_c[i] = 0;

         }
        
//      }
//    }
//  }
/*
  const EnzoConfig * enzo_config = enzo::config();

  // recompute the temperature
  EnzoComputeTemperature compute_temperature
    (enzo_config->ppm_density_floor,
     enzo_config->ppm_temperature_floor,
     enzo_config->ppm_mol_weight,
     enzo_config->physics_cosmology);

  compute_temperature.compute(enzo_block);
*/
  return;
}
void EnzoBlock::p_method_feedback_starss_end() 
{  
  EnzoMethodFeedbackSTARSS * method = static_cast<EnzoMethodFeedbackSTARSS*> (this->method());
  method->add_accumulate_fields(this);

  return;
}

void EnzoMethodFeedbackSTARSS::compute_ (Block * block)
{

  //----------------------------------------------------
  // some constants here that might get moved to parameters or something else
  const float SNII_ejecta_mass = 10.5;
  const float SNIa_ejecta_mass = 1.4;
  const float z_solar          = cello::metallicity_solar; //0.02 // Solar metal fraction (fire2)
  //-----------------------------------------------------

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  Particle particle = enzo_block->data()->particle();
  EnzoUnits * enzo_units = enzo::units();

  double current_time  = block->time();

  Field field = enzo_block->data()->field();

  enzo_float * d           = (enzo_float *) field.values("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");
  enzo_float * mf          = (enzo_float *) field.values("metal_density");

  // Obtain grid sizes and ghost sizes
  int mx, my, mz, gx, gy, gz, nx, ny, nz;
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  field.size(&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;

//  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
//  enzo_float cosmo_a = 1.0;

  const int rank = cello::rank();
/*
  if (cosmology) {
    enzo_float cosmo_dadt = 0.0;
    double dt    = block->dt();
    cosmology->compute_expansion_factor(&cosmo_a,&cosmo_dadt,current_time+0.5*dt);
    if (rank >= 1) hx *= cosmo_a;
    if (rank >= 2) hy *= cosmo_a;
    if (rank >= 3) hz *= cosmo_a;
  }
*/
  // apply feedback depending on particle type
  // for now, just do this for all star particles

  int it = particle.type_index("star");

  // AE: Copying this from MechStars_FeedbackRoutine

  int numSN = 0; // counter of SN events
  int count = 0; // counter of particles
 
/*  if (particle.num_particles(it) <= 0){
    // refresh
    cello::refresh(ir_feedback_)->set_active(enzo_block->is_leaf());
    enzo_block->refresh_start(ir_feedback_, CkIndex_EnzoBlock::p_method_feedback_starss_end());
    return;
  }
*/
  const int ia_m = particle.attribute_index (it, "mass");
  const int ia_x  = particle.attribute_index (it, "x");
  const int ia_y  = particle.attribute_index (it, "y");
  const int ia_z  = particle.attribute_index (it, "z");
  const int ia_vx = particle.attribute_index (it, "vx");
  const int ia_vy = particle.attribute_index (it, "vy");
  const int ia_vz = particle.attribute_index (it, "vz");
  const int ia_l = particle.attribute_index (it, "lifetime");
  const int ia_c = particle.attribute_index (it, "creation_time");
  const int ia_mf = particle.attribute_index (it, "metal_fraction");
  const int ia_sn = particle.attribute_index (it, "number_of_sn"); // name change?

  const int dm = particle.stride(it, ia_m);
  const int dp = particle.stride(it, ia_x);
  const int dv = particle.stride(it, ia_vx);
  const int dl = particle.stride(it, ia_l);
  const int dc = particle.stride(it, ia_c);
  const int dmf = particle.stride(it, ia_mf);
  const int dsn = particle.stride(it, ia_sn);

  const int nb = particle.num_batches(it);

  for (int ib=0; ib<nb; ib++){
    enzo_float *px=0, *py=0, *pz=0, *pvx=0, *pvy=0, *pvz=0;
    enzo_float *plifetime=0, *pcreation=0, *pmass=0, *pmetal=0, *psncounter=0;

    pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
    pmetal = (enzo_float *) particle.attribute_array(it, ia_mf, ib);

    px  = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py  = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz  = (enzo_float *) particle.attribute_array(it, ia_z, ib);
    pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
    pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
    pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

    plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
    pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib);

    psncounter = (enzo_float *) particle.attribute_array(it, ia_sn, ib);

    int np = particle.num_particles(it,ib);

    for (int ip=0; ip<np; ip++){
      // AE: Check and see if these actually differ.... (ask James?)
      int ipdp = ip*dp; // pos
      int ipdm = ip*dm; // mass
      int ipdv = ip*dv; // velocity
      int ipdl = ip*dl; // lifetime
      int ipdc = ip*dc; // creation time
      int ipdmf = ip*dmf; // metallicity
      int ipsn  = ip*dsn; // number of SNe counter

      if (pmass[ipdm] > 0.0 && plifetime[ipdl] > 0.0){
        const double age = (current_time - pcreation[ipdc]) * enzo_units->time() / cello::Myr_s;
        count++; // increment particles examined here

        // compute coordinates of central feedback cell
        // this must account for ghost zones
        double xcell = (px[ipdp] - xm) / hx + gx - 0.5;
        double ycell = (py[ipdp] - ym) / hy + gy - 0.5;
        double zcell = (pz[ipdp] - zm) / hz + gz - 0.5;

        int ix       = ((int) floor(xcell + 0.5));
        int iy       = ((int) floor(ycell + 0.5));
        int iz       = ((int) floor(zcell + 0.5));

        int index = INDEX(ix,iy,iz,mx,my); // 1D cell index of star position

        //if(pmass[ipdm] * enzo_units->mass_units() / cello::mass_solar
        //       < enzo_config->method_star_maker_minimum_star_mass){
        //
        //}
        // skipping continual formation routines here (May need later)

        int nSNII = 0, nSNIa = 0;
        double SNMassEjected = 0.0, SNMetalEjected = 0.0;

        /* determine how many supernova events */
        if (single_sn_){

          /* Determine SN events from rates (currently taken from Hopkins 2018) */

          determineSN(age, &nSNII, &nSNIa, pmass[ipdm] * enzo_units->mass() / cello::mass_solar,
                      enzo_units->time(), block->dt());

          numSN += nSNII + nSNIa;
#ifdef DEBUG_FEEDBACK_STARSS
          if (nSNII + nSNIa > 0){
            CkPrintf("MethodFeedbackSTARSS SNe %d %d level = %d age = %f\n", nSNII, nSNIa, block->level(), age);
          }
#endif
          /* AJE: Can I just change this to a single call to deposit feedback for
                  the total ejecta of SNe and winds ? */

          if (nSNII+nSNIa > 0){
            /* set feedback properties based on number and types of SN */
            double energySN = (nSNII+nSNIa) * 1.0e51;

            SNMassEjected = SNII_ejecta_mass * nSNII +
                            SNIa_ejecta_mass * nSNIa; // AE: split this in two channels for yields

            const double starZ = pmetal[ipdmf] / z_solar;

            // TODO:compute metal mass ejected here

            //SNMetalEjected = ;

            /* Fixed mass ejecta */


            this->deposit_feedback( block, energySN, SNMassEjected, SNMetalEjected,
                                    pvx[ipdv],pvy[ipdv],pvz[ipdv],
                                    px[ipdp],py[ipdp],pz[ipdp],
                                    ix, iy, iz, 0, nSNII, nSNIa, starZ); // removed P3

            // add counter for number of SNe for the particle
            psncounter[ipsn] += nSNII+nSNIa;



          } // if nSNII or nSNIa > 0
        } // if single SN

        // Now do stellar winds
        double windMass=0.0, windMetals=0.0, windEnergy=0.0;
        if (enzo_config->method_feedback_stellar_winds){


          const double starZ = pmetal[ipdmf] / z_solar;


          determineWinds(age, &windEnergy, &windMass, &windMetals,
                         pmass[ipdm] * pmass[ipdm] * enzo_units->mass() / cello::mass_solar,
                         starZ, enzo_units->time(), block->dt());

          /* Error messages go here */

          if (windMass > 0){
          #ifdef DEBUG_FEEDBACK_STARSS
            CkPrintf("STARSS_FB: Adding stellar winds...\n");
          #endif
            this->deposit_feedback( block, windEnergy, windMass, windMetals,
                                    pvx[ipdv],pvy[ipdv],pvz[ipdv],
                                    px[ipdp],py[ipdp],pz[ipdp],
                                    ix, iy, iz, 1, 0, 0, 0.0); // removed P3

          } // if wind mass > 0
        } // if winds

        pmass[ipdm] -= std::max(0.0,
                       (windMass + SNMassEjected) /
                       (enzo_units->mass()/cello::mass_solar));


        //


      } // if mass and lifetime > 0

    } // end particle loop
  } // end batch loop


  if (count > 0){
    CkPrintf("FeedbackSTARSS: Num FB particles = %d  Events = %d  FeedbackTime %e\n",
              count, numSN, 0.00);
  }

  // refresh
  cello::refresh(ir_feedback_)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_start(ir_feedback_, CkIndex_EnzoBlock::p_method_feedback_starss_end());
  return;
}




// ----------------------------------------------------------------------------


void EnzoMethodFeedbackSTARSS::deposit_feedback (Block * block,
                                              double ejectaEnergy,
                                              double ejectaMass,
                                              double ejectaMetals,
                                              const enzo_float up, const enzo_float vp, const enzo_float wp,
                                              const enzo_float xp, const enzo_float yp, const enzo_float zp,
                                              const int ix, const int iy, const int iz,
                                              const int winds, const int nSNII,
                                              const int nSNIA,
                                              enzo_float starZ)

 const throw(){
  /*
   This routine will create a cube of coupling particles, where we determine
      the feedback quantities.  The vertices of the cube are ~coupled particles
      and all have radius dx from the source particle.
      Each vertex particle will then be CIC deposited to the grid!
  */
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();
  EnzoUnits * enzo_units = enzo::units();
  double vunit = enzo_units->velocity();
  double rhounit = enzo_units->density();
  double lunit = enzo_units->length();
  double munit = enzo_units->mass();
  double tunit = enzo_units->time();
  double eunit = munit*vunit*vunit; // energy
  double Tunit = enzo_units->temperature();

  const EnzoConfig * enzo_config = enzo::config();
  bool AnalyticSNRShellMass = enzo_config->method_feedback_analytic_SNR_shell_mass;

  double tiny_number = 1e-20;

  // Obtain grid sizes and ghost sizes
  int mx, my, mz, gx, gy, gz, nx, ny, nz;
  double xm, ym, zm, xp_, yp_, zp_, hx, hy, hz;
  field.size(&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);
  enzo_block->data()->lower(&xm,&ym,&zm);
  enzo_block->data()->upper(&xp_,&yp_,&zp_);
  field.cell_width(xm,xp_,&hx,ym,yp_,&hy,zm,zp_,&hz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;

  double size = mx*my*mz;

  //EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  //enzo_float cosmo_a = 1.0;

  const int rank = cello::rank();
/*
  if (cosmology) {
    enzo_float cosmo_dadt = 0.0;
    double dt    = block->dt();
    double current_time = block->time();
    cosmology->compute_expansion_factor(&cosmo_a,&cosmo_dadt,current_time+0.5*dt);
    if (rank >= 1) hx *= cosmo_a;
    if (rank >= 2) hy *= cosmo_a;
    if (rank >= 3) hz *= cosmo_a;
  }
*/
  double cell_volume = hx*hy*hz * enzo_units->volume();

  enzo_float * d           = (enzo_float *) field.values("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");

  enzo_float * vx          = (enzo_float *) field.values("velocity_x");
  enzo_float * vy          = (enzo_float *) field.values("velocity_y");
  enzo_float * vz          = (enzo_float *) field.values("velocity_z");

  enzo_float * mf          = (enzo_float *) field.values("metal_density");

  const int index = INDEX(ix,iy,iz,mx,my);

  int stretch_factor = 1.0;  //1.5/sin(M_PI/10.0);  // How far should cloud particles be from their host
                               // in units of dx. Since the cloud forms a sphere shell, stretchFactor > 1 is not recommended

  const int nCouple = 26; // 3x3x3 cube minus central cell
  const double A = stretch_factor * hx;

  /*
      Make a cloud of coupling particles;
      each is on a grid with dx spacing
      from the host particles.
  */

  enzo_float CloudParticleVectorX[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  enzo_float CloudParticleVectorY[] = {1, 1, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, 0, -1, -1, -1, 1, 1, 1, 0, 0, -1, -1, -1, 0};
  enzo_float CloudParticleVectorZ[] = {1, 0, -1, 1, -1, 1, 0, -1, 0, 1, 0, -1, 1, -1, 1, 0, -1, 1, 0, -1, 1, -1, 1, 0, -1, 0};
  enzo_float weightsVector[nCouple];
  /* Set position of feedback cloud particles */

  enzo_float CloudParticlePositionX[nCouple];
  enzo_float CloudParticlePositionY[nCouple];
  enzo_float CloudParticlePositionZ[nCouple];

  /*all possible values of x,y,z with origin particle at x=y=z=0.0 */
  for (int cpInd = 0; cpInd < nCouple; cpInd++)
  {
      double norm = sqrt(CloudParticleVectorX[cpInd] * CloudParticleVectorX[cpInd] +
                        CloudParticleVectorY[cpInd] * CloudParticleVectorY[cpInd] +
                        CloudParticleVectorZ[cpInd] * CloudParticleVectorZ[cpInd]);
      double inv_norm = 1.0 / norm;
      double xbaMag = A * A * norm * norm;
      /* in this cloud, take the coupling particle position as 0.5, 0.5, 0.5 */
      // get position, adding the width of ghost zones. Adding ghost zone width
      // because cic_deposit.F  
      CloudParticlePositionX[cpInd] = xp - CloudParticleVectorX[cpInd]*inv_norm * A;
      CloudParticlePositionY[cpInd] = yp - CloudParticleVectorY[cpInd]*inv_norm * A;
      CloudParticlePositionZ[cpInd] = zp - CloudParticleVectorZ[cpInd]*inv_norm * A;
      weightsVector[cpInd] = 1.0; //0.5 * (1. - 1. / (1. + 1. / 4. / cello::pi / 26. / xbaMag / cello::pi));
      /* turn the vectors back into unit-vectors */
      CloudParticleVectorZ[cpInd] *= inv_norm;
      CloudParticleVectorY[cpInd] *= inv_norm;
      CloudParticleVectorX[cpInd] *= inv_norm;
  }
  float weightsSum = 0.0;
  for (int wind = 0; wind < nCouple; wind++)
  {
      weightsSum += weightsVector[wind];
  }
  for (int wind = 0; wind < nCouple; wind++)
  {
      weightsVector[wind] /= weightsSum;
      if (weightsVector[wind] == 0 || isnan(weightsVector[wind]))
      {
          ERROR("EnzoMethodFeedbackSTARSS::deposit_feedback()","NaN weight Vector!")
      }
  }

/* AJE LEFT OFF HERE */

//TODO: Add refresh here after actual deposition

     /* 
         transform to comoving with the star and transform velocities to momenta
         for easy momentum deposition, since the fluid variable that PPM pushes
         around is velocity
     */
  
  this->transformComovingWithStar(d,vx,vy,vz,up,vp,wp,mx,my,mz, 1);

  /* 
     Use averaged quantities across multiple cells so that deposition is stable.
     vmean is used to determine whether the supernova shell calculation should proceed:
     M_shell > 0 iff v_shell > v_gas 
  */ 
  double Z_mean=0, d_mean=0, n_mean=0, v_mean=0, mu_mean=0;
  for (int ix_ = ix-1; ix_ < ix+2; ix_++) {
    for (int iy_ = iy-1; iy_ < iy+2; iy_++) {
      for (int iz_ = iz-1; iz_ < iz+2; iz_++) {
        int ind = INDEX(ix_,iy_,iz_, mx,my);
        Z_mean += mf[ind] / d[ind];
        mu_mean += enzo_config->ppm_mol_weight; // TODO: make EnzoComputeMolecularWeight? 
        d_mean += d[ind];
      } 
    }
  }
  
  v_mean *= vunit / 27; // cm/s!
  Z_mean *= 1/(27 * cello::metallicity_solar); // Zsun
  mu_mean /= 27;
  d_mean *= rhounit/27; // g/cm^3 
  n_mean = d_mean / (cello::mass_hydrogen/mu_mean);

#ifdef DEBUG_FEEDBACK_STARSS
  CkPrintf("STARSS_FB: Zmean = %e Dmean = %e (%e) mu_mean = %e ", Z_mean, d_mean, d_mean / rhounit, mu_mean);
  CkPrintf("Nmean = %f vmean = %f\n", n_mean, v_mean/1e5);
#endif

  double Z_Zsun = Z_mean;
  double fz = std::min(2.0, pow(std::max(Z_Zsun,0.01),-0.14));

  // Cooling radius as in Hopkins, but as an average over cells
  
  double CoolingRadius = 28.4 * pow(std::max(0.1, n_mean), -3.0 / 7.0) * pow(ejectaEnergy / 1.0e51, 2.0 / 7.0) * fz;
  double coupledEnergy = ejectaEnergy;
  double cellwidth = hx*lunit / cello::pc_cm; // pc
  double dxRatio = cellwidth / CoolingRadius;

  /* 
     We want to couple one of four phases: free expansion, Sedov-taylor, shell formation, or terminal
     The first three phases are take forms from Kim & Ostriker 2015, the last from Cioffi 1988
  */

 
  // if we resolve free expansion, all energy is thermally coupled
  double p_free = sqrt(2*ejectaMass*cello::mass_solar*ejectaEnergy) / cello::mass_solar/1e5;
  double r_free = 2.75*pow(ejectaMass / 3 / n_mean, 1.0/3); // free expansion radius eq. 2
  //bool use_free false; // could just deposit free expansion into host cell, really...
  
  double t3_sedov = pow(std::max(r_free, std::max(cellwidth, r_free)) * cello::pc_cm / 
                    (5.0 * cello::pc_cm * pow(ejectaEnergy / 1e51 / n_mean, 1.0 / 5.0)), 5./ 2.);
  
  double p_sedov = 2.21e4 * pow(ejectaEnergy / 1e51, 4. / 5.) * 
                   pow(n_mean, 1. / 5.) * pow(t3_sedov, 3. / 5.); // eq 16

  // shell formation radius eq. 8
  double r_shellform = 22.6*pow(ejectaEnergy/1e51,0.29) * pow(n_mean,-0.42);
  double p_shellform = 3.1e5 * pow(ejectaEnergy / 1e51, 0.94) * pow(n_mean, -0.13); // p_sf = m_sf*v_sf eq 9,11

  /*
      terminal momentum
          We picked the analytic form from Thornton as it matched high-resolution SN we ran the best.
  */ 

  double pTerminal;
  if (Z_Zsun > 0.01) {
    // Thornton, 1998, M_s * V_s, eqs 22, 23, 33, 34
    pTerminal = 1.6272e5 * pow(ejectaEnergy/1e51, 13./14.) * pow(n_mean, -0.25) * pow(Z_Zsun, -0.36);
  }

  else {
    pTerminal = 8.3619e5 * pow(ejectaEnergy/1e51, 13./14.) * pow(n_mean, -0.25);
  } 

  // fading radius of a supernova, using gas energy of the host cell and ideal gas approximations

  double T = ge[index] / d[index] * Tunit;
  double cSound = sqrt(cello::kboltz*T/cello::mass_hydrogen) / 1e5; // km/s
  double r_fade = std::max(66.0*pow(ejectaEnergy/1e51, 0.32)*pow(n_mean, -0.37)*pow(cSound/10, -2.0/5.0), 
                           CoolingRadius * 1.5);
  double fadeRatio = cellwidth/r_fade;

#ifdef DEBUG_FEEDBACK_STARSS
  CkPrintf("STARSS_FB: Fading: T = %e; Cs = %e; R_f = %e; fadeR = %f\n", T, cSound, r_fade, fadeRatio);
#endif

  double coupledMomenta = 0.0;
  double eKinetic = 0.0;

  double cw_eff = cellwidth; // effective cell width couples to farther than just dx.
                             // theres a lot of numerical fudge factors here because of that.
                             // the actual coupling is between 2-3 dx, depending on position within the cell

  double dx_eff = cw_eff / CoolingRadius;
  double fader = cw_eff / r_fade;

  double shellVelocity = 413.0 * pow(n_mean, 1.0 / 7.0) * pow(Z_Zsun, 3.0 / 14.0) 
                               * pow(coupledEnergy / 1e51, 1.0 / 14.0);
                               
  double ratio_pds = cw_eff/r_shellform;

  shellVelocity *=  ratio_pds > 1 ? pow(dx_eff, -7.0 / 3.0) : 1; //km/s

  float beta =  std::max( 1.0, shellVelocity / std::max(1.0, cSound));
  float nCritical = 0.0038 *(pow(n_mean * T / 1e4, 7.0/9.0) 
                           * pow(beta, 14.0/9.0))/(pow(ejectaEnergy/1e51, 1.0/9.0) 
                           * pow(std::max(0.0001, Z_Zsun), 1.0/3.0)); // n/cc

  float r_merge = std::max(151.0 * pow((ejectaEnergy/1e51) / 
                      beta / beta / n_mean / T * 1e4, 1./3.), 1.5 * r_fade); // pc
  float merger = cw_eff / r_merge;
  bool faded = fader > 1;

#ifdef DEBUG_FEEDBACK_STARSS
  CkPrintf(
     "STARSS_FB: RADII: cell = %e, free = %e, shellform = %e, cooling = %e, fade = %e t_3=%e, R_m=%e\n",
     cellwidth, r_free, r_shellform, CoolingRadius, r_fade, t3_sedov, r_merge);
#endif

  if (! winds)
    {   // this calculation for SNe only
        if (cw_eff < r_free){
          coupledMomenta = 0.0;// Thermal coupling only at free expansion limit.
        #ifdef DEBUG_FEEDBACK_STARSS
          CkPrintf("STARSS_FB: modifying free phase: p = %e\n", coupledMomenta);
        #endif
        }
        if (r_free < cw_eff && dx_eff <= 1){
          coupledMomenta = std::min(p_sedov, pTerminal*dx_eff);
        #ifdef DEBUG_FEEDBACK_STARSS       
          CkPrintf("STARSS_FB: Coupling Sedov-Terminal phase: p = %e (ps = %e, pt = %e, dxe = %e)\n", 
                 coupledMomenta, p_sedov, pTerminal, dx_eff);
        #endif
        }
        if (dx_eff > 1){   
          coupledMomenta = pTerminal/ sqrt(std::min(1.5, dx_eff));
        #ifdef DEBUG_FEEDBACK_STARSS
          CkPrintf("STARSS_FB: Coupling Terminal phase: p = %e; dx_eff = %e\n", coupledMomenta, dx_eff);
        #endif
        }
        if (fader > 1 && enzo_config->method_feedback_fade_SNR){ 
            // high coupling during the fading regime leads to SNRs on the root-grid in 6-level AMR simulations!
            coupledMomenta = pTerminal * (1.0 - tanh(pow(fader * merger, 2.5)));
        #ifdef DEBUG_FEEDBACK_STARSS 
          CkPrintf("STARSS_FB: Coupling Fading phase: p = %e\n", coupledMomenta);
        #endif
        }
        // critical density to skip snowplough (remnant combines with ISM before radiative phase); eq 4.9 cioffi 1988
      #ifdef DEBUG_FEEDBACK_STARSS
        CkPrintf("STARSS_FB: Checking critical density metric..." 
                 "(n_mean = %e; N_Crit = %e; factors: %e %e %e; beta = %e/%e == %e; rmerge = %e)\n", 
                  n_mean, nCritical, pow(n_mean * T / 1e4, 7.0/9.0), 
                  pow(ejectaEnergy/1e51, 1.0/9.0), pow(fz, 1.0/3.0), 
                  shellVelocity , cSound, beta, r_merge);
      #endif

        if (n_mean <= 10.0 * nCritical){ // in high-pressure, low nb, p_t doesn't hold since there is essentially no radiative phase.
                                         // thermal energy dominates the evolution (Tang, 2005, doi 10.1086/430875 )
                                         // We inject 100% thermal energy to simulate this recombining with the ISM
                                         // and rely on the hydro and the thermal radiation to arrive at the right solution
          coupledMomenta = coupledMomenta * (1.0-tanh(pow(1.45*nCritical/n_mean, 6.5)));
        #ifdef DEBUG_FEEDBACK_STARSS
          CkPrintf("STARSS_FB: Adjusting for high-pressure low-n phase (thermal coupling: Nc = %e): p = %e\n", 
                   nCritical, coupledMomenta);
        #endif
        }
            
        if (T > 1e6 && coupledMomenta > 1e5){
          #ifdef DEBUG_FEEDBACK_STARSS
            CkPrintf("STARSS_FB: Coupling high momenta to very hot gas!! (p= %e, T= %e, n_c = %e)\n", coupledMomenta, T, nCritical);
          #endif
        }
    } // endif no winds

    else { // if winds
      coupledMomenta = sqrt(ejectaMass*cello::mass_solar * 0.5 * ejectaEnergy)/cello::mass_solar/1e5;
    }
  #ifdef DEBUG_FEEDBACK_STARSS
    CkPrintf("STARSS_FB: Calculated p = %e (sq_fact = %e; p_f = %e; p_t = %e; mcell = %e; mcpl = %e)\n",
       coupledMomenta, (d_mean * cell_volume/cello::mass_solar) / ejectaMass * 63, p_free, 
       pTerminal, d_mean * cell_volume /cello::mass_solar, ejectaMass/27.0);  
  #endif

  /*
    If resolution is in a range comparable to Rcool and
    Analytic SNR shell mass is on, adjust the shell mass
    upper range of applicability for shell mass is determined by
    local average gas velocity (v_shell ~ c_s = no shell)
  */

  double centralMass = 0.0;
  double centralMetals = 0.0;
  double maxEvacFraction = 0.75; // TODO: make this a parameter
  double shellMass = 0.0;
   
  if (coupledEnergy > 0 && AnalyticSNRShellMass && !winds) 
  {
    if (dx_eff < 1) 
      shellMass = std::min(1e8, coupledMomenta / shellVelocity); //Msun

    else if (dx_eff < 1) {
      if (Z_Zsun > 0.01)
        shellMass = 1.41e4*pow(ejectaEnergy/1e51, 6./7.) * pow(d_mean, -0.24) * pow(Z_Zsun, -0.27);
      else
        shellMass = 4.89e4*pow(ejectaEnergy/1e51, 6./7.) * pow(d_mean, -0.24);
    }
   
    if (shellMass > 0)
    {
      /* cant let host cells evacuate completely! Will eventually subtract shell mass from host cell,
         so need to make sure that we don't end up with negative masses.
         Shell mass will be evacuated from central cells by CIC a negative mass,
         so have to check that the neighbors can handle it too*/
      for (int ix_ = ix-1; ix_ <= ix+1; ix_++) {
        for (int iy_ = iy-1; iy_ <= iy+1; iy_++) {
          for (int iz_ = iz-1; iz_ <= iz+1; iz_++) {
            int flat = INDEX(ix_,iy_,iz_,mx,my);
            if (flat < size) {
              // only record if this cell would've been touched by CIC on star particle 
              // (estimated as the "blast interior")
              centralMass += d[flat];
              centralMetals += mf[flat]; // TODO: original code had mass += density here for some reason 

            } //endif flat < size 
          }
        }
      } // endfor ix_
      
      centralMass *= munit / cello::mass_solar; // put into solar units
  
      if (shellMass > maxEvacFraction * centralMass) {
      #ifdef DEBUG_FEEDBACK_STARSS
        CkPrintf("STARSS: Shell mass too high for host cells: Rescaling %e -> %e\n", 
                 shellMass, maxEvacFraction * centralMass);
      #endif
        shellMass = maxEvacFraction * centralMass;
      }  
    } // endif shellMass > 0
      
  } // endif coupledEnergy >0 && AnalyticSNRShellMass && !winds
  
  double shellMetals = std::min(maxEvacFraction*centralMetals * munit/cello::mass_solar, 
                                Z_Zsun * cello::metallicity_solar * shellMass); //Msun

#ifdef DEBUG_FEEDBACK_STARSS
  if (AnalyticSNRShellMass) {
    CkPrintf("STARSS_FB: Shell_m = %e Shell_z = %e shell_V= %e P = %e M_C = %e\n",
             shellMass, shellMetals, shellVelocity, coupledMomenta, centralMass);    
  }  
#endif

  double coupledMass = shellMass + ejectaMass;
     
  eKinetic = coupledMomenta * coupledMomenta / (2.0 *(coupledMass)) * cello::mass_solar * 1e10;
  if (eKinetic > (nSNII+nSNIA) * 1e51 && !winds){
  #ifdef DEBUG_FEEDBACK_STARSS
    CkPrintf("STARSS_FB: Rescaling high kinetic energy %e -> ", eKinetic);
  #endif
    coupledMomenta = sqrt(2.0 * (coupledMass*cello::mass_solar) * ejectaEnergy)/cello::mass_solar/1e5;
    eKinetic = coupledMomenta * coupledMomenta / (2.0 *(coupledMass) * cello::mass_solar) * cello::mass_solar * 1e10;
  #ifdef DEBUG_FEEDBACK_STARSS        
    CkPrintf("STARSS_FB:  %e; new p = %e\n", eKinetic, coupledMomenta);
  #endif
    }
    
#ifdef DEBUG_FEEDBACK_STARSS
  CkPrintf("STARSS_FB: Ekinetic = %e Mass = %e\n",
           eKinetic, d_mean * pow(lunit * hx, 3) / cello::mass_solar);
#endif
    if (eKinetic > 1e60 && winds)
    {
    #ifdef DEBUG_FEEDBACK_STARSS
      CkPrintf("STARSS_FB: winds Ekinetic = %e Mass = %e\n",
                eKinetic, d_mean * pow(lunit * hx, 3) / cello::mass_solar);
    #endif
      ERROR("EnzoMethodFeedbackSTARSS::deposit_feedback()","winds Ekinetic > reasonability!\n");
    }

    if (eKinetic > 1e60 && !winds)
    {
    #ifdef DEBUG_FEEDBACK_STARSS
      CkPrintf("STARSS_FB: Ekinetic = %e Mass = %e\n",
               eKinetic, d_mean * pow(lunit * hx, 3) / cello::mass_solar);
    #endif
      ERROR("EnzoMethodFeedbackSTARSS::deposit_feedback()","SNE Ekinetic > reasonability!\n");
    }

    double coupledGasEnergy = std::max(ejectaEnergy - eKinetic, 0.0);
  #ifdef DEBUG_FEEDBACK_STARSS
    CkPrintf("STARSS_FB: Coupled Gas Energy = %e\n", coupledGasEnergy);
  #endif
    if (dxRatio > 1.0 && !winds){ // if we apply this reduction to winds, 
                                  // then there is literally *no* effect, 
                                  // even at Renaissance resolution.
        coupledGasEnergy = (coupledGasEnergy * pow(dxRatio, -6.5));
    #ifdef DEBUG_FEEDBACK_STARSS
      CkPrintf("STARSS_FB: Reducing gas energy... GE = %e\n", coupledGasEnergy);
    #endif
    }

  double coupledMetals = 0.0, SNIAmetals = 0.0, SNIImetals = 0.0;//, P3metals = 0.0;
  if (winds) coupledMetals = ejectaMetals; // winds only couple to metals
  if (AnalyticSNRShellMass && !winds) coupledMetals += shellMetals;
  
  ejectaMetals += nSNIA * 1.4; // add ejecta from Type Ia SNe 
  ejectaMetals += nSNII * (1.91 + 0.0479 * std::max(starZ, 1.65)); // add ejecta from Type II

  coupledMetals += ejectaMetals;

#ifdef DEBUG_FEEDBACK_STARSS
  CkPrintf("STARSS_FB: Coupled Metals: %e %e %e\n", ejectaMetals, shellMetals, coupledMetals);
#endif

  if (!winds) coupledEnergy = std::min((nSNII + nSNIA) * 1e51, eKinetic);

  // Subtract shell mass from the central cells
  double minusRho=0, minusZ=0, msubtracted=0;
  double remainMass=shellMass*cello::mass_solar/munit, remainZ = shellMetals*cello::mass_solar/munit;
  bool massiveCell;
  if (shellMass > 0 && AnalyticSNRShellMass)
  {
    double zsubtracted=0;
    massiveCell=false;
    msubtracted=0;
    for (int ix_ = ix-1; ix_ <= ix+1; ix_++) {
      for (int iy_ = iy-1; iy_ <= iy+1; iy_++) {
        for (int iz_ = iz-1; iz_ <= iz+1; iz_++) {
          int flat = INDEX(ix_,iy_,iz_,mx,my);
          
          if ( !(flat > 0 && flat < size) ) continue; // TODO: do I even need this?
          
          double dpre = d[flat];
          double zpre = mf[flat];
          double pre_z_frac = zpre / dpre;

        #ifdef DEBUG_FEEDBACK_STARSS
          CkPrintf("STARSS: Baryon Prior: d_Msun = %e, mc = %e, ms = %e, m_z = %e, z = %e\n",
                   d[flat] * munit/cello::mass_solar, centralMass, shellMass, shellMetals, pre_z_frac);
        #endif
          d[flat] = std::max(dpre - remainMass/27.0, (1.0-maxEvacFraction)*dpre);
          minusRho    += dpre - d[flat];
          msubtracted += dpre - d[flat];
          
          mf[flat] = std::max(tiny_number, zpre - remainZ/27.0);
          minusZ      += zpre - mf[flat];
          zsubtracted += zpre - mf[flat];

        } // endfor iz_
      } // endfor iy_
    } // endfor ix_
    remainMass -= msubtracted;
    remainZ    -= zsubtracted;
  } // endif shellMass > 0 && AnalyticSNRShellMass

  minusRho *= munit/cello::mass_solar; //Msun
  minusZ   *= munit/cello::mass_solar; //Msun

  double minusM = minusRho;

  if (minusM != coupledMass - ejectaMass && shellMass > 0 && AnalyticSNRShellMass)
  {
  #ifdef DEBUG_FEEDBACK_STARSS
    CkPrintf("STARSS_FB: Of %e, only subtracted %e; rescaling the coupling mass\n", shellMass,
              minusRho);
  #endif
  
    double oldcouple = coupledMass;
    coupledMass = minusM + ejectaMass;
    coupledMetals = minusZ + ejectaMetals; // metal_density in Msun/code_volume

    // if we're in here, we found an inconsistent place where the mass within 
    // cannot support the expected shell.
    // the only real choice, to keep things from exploding in a bad way, 
    // is to rescale the momentum accordingly.
    // If this is triggered, dont expect the terminal momenta relations to hold up.

    if (coupledMomenta*coupledMomenta / 
        (2.0*coupledMass*cello::mass_solar)*cello::mass_solar * 1e10 > 1e51) {
      coupledMomenta = std::min(coupledMomenta, sqrt(2.0 * (coupledMass*cello::mass_solar)
                              * ejectaEnergy)/cello::mass_solar/1e5);

    #ifdef DEBUG_FEEDBACK_STARSS
      CkPrintf("STARSS_FB: rescaled momentum to %e (est KE = %e)\n", 
                coupledMomenta, coupledMomenta*coupledMomenta / (2*coupledMass) 
              * cello::mass_solar * 1e10);
    #endif
    }
  
  } // endif minusM != coupledMass - ejectaMass

#ifdef DEBUG_FEEDBACK_STARSS
  CkPrintf("STARSS_FB: Before unit conversions -- coupledEnergy = %e; coupledGasEnergy = %e;\n"
           "           coupledMass = %e; coupledMetals = %e; coupledMomenta = %e\n",
                       coupledEnergy, coupledGasEnergy, coupledMass, coupledMetals,
                       coupledMomenta);
  CkPrintf("STARSS_FB: Units -- munit = %e; vunit = %e; Tunit = %e\n"
           "                    eunit = %e; rhounit = %e; lunit = %e\n",
                       munit, vunit, Tunit, eunit, rhounit, lunit);
#endif
  //put everything back into code units before CIC
  //these values correspond to TOTAL (energy, mass, etc.)
  //over all the coupling particles.
  //Note that coupledMass and coupledMetals are densities.

  coupledEnergy /= eunit; // coupledEnergy is kinetic energy at this point
  coupledGasEnergy /= eunit;
  coupledMass /= (munit/cello::mass_solar);
  coupledMetals /= (munit/cello::mass_solar);
  coupledMomenta /= (munit/cello::mass_solar * vunit / 1e5 / sqrt(nCouple) ); 

  coupledEnergy += coupledGasEnergy;

  // Create SN coupling particles
//#ifdef DEBUG_FEEDBACK_STARSS
//  CkPrintf("STARSS_FB: Creating coupling particles...\n");
//#endif

  //this->createCouplingParticles(enzo_block, nCouple, coupledEnergy, coupledGasEnergy,
  //                             coupledMass, coupledMetals, coupledMomenta,
  //                              CloudParticlePositionX,CloudParticlePositionY,CloudParticlePositionZ,  
  //                              CloudParticleVectorX,CloudParticleVectorY,CloudParticleVectorZ);


  // transform specific energy to energy before deposition

//  for (int i=0; i<mx*my*mz; i++)
//  {
//    te[i] *= (d[i]*hx*hy*hz);
//    ge[i] *= (d[i]*hx*hy*hz);
//  }

  enzo_float coupledMass_list[nCouple], coupledMetals_list[nCouple];
  enzo_float coupledMomenta_x[nCouple], coupledMomenta_y[nCouple], coupledMomenta_z[nCouple];
  enzo_float coupledEnergy_list[nCouple], coupledGasEnergy_list[nCouple];

  std::fill_n(coupledMass_list  ,nCouple,  coupledMass/nCouple);
  std::fill_n(coupledMetals_list,nCouple,coupledMetals/nCouple);

  for (int n=0; n<nCouple; n++)
  {
    coupledMass_list[n] = coupledMass/nCouple;
    coupledMetals_list[n] = coupledMetals/nCouple;
    coupledMomenta_x[n] = coupledMomenta/nCouple*CloudParticleVectorX[n];
    coupledMomenta_y[n] = coupledMomenta/nCouple*CloudParticleVectorY[n];
    coupledMomenta_z[n] = coupledMomenta/nCouple*CloudParticleVectorZ[n];
    coupledEnergy_list[n] = coupledEnergy/nCouple;
    coupledGasEnergy_list[n] = coupledGasEnergy/nCouple;
  }

#ifdef DEBUG_FEEDBACK_STARSS
  double sumEnergy=0.0, sumMomenta=0.0, sumMass=0.0, sumMetals=0.0, sumInternal=0.0; 
  for (int i=0; i<mx*my*mz; i++) 
  {
    sumEnergy += te[i];
    sumInternal += ge[i];
    sumMomenta += vx[i]*vx[i];
    sumMomenta += vy[i]*vy[i];
    sumMomenta += vz[i]*vz[i];
    sumMass += d[i];
    sumMetals += mf[i]; 
  }
 
  sumMomenta = sqrt(sumMomenta);

  CkPrintf("STARSS_FB: Total before deposition -- Energy = %e;\n"
           "                                     Internal energy = %e;\n" 
           "                                     P = %e;\n"
           "                                     Mass (actually density) = %e\n"
           "                                     Metal Mass (actually density) = %e\n",
            sumEnergy, sumInternal, sumMomenta, sumMass, sumMetals);
  CkPrintf("STARSS_FB: source cell internal_energy = %e; total_energy = %e\n", ge[index], te[index]);
  CkPrintf("STARSS_FB: CIC depositing coupling particles...\n");
  CkPrintf("STARSS_FB: Depositing (code units) -- coupledEnergy = %e; coupledGasEnergy = %e;\n"
           "           coupledMass = %e; coupledMetals = %e; coupledMomenta = %e\n",
                       coupledEnergy, coupledGasEnergy, coupledMass, coupledMetals,
                       coupledMomenta);


#endif

  enzo_float left_edge[3] = {xm-gx*hx, ym-gy*hy, zm-gz*hz};
  //enzo_float left_edge[3] = {xm, ym, zm};
  // CIC deposit coupling particles
  //
  enzo_float * d_dep  = (enzo_float *) field.values("density_deposit");
  enzo_float * te_dep = (enzo_float *) field.values("total_energy_deposit");
  enzo_float * ge_dep = (enzo_float *) field.values("internal_energy_deposit");
  enzo_float * mf_dep = (enzo_float *) field.values("metal_density_deposit");
  enzo_float * vx_dep = (enzo_float *) field.values("velocity_x_deposit");
  enzo_float * vy_dep = (enzo_float *) field.values("velocity_y_deposit");
  enzo_float * vz_dep = (enzo_float *) field.values("velocity_z_deposit");

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMass_list, d_dep, &left_edge,
   &mx, &my, &mz, &hx, &A);

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMomenta_x, vx_dep, &left_edge,
   &mx, &my, &mz, &hx, &A);

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMomenta_y, vy_dep, &left_edge,
   &mx, &my, &mz, &hx, &A);

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMomenta_z, vz_dep, &left_edge,
   &mx, &my, &mz, &hx, &A);

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMetals_list, mf_dep, &left_edge,
   &mx, &my, &mz, &hx, &A);


  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledEnergy_list, te_dep, &left_edge,
   &mx, &my, &mz, &hx, &A);

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledGasEnergy_list, ge_dep, &left_edge,
   &mx, &my, &mz, &hx, &A);


  //copy values for active zones (include ghost zones if > 1 block???)
  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){
        int i = INDEX(ix,iy,iz,mx,my);
  
         d[i] += d_dep[i];
        mf[i] += mf_dep[i]; 
        double cell_mass = d[i]*hx*hy*hz;
        
        te[i] += te_dep[i] / cell_mass;
        ge[i] += ge_dep[i] / cell_mass;
        vx[i] += vx_dep[i];
        vy[i] += vy_dep[i];
        vz[i] += vz_dep[i];
 
      }
    }
  }

#ifdef DEBUG_FEEDBACK_STARSS
  sumEnergy=sumMomenta=sumMass=sumMetals=sumInternal=0.0; 
  for (int i=0; i<mx*my*mz; i++) 
  {
    sumEnergy += te[i]; // *d[i]*cell_volume/enzo_units->volume();
    sumInternal += ge[i];// *d[i]*cell_volume/enzo_units->volume();
    sumMomenta += vx[i]*vx[i];
    sumMomenta += vy[i]*vy[i];
    sumMomenta += vz[i]*vz[i];
    sumMass += d[i];
    sumMetals += mf[i];
  }
 
  sumMomenta = sqrt(sumMomenta);

  CkPrintf("STARSS_FB: Total after deposition -- Energy = %e;\n"
           "                                     Internal energy = %e;\n" 
           "                                     P = %e;\n"
           "                                     Mass (actually density) = %e;\n"
           "                                     Metal Mass (actually density) = %e\n",
            sumEnergy, sumInternal, sumMomenta, sumMass, sumMetals);
  //CkPrintf("STARSS_FB: Deleting coupling particles...\n");
#endif


  // transform back to specific energy
//  for (int i=0; i<mx*my*mz; i++)
//  {
//    te[i] /= (d[i]*hx*hy*hz);
//    ge[i] /= (d[i]*hx*hy*hz);
//  }


  //this->deleteCouplingParticles(enzo_block);

//#ifdef DEBUG_FEEDBACK_STARSS
  //CkPrintf("STARSS_FB: Updating velocity fields to account for momentum deposition...\n");
//#endif


  // transform velocities back to "lab" frame
  // convert velocity (actually momentum at the moment) field back to velocity 
  this->transformComovingWithStar(d,vx,vy,vz,up,vp,wp,mx,my,mz, -1);
  this->transformComovingWithStar(d_dep,vx_dep,vy_dep,vz_dep,up,vp,wp,mx,my,mz, -1);

  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){
        int i = INDEX(ix,iy,iz,mx,my);
  
         d_dep[i] = 0.0;
        mf_dep[i] = 0.0;
        te_dep[i] = 0.0;
        ge_dep[i] = 0.0;
        vx_dep[i] = 0.0;
        vy_dep[i] = 0.0;
        vz_dep[i] = 0.0;

 
      }
    }
  }



// TODO: Update multispecies fields within deposited cells to be consistent with added mass

#ifdef DEBUG_FEEDBACK_STARSS
  CkPrintf("STARSS_FB: Refreshing fields...\n");
#endif

//  if (enzo_config->field_prolong == "enzo") {
    // TODO: scale velocity fields by mass before refreshing. How 
    //       do I make sure neighboring blocks' velocity fields
    //       are also converted to momentum before refresh??
    //
    //       Refresh is only needed to catch rare-ish edge cases of stars exploding right
    //       on the border of two blocks. Assuming the gas mass of the ghost cell after deposition
    //       and its corresponding active cell in the neighboring black are close to eachother, not converting
    //       velocity->momentum and specific energy->energy shouldn't create a *huge* error. This algorithm only
    //       deposits ~Msun to a cell, so should be okay-ish. 
    //
    // currently done automatically for linear interpolation
    // but not EnzoProlong.
    //
    // TODO: Remove this if/when automatic scaling is
    //       added to EnzoProlong
    
//  }


  // refresh
//#endif
}



double EnzoMethodFeedbackSTARSS::timestep (Block * block) throw()
{
  // In general this is not needed, but could imagine putting timestep
  // limiters in situations where, for example, one would want
  // dt < star_lifetime (or something like that), especially if
  // important things happen throughout the star's lifetime.
  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();

  double dtStar = std::numeric_limits<double>::max();
  if (block->level() >= sf_minimum_level_){
    const double pSNmax = 0.0005408 * enzo_config->method_star_maker_minimum_star_mass *
                          block->dt() * enzo_units->time() / cello::Myr_s * 1.25;
    if (pSNmax > 1.0) dtStar = block->dt() * 1.0 / pSNmax;
  }

  return dtStar;
}


// #endif
