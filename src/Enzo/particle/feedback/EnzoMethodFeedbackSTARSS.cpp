
/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodFeedbackSTARSS.cpp
/// @author     Will Hicks (whicks@ucsd.edu)
//              Andrew Emerick (aemerick11@gmail.com)
//             
/// @date
/// @brief  Implements the STARSS model for stellar feedback
///         as implemented into Enzo as the MechStars model
///         by Azton Wells. This is a direct copy/paste of the
///         MechStars methods in Enzo, adapted for Enzo-E.


#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/particle/particle.hpp"

#include <time.h>

//#ifdef NOTDEFINED // for now... since not done coding

//#define DEBUG_FEEDBACK_STARSS
//#define DEBUG_FEEDBACK_STARSS_SN

// =============================================================================
// splice these off to a different file (later)
// TODO: Maybe create EnzoStarParticle class and add rate-calculating functions there?

int EnzoMethodFeedbackSTARSS::determineSN(double age_Myr, int* nSNII, int* nSNIA,
                double pmass_Msun, double tunit, float dt){

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
    if (this->supernovae_ && NEvents < 0)
    {
        /* age-dependent rates */
        if (age_Myr < 3.401)
        {
            RII = 0.0;
            RIA = 0.0;
        }
        if (3.401 <= age_Myr && age_Myr < 10.37)
        {
            RII = 5.408e-4;
            RIA = 0.0;
        }
        if (10.37 <= age_Myr && age_Myr < 37.53)
        {
            RII = 2.516e-4;
            RIA = 0.0;
        }
        if (37.53 <= age_Myr)
        {
            RII = 0.0;
            RIA = 5.3e-8+1.6e-5*exp(-0.5*pow((age_Myr-50.0)/10.0, 2));
        }
        /* rates -> probabilities */
        if (RII > 0){
            PII = RII * pmass_Msun / enzo_constants::Myr_s *tunit*dt;
            double random = double(mt())/double(mt.max());
            if (PII > 1.0){
                if (this->unrestricted_sn_) {
                    int round = (int)PII;
                    *nSNII = round;  // number of Type II SNe that will definitely go off
                    PII -= round; // probability for setting off one more supernova
                }
                else {
                  // Restrict particles to at most one supernova per timestep
                  // if unrestricted_sn = false
                  *nSNII = 1;
                  PII = 0; // set probability for setting off another supernova to zero
                }
            }
           
            if (random <= PII){
                *nSNII += 1;
            }
        }

        #ifdef DEBUG_FEEDBACK_STARSS_SN
          CkPrintf("MethodFeedbackSTARSS::determineSN() -- pmass_Msun = %f; age = %f; RII = %f; RIA = %f\n",
                    pmass_Msun, age_Myr, RII, RIA);
        #endif

        if (RIA > 0){
            PIA = RIA*pmass_Msun / enzo_constants::Myr_s *tunit*dt;
            float random = float(rand())/float(RAND_MAX);

            if (PIA > 1.0)
            {
                if (this->unrestricted_sn_) {
                    int round = int(PIA);
                    *nSNIA = round; // number of Type Ia SNe that will definitely go off
                    PIA -= round; // probability for setting off one more supernova
                }
                else {
                    *nSNIA = 1;
                    PIA = 0;
                }
            }
            
            if (!this->unrestricted_sn_ && *nSNII > 0) {
                // for unrestricted_sn = false, only set off a Type IA supernova
                // if nSNII == 0
                *nSNIA = 0;
                PIA = 0;
            }

            if (random < PIA){
                *nSNIA += 1;
            }
        }

        #ifdef DEBUG_FEEDBACK_STARSS_SN
          CkPrintf("MethodFeedbackSTARSS::determineSN() -- PII = %f; PIA = %f\n", PII, PIA);
        #endif
    }
 
        return 1;
}

int EnzoMethodFeedbackSTARSS::determineWinds(double age_Myr, double * eWinds, double * mWinds, double * zWinds,
                      double pmass_Msun, double metallicity_Zsun, double tunit, double dt) 
    {
    // age in Myr

    const EnzoConfig * enzo_config = enzo::config();
    bool oldEnough = (age_Myr < 0.0001)?(false):(true);
    double windE = 0,  wind_mass_solar = 0, windZ = 0.0;
    double wind_factor = 0.0;
    double e_factor = 0.0;

    if (pmass_Msun > 11 && oldEnough){

        if (0.001 < age_Myr && age_Myr < 1.0){
            wind_factor = 4.763 * std::min((0.01 + metallicity_Zsun), 1.0);
        }
        if (1 <= age_Myr && age_Myr < 3.5){
            wind_factor = 4.763*std::min(0.01+metallicity_Zsun, 1.0)* 
                pow(age_Myr, 1.45+0.08*std::min(log(metallicity_Zsun), 1.0));
        }
        if (3.5 <= age_Myr && age_Myr < 100){
            wind_factor = 29.4*pow(age_Myr/3.5, -3.25)+0.0042;
        }
        if (age_Myr < 100){
            double d = powl(1+age_Myr/2.5, 1.4);
            double a50 = powl(double(age_Myr)/10.0, 5.0);
            e_factor = 5.94e4 / d + a50 +4.83;
        }
        if (100 <= age_Myr){
            e_factor = 4.83;
            wind_factor = 0.42*pow(age_Myr/1000, -1.1)/(19.81/log(age_Myr));
        }
        wind_mass_solar = pmass_Msun * wind_factor; // Msun/Gyr
        wind_mass_solar = wind_mass_solar*dt*tunit/(1e3 * enzo_constants::Myr_s); // Msun

        if (wind_mass_solar > pmass_Msun){
            CkPrintf("Winds too large Mw = %e, Mp = %e age=%f, Z = %e\n",
                wind_mass_solar, pmass_Msun, age_Myr, metallicity_Zsun);
            wind_mass_solar = 0.125*pmass_Msun; // limit loss to huge if necessary.
        }
        windZ = std::max(enzo_constants::metallicity_solar, 0.016+0.0041*std::max(metallicity_Zsun, 1.65)+0.0118)*wind_mass_solar;
        windE = e_factor * 1e12 * wind_mass_solar;


    *mWinds = wind_mass_solar;
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
    // NOTE: This transforms the velocity field into a momentum density
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
    // back to "lab" frame. Convert momentum density field back to velocity
    for (int ind = 0; ind<size; ind++) {
      //if (density[ind] <= 10*1e-20) continue;
      if (density[ind] == 0) continue;
      double mult = 1/density[ind];
      velocity_x[ind] = velocity_x[ind]*mult + up;
      velocity_y[ind] = velocity_y[ind]*mult + vp;
      velocity_z[ind] = velocity_z[ind]*mult + wp;
    }
  }

}


double EnzoMethodFeedbackSTARSS::Window(double xd, double yd, double zd, double width) const throw()
{
    float wx = 0, wy = 0, wz = 0;
    if (std::abs(xd) <= width)
            wx = 1.0 - std::abs(xd)/width;
    if (std::abs(yd) <= width)
            wy = 1.0 - std::abs(yd)/width;
    if (std::abs(zd) <= width)
            wz = 1.0 - std::abs(zd)/width;
    return wx * wy * wz;
}

extern "C" void FORTRAN_NAME(cic_deposit)
  ( enzo_float (*px)[26], enzo_float (*py)[26], enzo_float (*pz)[26], const int * rank,
    const int * nCouple, enzo_float (*mass)[26], enzo_float * field, enzo_float (*left_edge)[3],
    int * mx, int * my, int * mz, double * hx, const double * cloudsize );
// =============================================================================
// =============================================================================
// =============================================================================

EnzoMethodFeedbackSTARSS::EnzoMethodFeedbackSTARSS(ParameterGroup p)
  : Method()
  , ir_feedback_(-1)
{
  FieldDescr * field_descr = cello::field_descr();
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  ASSERT("EnzoMethodFeedbackSTARSS::EnzoMethodFeedbackSTARSS",
         "untested without dual-energy formalism",
         ! enzo::fluid_props()->dual_energy_config().is_disabled());

  // required fields
  cello::define_field("density");
  cello::define_field("pressure");
  cello::define_field("total_energy");
  cello::define_field("internal_energy");
  cello::define_field("velocity_x");
  cello::define_field("velocity_y");
  cello::define_field("velocity_z");
  cello::define_field("metal_density");

  
  cello::define_field_in_group("metal_density","color");

  // Initialize refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  // previously a parameter was parsed called "Method:feedback:min_level",
  // but that never got used!

  supernovae_              = p.value_logical("supernovae",true);
  unrestricted_sn_         = p.value_logical("unrestricted_sn",true);
  stellar_winds_           = p.value_logical("stellar_winds",true);
  radiation_               = p.value_logical("radiation",true);
  analytic_SNR_shell_mass_ = p.value_logical("analytic_SNR_shell_mass",true);
  fade_SNR_                = p.value_logical("fade_SNR",true);

  // initialize NEvents parameter (mainly for testing). Sets off 'NEvents'
  // supernovae, with at most one supernova per star particle per cycle.
  NEvents                  = p.value_integer("NEvents",-1);

  // Initialize temporary fields
  i_d_dep  = cello::field_descr()->insert_temporary();
  i_te_dep = cello::field_descr()->insert_temporary();
  i_ge_dep = cello::field_descr()->insert_temporary();
  i_mf_dep = cello::field_descr()->insert_temporary();
  i_vx_dep = cello::field_descr()->insert_temporary();
  i_vy_dep = cello::field_descr()->insert_temporary();
  i_vz_dep = cello::field_descr()->insert_temporary();
  i_d_shell= cello::field_descr()->insert_temporary();

  i_d_dep_a  = cello::field_descr()->insert_temporary();
  i_te_dep_a = cello::field_descr()->insert_temporary();
  i_ge_dep_a = cello::field_descr()->insert_temporary();
  i_mf_dep_a = cello::field_descr()->insert_temporary();
  i_vx_dep_a = cello::field_descr()->insert_temporary();
  i_vy_dep_a = cello::field_descr()->insert_temporary();
  i_vz_dep_a = cello::field_descr()->insert_temporary();
  i_d_shell_a= cello::field_descr()->insert_temporary();

  // Deposition across grid boundaries is handled using refresh with set_accumulate=true.
  // The set_accumulate flag tells Cello to include ghost zones in the refresh operation,
  // and adds the ghost zone values from the "src" field to the corresponding active zone
  // values in the "dst" field.
 
  // I'm currently using a set of two initially empty fields for each field that SNe directly affect.
  // CIC deposition directly modifies the *_deposit fields (including ghost zones). The ghost zone
  // values are then sent to the *deposit_accumulate fields during the refresh operation. Values are then
  // copied back to the original field. 
  
  ir_feedback_ = add_refresh_();
  cello::simulation()->refresh_set_name(ir_feedback_,name()+":add");
  Refresh * refresh_fb = cello::refresh(ir_feedback_); 

  refresh_fb->set_accumulate(true);

  refresh_fb->add_field_src_dst
    (i_d_dep, i_d_dep_a);
  refresh_fb->add_field_src_dst
    (i_te_dep, i_te_dep_a);
  refresh_fb->add_field_src_dst
    (i_ge_dep, i_ge_dep_a);
  refresh_fb->add_field_src_dst
    (i_mf_dep, i_mf_dep_a);
  refresh_fb->add_field_src_dst
    (i_vx_dep, i_vx_dep_a);
  refresh_fb->add_field_src_dst
    (i_vy_dep, i_vy_dep_a);
  refresh_fb->add_field_src_dst
    (i_vz_dep, i_vz_dep_a);
  refresh_fb->add_field_src_dst
    (i_d_shell, i_d_shell_a);
 
  refresh_fb->set_callback(CkIndex_EnzoBlock::p_method_feedback_starss_end());

  return;
}

void EnzoMethodFeedbackSTARSS::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | supernovae_;
  p | unrestricted_sn_;
  p | stellar_winds_;
  p | radiation_;
  p | analytic_SNR_shell_mass_;
  p | fade_SNR_;
  p | NEvents;
  p | ir_feedback_;

  p | i_d_dep;
  p | i_te_dep;
  p | i_ge_dep;
  p | i_mf_dep;
  p | i_vx_dep;
  p | i_vy_dep;
  p | i_vz_dep;
  p | i_d_shell;

  p | i_d_dep_a;
  p | i_te_dep_a;
  p | i_ge_dep_a;
  p | i_mf_dep_a;
  p | i_vx_dep_a;
  p | i_vy_dep_a;
  p | i_vz_dep_a;
  p | i_d_shell_a;

  return;
}

void EnzoMethodFeedbackSTARSS::compute (Block * block) throw()
{

  if (block->is_leaf()){
    this->compute_(block);
  }

  else {
    block->compute_done();
  }
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
  
  enzo_float * d_dep  = (enzo_float *) field.values(i_d_dep);
  enzo_float * te_dep = (enzo_float *) field.values(i_te_dep);
  enzo_float * ge_dep = (enzo_float *) field.values(i_ge_dep);
  enzo_float * mf_dep = (enzo_float *) field.values(i_mf_dep);
  enzo_float * vx_dep = (enzo_float *) field.values(i_vx_dep);
  enzo_float * vy_dep = (enzo_float *) field.values(i_vy_dep);
  enzo_float * vz_dep = (enzo_float *) field.values(i_vz_dep);

  enzo_float * d_dep_a  = (enzo_float *) field.values(i_d_dep_a);
  enzo_float * te_dep_a = (enzo_float *) field.values(i_te_dep_a);
  enzo_float * ge_dep_a = (enzo_float *) field.values(i_ge_dep_a);
  enzo_float * mf_dep_a = (enzo_float *) field.values(i_mf_dep_a);
  enzo_float * vx_dep_a = (enzo_float *) field.values(i_vx_dep_a);
  enzo_float * vy_dep_a = (enzo_float *) field.values(i_vy_dep_a);
  enzo_float * vz_dep_a = (enzo_float *) field.values(i_vz_dep_a);

  enzo_float * d_shell   = (enzo_float *) field.values(i_d_shell);
  enzo_float * d_shell_a = (enzo_float *) field.values(i_d_shell_a);

  EnzoUnits * enzo_units = enzo::units();
  double cell_volume_code = hx*hy*hz;
  double cell_volume_cgs = cell_volume_code * enzo_units->volume();
  double rhounit = enzo_units->density();

  double maxEvacFraction = 0.75; // TODO: make this a parameter

  // multiply by this value to convert cell density (in code units)
  // to cell mass in Msun
  double rho_to_m = rhounit*cell_volume_cgs / enzo_constants::mass_solar;

  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){
        int i = INDEX(ix,iy,iz,mx,my);
        if (d_dep_a[i] != 0) { // if any deposition

          double d_old = d[i];
          
          // Could have a race condition here where if one particle updates the density of a cell,
          // that update won't get communicated to the other block until the end of the cycle.
          // In rare cases, this could result in a cell having a negative density because "centralMass"
          // and "centralMetals" in deposit_feedback() will be overpredicted going into the negative-mass
          // CiC for clearing out gas from the central cell. Catch this case by setting density
          // to (1-maxEvacFraction) * density[i] and metal_density to (1-maxEvacFraction) * metal_density[i] if
          // either go negative.

          if (d[i] + d_dep_a[i] < 0) {
            d[i] *= 1-maxEvacFraction;
          }
          else {
            d[i] += d_dep_a[i];
          }

          if (mf[i] + mf_dep_a[i] < 0) {
            mf[i] *= 1-maxEvacFraction;
          }
          else {
            mf[i] += mf_dep_a[i];
          }

          double d_new = d[i];
          double inv_dens_new = 1.0/d_new;

          double M_scale_tot = d_new / d_old;
          double M_scale_shell = d_shell_a[i]/d[i];

          // Here, te_dep_a and ge_dep_a are carrying "energy density" (not specific energy)
          // and vx_dep_a, vy_dep_a, and vy_dep_a are carrying velocity of the shell (not momentum density)
 
          te[i] = te[i] / M_scale_tot + std::max( (double) te_dep_a[i], 0.0) * inv_dens_new; 
          ge[i] = ge[i] / M_scale_tot + std::max( (double) ge_dep_a[i], 0.0) * inv_dens_new;
          vx[i] += vx_dep_a[i] * M_scale_shell;
          vy[i] += vy_dep_a[i] * M_scale_shell;
          vz[i] += vz_dep_a[i] * M_scale_shell;

          // rescale color fields to account for new densities
          EnzoMethodStarMaker::rescale_densities(enzo_block, i, M_scale_tot);
          // undo rescaling of metal_density field
          mf[i] /= M_scale_tot;

         }        
        
      }
    }
  }

  deallocate_temporary_(enzo_block);
  return;
}

void EnzoBlock::p_method_feedback_starss_end() 
{  
  EnzoMethodFeedbackSTARSS * method = static_cast<EnzoMethodFeedbackSTARSS*> (this->method());
  method->add_accumulate_fields(this);

  compute_done();
  return;
}

void EnzoMethodFeedbackSTARSS::compute_ (Block * block)
{

  //----------------------------------------------------
  // some constants here that might get moved to parameters or something else
  const float SNII_ejecta_mass_Msun = 10.5;
  const float SNIa_ejecta_mass_Msun = 1.4;
  const float z_solar          = enzo_constants::metallicity_solar; //Solar metal fraction (0.012)
  //-----------------------------------------------------

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  Particle particle = enzo_block->data()->particle();
  EnzoUnits * enzo_units = enzo::units();

  double munit = enzo_units->mass();
  double tunit = enzo_units->time();
  double lunit = enzo_units->length();

  double current_time  = block->state()->time();

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

  const int rank = cello::rank();

  // apply feedback depending on particle type
  // for now, just do this for all star particles

  int it = particle.type_index("star");

  // AE: Copying this from MechStars_FeedbackRoutine

  int numSN = 0; // counter of SN events
  int count = 0; // counter of particles

  allocate_temporary_(enzo_block);

  // initialize temporary fields as zero
  enzo_float * d_dep  = (enzo_float *) field.values(i_d_dep);
  enzo_float * te_dep = (enzo_float *) field.values(i_te_dep);
  enzo_float * ge_dep = (enzo_float *) field.values(i_ge_dep);
  enzo_float * mf_dep = (enzo_float *) field.values(i_mf_dep);
  enzo_float * vx_dep = (enzo_float *) field.values(i_vx_dep);
  enzo_float * vy_dep = (enzo_float *) field.values(i_vy_dep);
  enzo_float * vz_dep = (enzo_float *) field.values(i_vz_dep);

  enzo_float * d_dep_a  = (enzo_float *) field.values(i_d_dep_a);
  enzo_float * te_dep_a = (enzo_float *) field.values(i_te_dep_a);
  enzo_float * ge_dep_a = (enzo_float *) field.values(i_ge_dep_a);
  enzo_float * mf_dep_a = (enzo_float *) field.values(i_mf_dep_a);
  enzo_float * vx_dep_a = (enzo_float *) field.values(i_vx_dep_a);
  enzo_float * vy_dep_a = (enzo_float *) field.values(i_vy_dep_a);
  enzo_float * vz_dep_a = (enzo_float *) field.values(i_vz_dep_a);

  enzo_float * d_shell   = (enzo_float *) field.values(i_d_shell);
  enzo_float * d_shell_a = (enzo_float *) field.values(i_d_shell_a);


  for (int i=0; i<mx*my*mz; i++){
    d_dep [i] = 0;    
    te_dep[i] = 0;   
    ge_dep[i] = 0;
    mf_dep[i] = 0;
    vx_dep[i] = 0;
    vy_dep[i] = 0;
    vz_dep[i] = 0;
    d_shell[i] = 0;

    d_dep_a [i] = 0;
    te_dep_a[i] = 0;
    ge_dep_a[i] = 0;
    mf_dep_a[i] = 0;
    vx_dep_a[i] = 0;
    vy_dep_a[i] = 0;
    vz_dep_a[i] = 0;
    d_shell_a[i] = 0;
  }

  double cell_volume = hx*hy*hz;

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
  const int ia_lum = particle.attribute_index (it, "luminosity");

  const int dm = particle.stride(it, ia_m);
  const int dp = particle.stride(it, ia_x);
  const int dv = particle.stride(it, ia_vx);
  const int dl = particle.stride(it, ia_l);
  const int dc = particle.stride(it, ia_c);
  const int dmf = particle.stride(it, ia_mf);
  const int dsn = particle.stride(it, ia_sn);
  const int dlum = particle.stride(it, ia_lum);

  const int nb = particle.num_batches(it);

  for (int ib=0; ib<nb; ib++){
    enzo_float *px=0, *py=0, *pz=0, *pvx=0, *pvy=0, *pvz=0;
    enzo_float *plifetime=0, *pcreation=0, *pmass=0, *pmetal=0, *psncounter=0, *plum=0;

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

    plum = (enzo_float *) particle.attribute_array(it, ia_lum, ib);

    int np = particle.num_particles(it,ib);

    for (int ip=0; ip<np; ip++){
      int ipdp = ip*dp; // pos
      int ipdm = ip*dm; // mass
      int ipdv = ip*dv; // velocity
      int ipdl = ip*dl; // lifetime
      int ipdc = ip*dc; // creation time
      int ipdmf = ip*dmf; // metallicity
      int ipsn  = ip*dsn; // number of SNe counter
      int ipdlum = ip*dlum; // particle luminosity counter

      double pmass_solar = pmass[ipdm] * munit/enzo_constants::mass_solar;

      if (pmass[ipdm] > 0.0 && plifetime[ipdl] > 0.0){
        const double age = (current_time - pcreation[ipdc]) * enzo_units->time() / enzo_constants::Myr_s;
        count++; // increment particles examined here

        // compute coordinates of central feedback cell
        // this must account for ghost zones

        int ix       = (int) floor((px[ipdp] - xm) / hx + gx);
        int iy       = (int) floor((py[ipdp] - ym) / hy + gy);
        int iz       = (int) floor((pz[ipdp] - zm) / hz + gz);

        int index = INDEX(ix,iy,iz,mx,my); // 1D cell index of star position

        int nSNII = 0, nSNIa = 0;
        double SNMassEjected = 0.0, SNMetalEjected = 0.0;

        if (supernovae_){ 

          /* Determine number of SN events from rates (currently taken from Hopkins 2018) */

          determineSN(age, &nSNII, &nSNIa, pmass_solar,
                      tunit, block->state()->dt());

          numSN += nSNII + nSNIa;

          #ifdef DEBUG_FEEDBACK_STARSS
            if (nSNII + nSNIa > 0){
              CkPrintf("MethodFeedbackSTARSS SNe %d %d level = %d age = %f\n", nSNII, nSNIa, block->level(), age);
          }
          #endif

          if (nSNII+nSNIa > 0){
            /* set feedback properties based on number and types of SN */
            double energySN = (nSNII+nSNIa) * 1.0e51;

            SNMassEjected = SNII_ejecta_mass_Msun * nSNII +
                            SNIa_ejecta_mass_Msun * nSNIa; // AE: split this in two channels for yields

            const double starZ = pmetal[ipdmf] / z_solar;

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
        if (this->stellar_winds_){

          const double starZ = pmetal[ipdmf] / z_solar;

          determineWinds(age, &windEnergy, &windMass, &windMetals,
                         pmass_solar,
                         starZ, tunit, block->state()->dt());

          if (windMass > 0){
            #ifdef DEBUG_FEEDBACK_STARSS
              CkPrintf("STARSS_FB: Adding stellar winds...\n");
            #endif
            this->deposit_feedback( block, windEnergy, windMass, windMetals,
                                    pvx[ipdv],pvy[ipdv],pvz[ipdv],
                                    px[ipdp],py[ipdp],pz[ipdp],
                                    ix, iy, iz, 1, 0, 0, starZ); // removed P3

          } // if wind mass > 0
        } // if winds

        // windMass and SNMassEjected are in units of Msun
        pmass[ipdm] -= std::max(0.0,
                       (windMass + SNMassEjected) /
                       (munit/enzo_constants::mass_solar));


        // ionizing radiation
        if (this->radiation_) {
          double Psi_ion;
          if (age < 3.5) {
              Psi_ion = 500; // units of Lsun/Msun
          }
          if (age >= 3.5 && age <= 25){
              Psi_ion = 60. * pow(age/3.5, -3.6) + 470 * pow(age/3.5, 0.045-1.82*std::log(age));
          }
          double lum_unit = munit * lunit * lunit / (tunit*tunit*tunit);
          plum[ipdlum] = Psi_ion * pmass_solar * enzo_constants::luminosity_solar / lum_unit; // erg/s 
        } // if radiation 


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
                                              const double starZ)

 const throw(){
  /*
   This routine will create a cube of coupling particles, where we determine
      the feedback quantities.  The vertices of the cube are ~coupled particles
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
  double eunit = munit*vunit*vunit; // energy units (NOTE: not specific energy) 
  
  const EnzoConfig * enzo_config = enzo::config();

  bool AnalyticSNRShellMass = this->analytic_SNR_shell_mass_;

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

  int size = mx*my*mz;

  const int rank = cello::rank();

  double cell_volume_code = hx*hy*hz;
  double cell_volume = cell_volume_code*enzo_units->volume();

  // conversion from code_density to mass in Msun
  double rho_to_m = rhounit*cell_volume / enzo_constants::mass_solar;

  enzo_float * d           = (enzo_float *) field.values("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");

  enzo_float * vx          = (enzo_float *) field.values("velocity_x");
  enzo_float * vy          = (enzo_float *) field.values("velocity_y");
  enzo_float * vz          = (enzo_float *) field.values("velocity_z");

  enzo_float * mf          = (enzo_float *) field.values("metal_density");

  enzo_float * temperature = (enzo_float *) field.values("temperature");

  enzo_float * dHI    = field.is_field("HI_density") ? 
          (enzo_float*) field.values("HI_density") : NULL;
  enzo_float * dHII   = field.is_field("HII_density") ? 
          (enzo_float*) field.values("HII_density") : NULL;
  enzo_float * dHeI   = field.is_field("HeI_density") ? 
          (enzo_float*) field.values("HeI_density") : NULL;
  enzo_float * dHeII  = field.is_field("HeII_density") ? 
          (enzo_float*) field.values("HeII_density") : NULL;
  enzo_float * dHeIII = field.is_field("HeIII_density") ? 
          (enzo_float*) field.values("HeIII_density") : NULL;
  enzo_float * d_el   = field.is_field("e_density") ?
          (enzo_float*) field.values("e_density") : NULL;
 
  enzo_float * dH2I   = field.is_field("H2I_density") ? 
          (enzo_float*) field.values("H2I_density") : NULL;
  enzo_float * dH2II  = field.is_field("H2II_density") ? 
          (enzo_float*) field.values("H2II_density") : NULL;
  enzo_float * dHM    = field.is_field("HM_density") ? 
          (enzo_float*) field.values("HM_density") : NULL;

  enzo_float * dDI    = field.is_field("DI_density") ? 
         (enzo_float *) field.values("DI_density") : NULL;
  enzo_float * dDII   = field.is_field("DII_density") ? 
         (enzo_float *) field.values("DII_density") : NULL;
  enzo_float * dHDI   = field.is_field("HDI_density") ? 
         (enzo_float *) field.values("HDI_density") : NULL;

  // "deposit" fields that hold running total for all exploding particles in this block
  enzo_float *  d_dep_tot = (enzo_float *) field.values(i_d_dep);
  enzo_float * te_dep_tot = (enzo_float *) field.values(i_te_dep);
  enzo_float * ge_dep_tot = (enzo_float *) field.values(i_ge_dep);
  enzo_float * mf_dep_tot = (enzo_float *) field.values(i_mf_dep);
  enzo_float * vx_dep_tot = (enzo_float *) field.values(i_vx_dep);
  enzo_float * vy_dep_tot = (enzo_float *) field.values(i_vy_dep);
  enzo_float * vz_dep_tot = (enzo_float *) field.values(i_vz_dep);

  // holds just shell densities (used for refresh+accumulate)
  enzo_float * d_shell   = (enzo_float *) field.values(i_d_shell);

  // allocate another set of temporary deposit fields for this event 
  enzo_float *  d_dep = new enzo_float[size]; 
  enzo_float * te_dep = new enzo_float[size];
  enzo_float * ge_dep = new enzo_float[size];
  enzo_float * mf_dep = new enzo_float[size];
  enzo_float * vx_dep = new enzo_float[size];
  enzo_float * vy_dep = new enzo_float[size];
  enzo_float * vz_dep = new enzo_float[size];

  // initialize temporary deposit fields as zero
  // Not doing so gives in screwy results
  for (int i=0; i<size; i++) {
    d_dep [i] = 0;
    te_dep[i] = 0;
    ge_dep[i] = 0;
    mf_dep[i] = 0;
    vx_dep[i] = 0;
    vy_dep[i] = 0;
    vz_dep[i] = 0;
  }



  const int index = INDEX(ix,iy,iz,mx,my);

  int stretch_factor = 1.0; // put coupling particles one cell-width away from star particle
  const int nCouple = 26; // 3x3x3 cube minus central cell
  const double A = stretch_factor * hx;

  /*
      Make a cloud of coupling particles;
      each is on a grid with dx spacing
      from the host particles.
  */

  enzo_float CloudParticleVectorX[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0,  0, 0,  0,  0,  0,  0, 1, 1,  1, 1,  1,  1,  1,  1, 1};
  enzo_float CloudParticleVectorY[] = { 1,  1,  1,  0,  0, -1, -1, -1,  0, 1, 1,  1, 0,  0, -1, -1, -1, 1, 1,  1, 0,  0, -1, -1, -1, 0};
  enzo_float CloudParticleVectorZ[] = { 1,  0, -1,  1, -1,  1,  0, -1,  0, 1, 0, -1, 1, -1,  1,  0, -1, 1, 0, -1, 1, -1,  1,  0, -1, 0};

  enzo_float CloudParticlePositionX[nCouple];
  enzo_float CloudParticlePositionY[nCouple];
  enzo_float CloudParticlePositionZ[nCouple];

  for (int cpInd = 0; cpInd < nCouple; cpInd++)
  {
      double norm = sqrt(CloudParticleVectorX[cpInd] * CloudParticleVectorX[cpInd] +
                         CloudParticleVectorY[cpInd] * CloudParticleVectorY[cpInd] +
                         CloudParticleVectorZ[cpInd] * CloudParticleVectorZ[cpInd]);
      double inv_norm = 1.0 / norm;
      double xbaMag = A * A * norm * norm;

      // get position 
      CloudParticlePositionX[cpInd] = xp + CloudParticleVectorX[cpInd] * A;
      CloudParticlePositionY[cpInd] = yp + CloudParticleVectorY[cpInd] * A;
      CloudParticlePositionZ[cpInd] = zp + CloudParticleVectorZ[cpInd] * A;

      /* turn the vectors into unit-vectors */
      CloudParticleVectorZ[cpInd] *= inv_norm;
      CloudParticleVectorY[cpInd] *= inv_norm;
      CloudParticleVectorX[cpInd] *= inv_norm;
  }

     /* 
         transform to comoving with the star and transform velocities to momenta
         for easy momentum deposition, since the fluid variable that PPM pushes
         around is velocity
     */
  
  this->transformComovingWithStar(d,vx,vy,vz,up,vp,wp,mx,my,mz, 1);
  this->transformComovingWithStar(d_shell,vx_dep_tot,vy_dep_tot,vz_dep_tot,up,vp,wp,mx,my,mz, 1);

  const GrackleChemistryData * grackle_chem = enzo::grackle_chemistry();
  const int primordial_chemistry = (grackle_chem == nullptr) ?
    0 : grackle_chem->get<int>("primordial_chemistry");

  /* 
     Use averaged quantities across multiple cells so that deposition is stable.
     vmean is used to determine whether the supernova shell calculation should proceed:
     M_shell > 0 iff v_shell > v_gas 
  */ 
  double Z_mean=0, d_mean=0, n_mean=0, v_mean=0, mu_mean=0;
  const double dflt_mu = static_cast<double>(enzo::fluid_props()->mol_weight());

  for (int ix_ = ix-1; ix_ < ix+2; ix_++) {
    for (int iy_ = iy-1; iy_ < iy+2; iy_++) {
      for (int iz_ = iz-1; iz_ < iz+2; iz_++) {
        int ind = INDEX(ix_,iy_,iz_, mx,my);
        Z_mean += mf[ind] / d[ind];
        // TODO: make EnzoComputeMolecularWeight, and access mu_field here?

        double mu;
        if (primordial_chemistry > 0) {
          mu = d_el[ind] + dHI[ind] + dHII[ind] + 0.25*(dHeI[ind]+dHeII[ind]+dHeIII[ind]);

          if (primordial_chemistry > 1) {
            mu += dHM[ind] + 0.5*(dH2I[ind]+dH2II[ind]);
          }
          if (primordial_chemistry > 2) {
            mu += 0.5*(dDI[ind] + dDII[ind]) + dHDI[ind]/3.0;
          }
          mu /= d[ind];
        } else {
          mu = dflt_mu;
          // in an older version, mu = dflt_mu/d[ind], but I think that was a
          // typo
        }

        mu_mean += mu;
        d_mean += d[ind];
      } 
    }
  }
  
  v_mean *= vunit / 27; // cm/s
  Z_mean *= 1/(27 * enzo_constants::metallicity_solar); // Zsun
  mu_mean /= 27;
  d_mean *= rhounit/27; // g/cm^3
  n_mean = d_mean / (enzo_constants::mass_hydrogen*mu_mean);

  #ifdef DEBUG_FEEDBACK_STARSS
    CkPrintf("STARSS_FB: Zmean = %e Dmean = %e (%e) mu_mean = %e ", Z_mean, d_mean, d_mean / rhounit, mu_mean);
    CkPrintf("Nmean = %f vmean = %f\n", n_mean, v_mean/1e5);
  #endif

  double Z_Zsun = Z_mean;
  double fz = std::min(2.0, pow(std::max(Z_Zsun,0.01),-0.14));

  // Cooling radius as in Hopkins, but as an average over cells
  double CoolingRadius = 28.4 * pow(std::max(0.01, n_mean), -3.0 / 7.0) * pow(ejectaEnergy / 1.0e51, 2.0 / 7.0) * fz; // pc
  double coupledEnergy = ejectaEnergy;
  double cellwidth_pc = hx*lunit / enzo_constants::pc_cm; // pc
  double dxRatio = cellwidth_pc / CoolingRadius;

  /* 
     We want to couple one of four phases: free expansion, Sedov-taylor, shell formation, or terminal
     The first three phases are take forms from Kim & Ostriker 2015, the last from Cioffi 1988
  */

 
  // if we resolve free expansion, energy is mostly thermally coupled
  double p_free = 1.73e4*sqrt(ejectaMass * ejectaEnergy/1e51 / 3.0); 
  double r_free = 2.75*pow(ejectaMass / 3 / n_mean, 1.0/3); // free expansion radius eq. 2
  
  double t3_sedov = pow(std::max(r_free, std::max(cellwidth_pc, r_free)) * enzo_constants::pc_cm / 
                    (5.0 * enzo_constants::pc_cm * pow(ejectaEnergy / 1e51 / n_mean, 1.0 / 5.0)), 5./ 2.);
  
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
    // Thornton, 1998, sqrt(2 * M_R * E_R), eqs 15, 21, 26, 32
    pTerminal = 2.7298e5 * pow(ejectaEnergy/1e51, 13./14.) * pow(n_mean, -0.12) * pow(Z_Zsun, -0.14); // Msun km/s
  }

  else {
    pTerminal = 8.5721e5 * pow(ejectaEnergy/1e51, 13./14.) * pow(n_mean, -0.12);
  } 

  // compute the temperature (returns temperature in Kelvin)
  EnzoComputeTemperature compute_temperature(enzo::fluid_props(),
                                             enzo_config->physics_cosmology);

  compute_temperature.compute(enzo_block);

  double T = temperature[index];

  double gamma = 5.0/3.0;
  double cSound = sqrt(gamma * 
                       enzo_constants::kboltz*T/(mu_mean*enzo_constants::mass_hydrogen)) / 1e5; // km/s

  // fading radius of a supernova, using gas energy of the host cell and ideal gas approximations
  double r_fade = std::max(66.0*pow(ejectaEnergy/1e51, 0.32)*pow(n_mean, -0.37)*pow(cSound/10, -2.0/5.0), 
                           CoolingRadius * 1.5);
  double fadeRatio = cellwidth_pc/r_fade;

  #ifdef DEBUG_FEEDBACK_STARSS
    CkPrintf("STARSS_FB: Fading: T = %e; Cs = %e; R_f = %e; fadeR = %f\n", T, cSound, r_fade, fadeRatio);
  #endif

  double coupledMomenta = 0.0;
  double eKinetic = 0.0;

  double cw_eff = cellwidth_pc;

  double dx_eff = cw_eff / CoolingRadius; // our resolution with respect to the cooling radius
  double fader = cw_eff / r_fade;

  double shellVelocity = 413.0 * pow(n_mean, 1.0 / 7.0) * pow(Z_Zsun, 3.0 / 14.0) 
                               * pow(coupledEnergy / 1e51, 1.0 / 14.0);
                               
  double ratio_pds = cw_eff/r_shellform;

  shellVelocity *=  ratio_pds > 1 ? pow(dx_eff, -7.0 / 3.0) : 1; //km/s

  float beta =  std::min(70.0, std::max( 1.0, shellVelocity / std::max(1.0, cSound)));

  // critical density to skip snowplough (remnant combines with ISM before radiative phase); eq 4.9 cioffi 1988
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
       cellwidth_pc, r_free, r_shellform, CoolingRadius, r_fade, t3_sedov, r_merge);
  #endif

  if (! winds)
    {   // this calculation for SNe only

        // free-expansion phase
        if (cw_eff < r_free){
          // Most of the energy should be deposited as thermal here. Using
          // Sedov momentum multiplied by dx/R_free for a smooth transition between the two phases.
          // In the limit of dx = R_free, coupledMomenta(free expansion) = coupledMomenta(sedov).
          coupledMomenta = std::min(p_sedov, pTerminal * dx_eff) * cw_eff/r_free;
          #ifdef DEBUG_FEEDBACK_STARSS
            CkPrintf("STARSS_FB: Coupling Free Expansion phase: p = %e\n", coupledMomenta);
          #endif
        }
 
        // Sedov phase
        if (r_free < cw_eff && dx_eff <= 1){
          coupledMomenta = std::min(p_sedov, pTerminal * dx_eff);
          #ifdef DEBUG_FEEDBACK_STARSS       
            CkPrintf("STARSS_FB: Coupling Sedov-Terminal phase: p = %e (ps = %e, pt = %e, dxe = %e)\n", 
                      coupledMomenta, p_sedov, pTerminal, dx_eff);
          #endif
        }

        // terminal phase
        if (dx_eff > 1){
          // Scaling pTerminal by tanh(dx/R_cool) here to account for error in determining where exactly
          // the terminal phase begins. If we're straddling the edge between Sedov and terminal phases
          // (i.e. when dx_eff = R_cool approximately), there could be some extra pressure-driven
          // expansion even though the algorithm detects that we're in the terminal phase.
          // tanh(dx/R_cool) -> 1 as dx increases
          coupledMomenta = pTerminal * tanh(dx_eff);
          #ifdef DEBUG_FEEDBACK_STARSS
            CkPrintf("STARSS_FB: Coupling Terminal phase: p = %e; dx_eff = %e\n", coupledMomenta, dx_eff);
          #endif
        }

        // fading phase
        if (faded && this->fade_SNR_){
          coupledMomenta *= (1.0 - tanh(pow(1.25*fader, 2)));
          #ifdef DEBUG_FEEDBACK_STARSS 
            CkPrintf("STARSS_FB: Coupling Fading phase: p = %e\n", coupledMomenta);
          #endif
        }

        #ifdef DEBUG_FEEDBACK_STARSS
        CkPrintf("STARSS_FB: Checking critical density metric..." 
                 "(n_mean = %e; N_Crit = %e; factors: %e %e %e; beta = %e/%e == %e; rmerge = %e)\n", 
                  n_mean, nCritical, pow(n_mean * T / 1e4, 7.0/9.0), 
                  pow(ejectaEnergy/1e51, 1.0/9.0), pow(fz, 1.0/3.0), 
                  shellVelocity , cSound, beta, r_merge);
        #endif
           
        if (T > 1e6 && coupledMomenta > 1e5) {
          #ifdef DEBUG_FEEDBACK_STARSS
            CkPrintf("STARSS_FB: Coupling high momenta to very hot gas!! (p= %e, T= %e, n_c = %e)\n", coupledMomenta, T, nCritical);
          #endif
        }
    } // endif no winds

    else { // if winds
      coupledMomenta = sqrt(ejectaMass*enzo_constants::mass_solar * 0.5 * ejectaEnergy)/enzo_constants::mass_solar/1e5;
    }

    #ifdef DEBUG_FEEDBACK_STARSS
      CkPrintf("STARSS_FB: Calculated p = %e (sq_fact = %e; p_f = %e; p_t = %e; mcell = %e; mcpl = %e)\n",
         coupledMomenta, (d_mean * cell_volume/enzo_constants::mass_solar) / ejectaMass * 63, p_free, 
         pTerminal, d_mean * cell_volume /enzo_constants::mass_solar, ejectaMass/27.0);  
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
    if (dx_eff < 1) // free expansion
       shellMass = 4*cello::pi/3 * d_mean * pow(cw_eff*enzo_constants::pc_cm,3) / enzo_constants::mass_solar; 

    else if (dx_eff > 1) {
      if (Z_Zsun > 0.01)
        shellMass = 1.41e4*pow(ejectaEnergy/1e51, 6./7.) * pow(d_mean, -0.24) * pow(Z_Zsun, -0.27);
      else
        shellMass = 4.89e4*pow(ejectaEnergy/1e51, 6./7.) * pow(d_mean, -0.24);
    }
   
    if (shellMass > 0)
    {
      // get total mass in cells affected by SN       
      for (int ix_ = ix-1; ix_ <= ix+1; ix_++) {
        for (int iy_ = iy-1; iy_ <= iy+1; iy_++) {
          for (int iz_ = iz-1; iz_ <= iz+1; iz_++) {
            int flat = INDEX(ix_,iy_,iz_,mx,my);

            // cell left edges for CiC (starting from ghost zones)
            double xcell = xm + (ix_+0.5 - gx)*hx; 
            double ycell = ym + (iy_+0.5 - gy)*hy;
            double zcell = zm + (iz_+0.5 - gz)*hz;

            double window = Window(xp-xcell, yp-ycell, zp-zcell, hx); // CiC fraction

            if (window <= 0.0) continue;

            centralMass += window * d[flat];
            centralMetals += window * mf[flat];
          }
        }
      } // endfor ix_
      
      centralMass *= rho_to_m; // Mass in Msun
      centralMetals *= rho_to_m;

      // Can't let host cells evacuate completely! 
      // Enforce maximum shellMass 
      shellMass = std::min(shellMass, maxEvacFraction*centralMass);

     } // endif shellMass > 0
      
  } // endif coupledEnergy >0 && AnalyticSNRShellMass && !winds
  
  double coupledMass = shellMass + ejectaMass;
  eKinetic = coupledMomenta * coupledMomenta / (2.0 *(coupledMass)) * enzo_constants::mass_solar * 1e10;
  
  double shellMetals = std::min(maxEvacFraction*centralMetals, 
                                Z_Zsun * enzo_constants::metallicity_solar * shellMass); //Msun

#ifdef DEBUG_FEEDBACK_STARSS
  if (AnalyticSNRShellMass) {
    CkPrintf("STARSS_FB: Shell_m = %e Shell_z = %e shell_V= %e P = %e M_C = %e\n",
             shellMass, shellMetals, shellVelocity, coupledMomenta, centralMass);    
  }  
#endif
 
  #ifdef DEBUG_FEEDBACK_STARSS
    CkPrintf("STARSS_FB: Ekinetic = %e Mass = %e\n",
             eKinetic, d_mean * pow(lunit * hx, 3) / enzo_constants::mass_solar);
  #endif

    if (eKinetic > 1e60)
    {
      #ifdef DEBUG_FEEDBACK_STARSS
        CkPrintf("STARSS_FB: winds = %d, Ekinetic = %e Mass = %e\n",
                  winds, eKinetic, d_mean * pow(lunit * hx, 3) / enzo_constants::mass_solar);
      #endif
      ERROR("EnzoMethodFeedbackSTARSS::deposit_feedback()","Ekinetic > 1e60 erg!\n");
    }

    double coupledGasEnergy = std::max(ejectaEnergy - eKinetic, 0.0);
    #ifdef DEBUG_FEEDBACK_STARSS
      CkPrintf("STARSS_FB: Coupled Gas Energy = %e\n", coupledGasEnergy);
    #endif

    if (dxRatio > 1.0 && !winds){
      // Fade internal energy deposited if our resolution is worse
      // than a cooling radius.
      //
      // Don't need to do this for stellar winds because winds are 
      // a lot weaker than SNe, so it makes essentially no difference.
      coupledGasEnergy = (coupledGasEnergy * pow(dxRatio, -6.5));

      #ifdef DEBUG_FEEDBACK_STARSS
        CkPrintf("STARSS_FB: Reducing gas energy... GE = %e\n", coupledGasEnergy);
      #endif
    }

  double coupledMetals = 0.0;//, SNIAmetals = 0.0, SNIImetals = 0.0, P3metals = 0.0;
  if (winds) coupledMetals = ejectaMetals; // winds only couple to metals
  if (AnalyticSNRShellMass && !winds) coupledMetals += shellMetals;
  
  ejectaMetals += nSNIA * 1.4; // add ejecta from Type Ia SNe 
  ejectaMetals += nSNII * (1.91 + 0.0479 * std::max(starZ, 1.65)); // add ejecta from Type II

  coupledMetals += ejectaMetals;

  #ifdef DEBUG_FEEDBACK_STARSS
    CkPrintf("STARSS_FB: Coupled Metals: %e %e %e\n", ejectaMetals, shellMetals, coupledMetals);
  #endif

  if (!winds) coupledEnergy = std::min(coupledEnergy, eKinetic);

  // Subtract shell mass from the central cells
  double minusRho=0, minusZ=0, msubtracted=0;
  double remainMass=shellMass/rho_to_m, remainZ = shellMetals/rho_to_m;
  if (shellMass > 0 && AnalyticSNRShellMass)
  {
    double zsubtracted=0;
    msubtracted=0;

    for (int ix_ = ix-1; ix_ <= ix+1; ix_++) {
      for (int iy_ = iy-1; iy_ <= iy+1; iy_++) {
        for (int iz_ = iz-1; iz_ <= iz+1; iz_++) {
          int flat = INDEX(ix_,iy_,iz_,mx,my);

          // cell left edges for CiC (starting from ghost zones)
          double xcell = xm + (ix_+0.5 - gx)*hx; 
          double ycell = ym + (iy_+0.5 - gy)*hy;
          double zcell = zm + (iz_+0.5 - gz)*hz;

          double window = Window(xp-xcell, yp-ycell, zp-zcell, hx); // CiC fraction
          if (window <= 0) continue;

          double dpre = d[flat];
          double zpre = mf[flat];
          double pre_z_frac = zpre / dpre;

          // subtract values from the "deposit" fields
          d_dep[flat] -= std::min(window * remainMass, maxEvacFraction*dpre);

          minusRho    += -1*d_dep[flat];
          msubtracted += -1*d_dep[flat];
         
          mf_dep[flat] -= std::min(window * remainZ, maxEvacFraction*zpre); 

          minusZ      += -1*mf_dep[flat];
          zsubtracted += -1*mf_dep[flat];

        } // endfor iz_
      } // endfor iy_
    } // endfor ix_
    remainMass -= msubtracted;
    remainZ    -= zsubtracted;
  } // endif shellMass > 0 && AnalyticSNRShellMass

  minusRho *= rho_to_m; //Msun
  minusZ   *= rho_to_m; //Msun

  double minusM = minusRho;
  if (minusM != coupledMass - ejectaMass && shellMass > 0 && AnalyticSNRShellMass)
  {
    #ifdef DEBUG_FEEDBACK_STARSS
      CkPrintf("STARSS_FB: Of %e, only subtracted %e; rescaling the coupling mass\n", shellMass,
                minusRho);
    #endif
  
    coupledMass = minusM + ejectaMass;
    coupledMetals = minusZ + ejectaMetals; // metal_density in Msun/code_volume

    // if we're in here, we found an inconsistent place where the mass within 
    // cannot support the expected shell.
    // the only real choice, to keep things from exploding in a bad way, 
    // is to rescale the momentum accordingly.
    // If this is triggered, dont expect the terminal momenta relations to hold up.

    if (coupledMomenta*coupledMomenta / 
        (2.0*coupledMass*enzo_constants::mass_solar)*enzo_constants::mass_solar * 1e10 > 1e51) {
      coupledMomenta = std::min(coupledMomenta, sqrt(2.0 * (coupledMass*enzo_constants::mass_solar)
                              * ejectaEnergy)/enzo_constants::mass_solar/1e5);

      #ifdef DEBUG_FEEDBACK_STARSS
        CkPrintf("STARSS_FB: rescaled momentum to %e (est KE = %e)\n", 
                  coupledMomenta, coupledMomenta*coupledMomenta / (2*coupledMass) 
                * enzo_constants::mass_solar * 1e10);
      #endif
    }
  
  } // endif minusM != coupledMass - ejectaMass

  coupledEnergy += coupledGasEnergy;

  #ifdef DEBUG_FEEDBACK_STARSS
      CkPrintf("STARSS_FB: Depositing -- Total Energy = %1.4e; Gas Energy = %1.4e;"
               " Mass = %1.4e; Metal Mass = %1.4e; Momentum = %1.4e\n",
               coupledEnergy, coupledGasEnergy, coupledMass, coupledMetals, coupledMomenta);
  #endif

  // put everything back into code units before CIC
  // these values correspond to TOTAL (energy, mass, momentum)
  // over all the coupling particles.

  coupledEnergy /= (eunit * cell_volume_code); // energy -> energy density
  coupledGasEnergy /= (eunit * cell_volume_code);
  coupledMass /= rho_to_m; // mass -> density
  coupledMetals /= rho_to_m;
  coupledMomenta /= (rho_to_m * vunit / 1e5); // momentum -> momentum_density 


  // Create arrays to pass into CiC routine
  enzo_float coupledMass_list[nCouple], coupledMetals_list[nCouple];
  enzo_float coupledMomenta_x[nCouple], coupledMomenta_y[nCouple], coupledMomenta_z[nCouple];
  enzo_float coupledEnergy_list[nCouple], coupledGasEnergy_list[nCouple];
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
  enzo_float left_edge[3] = {xm-gx*hx, ym-gy*hy, zm-gz*hz};

  // CiC deposit mass/energy/momentum
  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMass_list, d_dep, &left_edge,
   &mx, &my, &mz, &hx, &A);

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMomenta_x, vx_dep, &left_edge,
   &mx, &my, &mz, &hx, &hx);

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMomenta_y, vy_dep, &left_edge,
   &mx, &my, &mz, &hx, &hx);

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMomenta_z, vz_dep, &left_edge,
   &mx, &my, &mz, &hx, &hx);

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

  FORTRAN_NAME(cic_deposit)
  (&CloudParticlePositionX, &CloudParticlePositionY,
   &CloudParticlePositionZ, &rank, &nCouple, &coupledMass_list, d_shell, &left_edge,
   &mx, &my, &mz, &hx, &A);

  // copy deposited quantites to original fields
  for (int i=0; i<size; i++){
    double d_old = d[i]; 
    d[i] += d_dep[i];
    double d_new = d[i];

    double cell_mass = d_new*cell_volume_code;
    double M_scale = d_new / d_old;

    mf[i] += mf_dep[i]; 

    // need to rescale specific energies to account for added mass
    te[i] = te[i]/M_scale + te_dep[i]/d_new;
    ge[i] = ge[i]/M_scale + ge_dep[i]/d_new;

    vx[i] += vx_dep[i];
    vy[i] += vy_dep[i];
    vz[i] += vz_dep[i];
    // Rescale color fields to account for new densities.
    // Don't need to rescale metal_density because we already deposited
    // into the metal_density field.

    EnzoMethodStarMaker::rescale_densities(enzo_block, i, d_new/d_old);

    // undo rescaling of metal_density
    mf[i] /= (d_new/d_old);

    // add deposited quantities to fields that track depositions
    // of all star particles in the block this cycle
     d_dep_tot[i] += d_dep[i];
    mf_dep_tot[i] += mf_dep[i];
    te_dep_tot[i] += te_dep[i];
    ge_dep_tot[i] += ge_dep[i];
    vx_dep_tot[i] += vx_dep[i];
    vy_dep_tot[i] += vy_dep[i];
    vz_dep_tot[i] += vz_dep[i];
  }
  // transform velocities back to "lab" frame
  // convert velocity (actually momentum density at the moment) field back to velocity 
  this->transformComovingWithStar(d,vx,vy,vz,up,vp,wp,mx,my,mz, -1);
  this->transformComovingWithStar(d_shell,vx_dep_tot,vy_dep_tot,vz_dep_tot,up,vp,wp,mx,my,mz, -1);

  // clear temporary fields 
  delete []  d_dep;
  delete [] mf_dep;
  delete [] te_dep;
  delete [] ge_dep;
  delete [] vx_dep;
  delete [] vy_dep;
  delete [] vz_dep;
}



double EnzoMethodFeedbackSTARSS::timestep (Block * block) throw()
{
  // In general this is not needed, but could imagine putting timestep
  // limiters in situations where, for example, one would want
  // dt < star_lifetime (or something like that), especially if
  // important things happen throughout the star's lifetime.
  EnzoUnits * enzo_units = enzo::units();

  EnzoMethodStarMaker* starmaker_method =
    (EnzoMethodStarMaker*)(enzo::problem()->method("star_maker"));
  ASSERT("EnzoMethodFeedbackSTARSS::timestep",
         "requires \"star_maker\" method", starmaker_method != nullptr);

  double dtStar = std::numeric_limits<double>::max();
  if (block->level() >= starmaker_method->sf_minimum_level()){
    const auto dt   = block->state()->dt();
    const auto time = block->state()->time();
    const double pSNmax = 0.0005408 * starmaker_method->minimum_star_mass() *
      dt * time / enzo_constants::Myr_s * 1.25;
    if (pSNmax > 1.0) dtStar = dt * 1.0 / pSNmax;
  }

  return dtStar;
}


// #endif
