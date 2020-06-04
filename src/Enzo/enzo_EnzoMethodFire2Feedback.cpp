
/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodFire2Feedack.cpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
/// @date
/// @brief  Implements the FIRE2 model for stellar feedback
///         as implemented into Enzo as the MechStars model
///         by Azton Wells. This is intended to be a direct copy-over of the
///         MechStars methods.
///         TODO: List differences / changes here
///


#include "cello.hpp"
#include "enzo.hpp"

#include <time.h>

#ifdef NOTDEFINED // for now... since not done coding

// =============================================================================
// splice these off to a different file (later)

int determineSN(float age, int* nSNII, int* nSNIA,
                float massmass_solar, float TimeUnits, float dt){

    const EnzoConfig * enzo_config = enzo::config();

    if (NEvents > 0){
        *nSNII = 1;
        NEvents -= 1;
        return SUCCESS;
    }
    /* else, calculate SN rate, probability and determine number of events */
    int seed = time(NULL):
    *nSNII = 0;
    *nSNIA = 0;
    float RII=0, RIA=0, PII=0, PIA=0, random = 0;
    if (enzo_config->method_feedback_fire2_single_sn && NEvents < 0)
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
            srand(seed);
        // printf("Zcpl = %e", zCouple);
            PII = RII * massmass_solar / cello::Myr_s *TimeUnits*dt;
            // printf("PII =%f\n %f %e %f\n", PII, RII, massmass_solar, age);
            random = float(rand())/float(RAND_MAX);
            if (PII > 1.0 && enzo_config->method_feedback_fire2_unrestricted_SN){
                int round = (int)PII;
                *nSNII = round;
                PII -= round;
            }
            if (PII > 1.0 && !enzo_config->method_feedback_fire2_unrestricted_SN){
                ERROR("determineSN: PII too large!");
            }
            int psn = *nSNII;
            if (random < PII){
                *nSNII = psn+1;
            }
        }
        // printf("RANDOM = %f\n", random);
        // printf("N SNII=%d\n",*nSNII);

        if (RIA > 0){
            srand(seed);
            PIA = RIA*massmass_solar / cello:Myr_s *TimeUnits*dt;
            float random = float(rand())/float(RAND_MAX);

            if (PIA > 1.0 && enzo_config->method_feedback_fire2_unrestricted_SN)
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
    }
        return 1;
}

int determineWinds(float age, float* eWinds, float* mWinds, float* zWinds,
                   float massMsun, float zZsun, float TimeUnits, float dt){
    /* AJE: Disclaimer!!!
       The version of this function implemented by Azton matches the rates given
       in the FIRE2 method paper but does not seem to match up with the rates
       given in the actual FIRE code... below is the latter, but this might be
       an issue to use if its non-pubic (unless this is a mistake in the publication?)
    */


    // age in Myr

    const EnzoConfig * enzo_config = enzo::config();

    float windE = 0,  windM = 0, windZ = 0.0;
    float wind_factor = 0.0;
    float e_factor = 0.0;


    if (enzo_config->method_feedback_fire2_stellarwinds){

      if(zZsun > 3.0) zZsun = 3.0;
      if(zZsun < 0.01) zZsun = 0.01;

      if (age <= 1.0){
        wind_factor = 11.6846;
      } else {
        if (age <= 3.5){
          wind_factor = 11.6846 * zZsun *
                        pow(10.0, 1.838*(0.79+log10(zZsun))*(log10(age)));
        } else {
          if (age <= 100.0){
            wind_factor = 72.1215*pow(0.001*age/0.0035,-3.25)+0.0103;
          } else {
            wind_factor = 1.03*pow(0.001*age,-1.1)/(12.9-log(age*0.001));
          } // age < 100
        } // age < 3.5

      } // age < 1.0

      if (age <= 100.0){
        e_factor = 0.0013+16.0/(1.0+pow(age*0.001/0.0025,1.4)+pow(age*0.001/0.01,5.0));
      } else{
        e_factor = 0.0013;
      }

      wind_factor *= (enzo_config->method_feedback_fire2_gasreturnfraction *\
                     0.001 * (dt * (TimeUnits / cello::Myr_s)); // 0.001 converts from Msun/Gyr to Msun/Myr
      wind_factor = 1.0 - exp(-wind_factor);
      wind_factor = std::min(1.0, wind_factor);
      wind_factor *= 1.4 * 0.291175;

      double n_wind = (double)floor(wind_factor/rfrac); wind_factor -= n_wind *rfrac;
      windM = n_wind * rfrac;
      if( float(rand())/float(RAND_MAX) < wind_factor / rfrac){windM += rfrac};

      windM *= massMsun; // convert to total mas from fraction
      windZ = (0.016+0.0041+0.0118) * windM; // C+N+O yields with no Z scaling (total masss of metals)
      windE = e_factor * 1e12 * windM
    }


    *mWinds = windM;
    *zWinds = windZ;
    *eWinds = windE;

    return 1;
/*

        if (enzo_config->method_feedback_fire2_stellarwinds &&
            oldEnough && massMsun > 11){

            if (0.001 < age && age < 1.0){
                wind_factor =4.763 * min((0.01 + zZsun), 1.0) ;
            }
            if (1 <= age && age < 3.5){
                wind_factor = 4.763*min(0.01+zZsun, 1.0)*
                    pow(age, 1.45+0.08*min(log(zZsun), 1.0));
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
            windM = massMsun * wind_factor; //Msun/Gyr
            windM = windM*dtFixed*TimeUnits/3.1557e16; //Msun
            // if (debug) printf("First winds mass = %e\nFrom wf = %f, dt=%f Z = %e\n", windM, wind_factor, dtFixed, zZsun);
            //printf("eFactor = %f age = %f\n", e_factor, age);
            if (windM > massMsun){
                printf("Winds too large Mw = %e, Mp = %e age=%f, Z = %e\n",
                    windM, massMsun, age, zZsun);
                windM = 0.125*massMsun; // limit loss to huge if necessary.
            }
            windZ = max(0.02, 0.016+0.0041*max(zZsun, 1.65)+0.0118)*windM;
            windE = e_factor * 1e12 * windM;
            //fprintf(stdout, "Age = %e Ewinds = %e Mwinds = %e Zwinds = %e  Zsun = %e\n",
              //          age, windE, windM, windZ, zZsun);
            *mWinds = windM;
            *zWinds = windZ;
            *eWinds = windE;
        }

        return SUCCESS;
    }
*/
}

// =============================================================================
// =============================================================================
// =============================================================================

EnzoMethodFire2Feedack::EnzoMethodFire2Feedack
()
  : Method()
{
  FieldDescr * field_descr = cello::field_descr();
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  // required fields
  this->required_fields_ = std::vector<std::string>
                           {"density","pressure","internal_energy",
                            "total_energy"};

  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  sf_minimum_level_ = enzo_config->method_feedback_fire2_minimum_level;
  single_sn_        = enzo_config->method_feedback_fire2_single_sn;

  return;
}

void EnzoMethodFire2Feedack::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | sf_minimum_level_;
  p | single_sn_;

  return;
}

void EnzoMethodFire2Feedack::compute (Block * block) throw()
{

  if (block->is_leaf()){
    this->compute_(block);
  }

  block->compute_done();

  return;
}

void EnzoMethodFire2Feedback::compute_ (Block * block)
{

  //----------------------------------------------------
  // some constants here that might get moved to parameters or something else
  const float SNII_ejecta_mass = 10.5;
  const float SNIa_ejecta_mass = 1.4;
  const float z_solar          = 0.02; // Solar metal fraction (fire2)
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

  // apply feedback depending on particle type
  // for now, just do this for all star particles

  int it = particle.type_index("star");
  int count = 0;

  // AE: Copying this from MechStars_FeedbackRoutine

  int numSN = 0; // counter of SN events
  int count = 0; // counter of particles

  if (particle.num_particles(it) <= 0) return; // nothing to do here

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
    enzo_float *plifetime=0, *pcreation=0, *pmass=0, *pmetal=0;

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

    psncounter = (enzo_int *) particle.attribute_array(it, ia_sn, ib);

    int np = particle.num_particles(it,ib);

    for (int ip=0; ip<np; ip++){
      // AE: Check and see if these actually differ.... (ask James?)
      int ipdp = ip*dp; // pos
      int ipdm = ip*dm; // mass
      int ipdv = ip*dv; // velocity
      int ipdl = ip*dl; // lifetime
      int ipdc = ip*dc; // creation time
      int ipdmf = ip*dmf; // metallicity
      int ipsn  = ip*dsn; // number of SNe coutner

      if (pmass[ipdm] > 0.0 && plifetime[ipdl] > 0.0){
        const double age = (current_time - pcreation[ipdc]) * enzo_units->time() / cello::Myr_s;
        count++; // increment particles examined here

        // compute coordinates of central feedback cell
        // this must account for ghost zones
        double xcell = (xpos - xm) / hx + gx - 0.5;
        double ycell = (ypos - ym) / hy + gy - 0.5;
        double zcell = (zpos - zm) / hz + gz - 0.5;

        int ix       = ((int) floor(xcell + 0.5));
        int iy       = ((int) floor(ycell + 0.5));
        int iz       = ((int) floor(zcell + 0.5));

        // do cell shifting here if we need to... assume not for now ...

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

          determineSN(age, &nSNII, &nSNIa, pmass[ipdm] * enzo_units->mass_units() / cello::mass_solar,
                      enzo_units->time_units(), block->dt());

          numSN += nSNII + nSNIa;

#ifdef FEEDBACK_DEBUG
          if (numSN){
            CkPrintf("Fire2Feedback SNe %d %d level = %d age = %f\n", nSNII, nSNIa, level, age);
          }
#endif

          /* AJE: Can I just change this to a single call to deposit feedback for
                  the total ejecta of SNe and winds ? */

          if (numSN > 0){
            /* set feedback properties based on number and types of SN */
            double energySN = (numSN) * 1.0e51;

            SNMassEjected = SNII_ejecta_mass * nSNII +
                            SNIa_ejecta_mass * nSNIa; // AE: split this in two channels for yields

            const double starZ = pmetal[ipdmf] / z_solar;

            // compute metal mass ejected here

            SNMetalEjected = ;

            /* Fixed mass ejecta */

            this->deposit_feedback( energySN, SNMassEjected, SNMetalEjected,
                                    pvx[ipdv],pvy[ipdv],pvz[ipdv],
                                    px[ipdp],py[ipdp],pz[ipdp],
                                    ix, iy, iz, mx, my, mz, 0, nSNII, nSNIa, starZ); // removed P3

            // add counter for number of SNe for the particle
            psncounter[ipsn] += numSN;



          } // if nSNII or nSNIa > 0
        } // if single SN

        // Now do stellar winds
        if (enzo_config->method_feedback_fire2_stellarwinds){

          const double starZ = pmetal[ipdmf] / z_solar;

          determineWinds(age, &windEnergy, &windMass, &windMetals,
                         pmass[ipdm] * pmass[ipdm] * enzo_units->mass_units() / cello::mass_solar,
                         starZ, enzo_units->time_units(), block->dt());

          /* Error messages go here */

          if (windMass > 0){
            this->deposit_feedback( windEnergy, windMass, windMetals,
                                    pvx[ipdv],pvy[ipdv],pvz[ipdv],
                                    px[ipdp],py[ipdp],pz[ipdp],
                                    ix, iy, iz, mx, my, mz, 1, 0, 0, 0.0); // removed P3

          } // if wind mass > 0
        } // if winds

        pmass[ipdm] -= std::max(0.0,
                       (windMass + SNMassEjected) *enzo_units->mass_units()/cello::mass_solar);

        //


      } // if mass and lifetime > 0

    } // end particle loop
  } // end batch loop


  if (count > 0){
    CkPrintf("Fire2Feedback: Num FB particles = %d  Events = %d  FeedbackTime %e\n",
              count, numSN, 0.00);
  }

  return;
}




// ----------------------------------------------------------------------------


void EnzoMethodFire2Feedack::deposit_feedbak (EnzoBlock * enzo_block,
                                              const float ejectaEnergy,
                                              const float ejectaMass,
                                              const float ejectaMetals,
                                              const enzo_float up, const enzo_float vp, const enzo_float wp,
                                              const enzo_float xp, const enzo_float yp, const enzo_float zp,
                                              const int ix, const int iy, const int iz,
                                              const int wind_mode, const int nSNII,
                                              const int nSNIa,
                                              enzo_float starZ);

) const throw(){
  /*
   This routine will create an isocahedron of coupling particles, where we determine
      the feedback quantities.  The vertices of the isocahedron are ~coupled particles
      and all have radius dx from the source particle.
      Each vertex particle will then be CIC deposited to the grid!
  */

  Field field = enzo_block->data()->field();
  // Obtain grid sizes and ghost sizes
  int mx, my, mz, gx, gy, gz, nx, ny, nz;
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  field.size(&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);
  enzo_block->data()->lower(&xm,&ym,&zm);
  enzo_block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;

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

  const float phi  = (1.0 + sqrt(5.0)) * 0.5 ; // golden ratio
  const float iphi = 1.0 / phi;

  /* A DODECAHEDRON+ISOCAHEDRON */
  const int nCouple = 26;
  const float A = stretchFactor * dx;

  /*
      Make a cloud of coupling particles;
      each is on a grid with dx spacing
      from the host particles.
  */

  enzo_float CloudParticleVectorX[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  // {1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0,
  // iphi, iphi, -iphi, -iphi, phi, phi, -phi, -phi,
  // 0, 0, 0, 0, 1, 1, -1, -1, phi,-phi, phi,-phi};
  enzo_float CloudParticleVectorY[] = {1, 1, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, 0, -1, -1, -1, 1, 1, 1, 0, 0, -1, -1, -1, 0};
  // {1,1,-1,-1, 1, 1, -1, -1, iphi, iphi, -iphi, -iphi,
  // phi, -phi, phi,-phi, 0, 0, 0, 0,1, 1, -1, -1,
  // phi, -phi, -phi, phi, 0, 0, 0, 0};
  enzo_float CloudParticleVectorZ[] = {1, 0, -1, 1, -1, 1, 0, -1, 0, 1, 0, -1, 1, -1, 1, 0, -1, 1, 0, -1, 1, -1, 1, 0, -1, 0};
  // {1,-1, 1,-1, 1,-1, 1,-1, phi,-phi, phi,-phi,
  // 0, 0, 0, 0, iphi, -iphi, iphi, -iphi,
  // phi, -phi, -phi, phi, 0, 0, 0, 0, 1, 1, -1, -1};
  enzo_float weightsVector[nCouple];
  /* Set position of feedback cloud particles */

  enzo_float CloudParticlePositionX[nCouple];
  enzo_float CloudParticlePositionY[nCouple];
  enzo_float CloudParticlePositionZ[nCouple];

  /*all possible values of x,y,z with origin particle at x=y=z=0.0 */
  for (int cpInd = 0; cpInd < nCouple; cpInd++)
  {
      float norm = sqrt(CloudParticleVectorX[cpInd] * CloudParticleVectorX[cpInd] +
                        CloudParticleVectorY[cpInd] * CloudParticleVectorY[cpInd] +
                        CloudParticleVectorZ[cpInd] * CloudParticleVectorZ[cpInd]);
      float inv_norm = 1.0 / norm;
      float xbaMag = A * A * norm * norm;
      /* in this cloud, take the coupling particle position as 0.5, 0.5, 0.5 */
      CloudParticlePositionX[cpInd] = xp - CloudParticleVectorX[cpInd]*inv_norm * A;
      CloudParticlePositionY[cpInd] = yp - CloudParticleVectorY[cpInd]*inv_norm * A;
      CloudParticlePositionZ[cpInd] = zp - CloudParticleVectorZ[cpInd]*inv_norm * A;
      weightsVector[cpInd] = 0.5 * (1. - 1. / (1. + 1. / 4. / cello::pi / 26. / xbaMag / cello::pi));
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
          ERROR("EnzoMethodFire2Feedack: deposit_feedback: NaN weight Vector!")
      }
  }

/* AJE LEFT OFF HERE */


  return;
}



double EnzoMethodFire2Feedack::timestep (Block * block) const throw()
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
    if (pSNmax > 1.0) dtStar = dt * 1.0 / pSNmax;
  }

  return dtStar;
}


#endif
