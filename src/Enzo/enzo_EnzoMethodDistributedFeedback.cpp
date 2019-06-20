/// See LICENSE_CELLO file for license and copyright information

///
///
///
///

#include "cello.hpp"
#include "enzo.hpp"


double s49Lookup(double m){
  // Returns a broken power-law approximation to the ionizing luminosity emitted
  // by a main sequence star of mass m, normalized to 10^49 photon/s
  // See Parravano et al. 2003 ApJ 584 797 (doi:10.1086/345807)
  if (m<5)
    {
      return 0;
    }
  else if (m<7)
    {
      return (2.23e-15*std::pow(m, 11.5));
    }
  else if (m<12)
    {
      return(3.69e-13*std::pow(m, 8.87));
    }
  else if (m<20)
    {
      return(4.8e-12*std::pow(m, 7.85));
    }
  else if (m<30)
    {
      return(3.12e-8*std::pow(m, 4.91));
    }
  else if (m<40)
    {
      return(2.8e-5*std::pow(m, 2.91));
    }
  else if (m<60)
    {
      return(3.49e-4*std::pow(m, 2.23));
    }
  else
    {
      return(2.39e-3*std::pow(m, 1.76));
    }
}


double s99_wind_mass(const double & t){
  // returns s99 wind mass in solar masses / yr per 1.0E6 Msun cluster
  // as a function of delay time (t) in units of 10^7 years
  //
  double m = 0.0;

  const double tp2 = t*t;
  const double tp4 = tp2*tp2;

  if (t <= 0.29) {
    m = m- 2.72379457924;
    m = m + 9.03549102928*t;
    m = m - 349.073894935*tp2;
    m = m + 6133.48804337*tp2*t;
    m = m - 45526.5891824*tp4;
    m = m + 160159.422053*tp4*t;
    m = m - 254942.557778*tp4*tp2;
    m = m + 133227.581992*tp4*tp2*t;
  } else if (t > .29 && t <= .447) {
    m = m + 1024395.67006;
    m = m - 22826983.8411*t;
    m = m + 221695566.585*tp2;
    m = m - 1225731636.28*tp2*t;
    m = m + 4219889881.74*tp4;
    m = m - 9263931021.04*tp4*t;
    m = m + 12664879027.3*tp4*tp2;
    m = m - 9858823353.65*tp4*tp2*t;
    m = m + 3345877714.47*tp4*tp4;
  } else if (t>= 0.447 && t <= 3.24) {
    m = m - 69.9540656568;
    m = m + 540.489990705*t;
    m = m - 1785.45996729*tp2;
    m = m + 3267.21333796*tp2*t;
    m = m - 3721.01017711*tp4;
    m = m + 2776.66619261*tp4*t;
    m = m - 1381.69750895*tp4*tp2;
    m = m + 454.445793525*tp4*tp2*t;
    m = m - 94.8530920783*tp4*tp4;
    m = m + 11.3775633348*tp4*tp4*t;
    m = m - 0.597112042033*tp4*tp4*tp2;
  } else if (t >= 3.24 && t <= 7.28) {
    m = m + 14.8304028947;
    m = m - 16.0726787694*t;
    m = m + 4.55554172673*tp2;
    m = m - 0.522521195039*tp2*t;
    m = m + 0.0214246750892*tp4;
  } else if (t >= 7.28 && t <= 7.8) {
    m = m + 3533.50268935;
    m = m - 1430.45985176*t;
    m = m + 192.934935096*tp2;
    m = m - 8.67530777079*tp2*t;
  } else if (t >= 7.8) {
    m = m - 3.5441266853;
    m = m - 0.0101229134154*t;
  } // End of delay time cases

  // Find mass ejection rate in solar masses per year
  m = std::pow(10.0,m);

  return m;
}

double s99_wind_energy(const double &t, const double &tsoon){
  // computes wind energy in cm^2/s^2 of 1.0E6 Msun cluster
  // at some time t (10^7 yr)
  double wind_energy = 0.0;
  const double tp2 = t*t;
  const double tp4 = tp2*tp2;

    // If the most massive star in this particle has a delay
    // time greater than 22 Myr (i.e. it is low mass)
    // turn on 10 km/s wind outflows.
  if ((tsoon < 0) or (tsoon > 2.20)){
    wind_energy = 12.6540102838;

  } else {

  if ( t < 0.29) {
    wind_energy = wind_energy + 16.7458830136;
    wind_energy = wind_energy + 2.87170884625*t;
    wind_energy = wind_energy - 260.188160495*tp2;
    wind_energy = wind_energy + 7588.41970548*tp2*t;
    wind_energy = wind_energy - 109128.119673*tp4;
    wind_energy = wind_energy + 834565.297424*tp4*t;
    wind_energy = wind_energy - 3488638.06781*tp4*tp2;
    wind_energy = wind_energy + 7534432.3913*tp4*tp2*t;
    wind_energy = wind_energy - 6577304.157*tp4*tp4;
  } else if( t >= 0.29 && t <= 0.518) {
    wind_energy = wind_energy-8756108.90226;
    wind_energy = wind_energy+249183578.797*t;
    wind_energy = wind_energy-3210919258.32*tp2;
    wind_energy = wind_energy+24729596291.2*tp2*t;
    wind_energy = wind_energy-126485776877.0*tp4;
    wind_energy = wind_energy+451123199767.0*tp4*t;
    wind_energy = wind_energy-1.14486556168E12*tp4*tp2;
    wind_energy = wind_energy+2.067395301E12*tp4*tp2*t;
    wind_energy = wind_energy-2.60335249368E12*tp4*tp4;
    wind_energy = wind_energy+2.17720380795E12*tp4*tp4*t;
    wind_energy = wind_energy-1.08835588273E12*tp4*tp4*tp2;
    wind_energy = wind_energy+246366864343.0*tp4*tp4*tp2*t;
 } else if(t >= 0.518 && t <= 2.0) {
    wind_energy = wind_energy+300.659606389;
    wind_energy = wind_energy-2175.28137376*t;
    wind_energy = wind_energy+7038.17965731*tp2;
    wind_energy = wind_energy-12640.7809456*tp2*t;
    wind_energy = wind_energy+13818.3936865*tp4;
    wind_energy = wind_energy-9434.54106014*tp4*t;
    wind_energy = wind_energy+3935.34399667*tp4*tp2;
    wind_energy = wind_energy-918.140140181*tp4*tp2*t;
    wind_energy = wind_energy+91.8268783049*tp4*tp4;
  } else if(t >= 2.0 && t <= 3.23) {
    wind_energy = wind_energy+1955.41904193;
    wind_energy = wind_energy-4288.5933438*t;
    wind_energy = wind_energy+3935.44109106*tp2;
    wind_energy = wind_energy-1921.4747372*tp2*t;
    wind_energy = wind_energy+526.42694795*tp4;
    wind_energy = wind_energy-76.729393462*tp4*t;
    wind_energy = wind_energy+4.64818353202*tp4*tp2;
  } else if(t >= 3.23) {
    wind_energy = wind_energy+12.6540102838;
  } // End of delay time cases.
} // end tsoon check

  // in cm^2/^2
  wind_energy = std::pow(10.0,wind_energy); // **(wind_energy -2.0*log10(v1))

  return wind_energy;
}

double s99_sn_mass(const double &t){
  // ejecta mass of supernova in Msun from s99

  double m_ejSN = 0.0;
  const double tp2 = t*t;
  const double tp4 = tp2*tp2;

  if (t < 0.513) {
    m_ejSN = m_ejSN + 3.40965833751;
    m_ejSN = m_ejSN - 16.0346449798*t;
    m_ejSN = m_ejSN + 31.5091825735*tp2;
    m_ejSN = m_ejSN - 21.3218283568*tp2*t;
  } else if(.513 <= t && t <= .918) {
    m_ejSN = m_ejSN - 314538.854117;
    m_ejSN = m_ejSN + 4453582.08399*t;
    m_ejSN = m_ejSN - 28218211.3741*tp2;
    m_ejSN = m_ejSN + 105370186.068*tp2*t;
    m_ejSN = m_ejSN - 256824281.305*tp4;
    m_ejSN = m_ejSN + 426986197.681*tp4*t;
    m_ejSN = m_ejSN - 490461521.485*tp4*tp2;
    m_ejSN = m_ejSN + 384394390.035*tp4*tp2*t;
    m_ejSN = m_ejSN - 196752045.251*tp4*tp4;
    m_ejSN = m_ejSN + 59399337.5861*tp4*tp4*t;
    m_ejSN = m_ejSN - 8033095.66643*tp4*tp4*tp2;
  } else if (0.918 <= t && t <= 3.23) {
    m_ejSN = m_ejSN + 1.74261906723;
    m_ejSN = m_ejSN - 0.92589554122*t;
    m_ejSN = m_ejSN + 0.551250718292*tp2;
    m_ejSN = m_ejSN - 0.220085806978*tp2*t;
    m_ejSN = m_ejSN + 0.0510662546479*tp4;
    m_ejSN = m_ejSN - 0.00504400687495*tp4*t;
  } else if (t >= 3.23) {
    m_ejSN = m_ejSN + 2.67991943387;
    m_ejSN = m_ejSN - 0.461075452846*t;
    m_ejSN = m_ejSN - 0.0326899620754*tp2;
  } // End of delay time cases

  m_ejSN = std::pow(10.0, m_ejSN);

  return m_ejSN;
}

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
  stencil_rad_               = ( (int) ((stencil_ - 1) / 2.0));
  number_of_feedback_cells_  = stencil_ * stencil_ * stencil_;
  shift_cell_center_         = enzo_config->method_feedback_shift_cell_center;


  use_ionization_feedback_   = enzo_config->method_feedback_use_ionization_feedback;
  // Do error checking here to make sure all required
  // fields exist..
  // NOTE: Good idea - make ALL method objects have a 'required fields' list
  //       defined in the header file and { loop through this in a more general
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

  double current_time  = block->time();

  // apply feedback depending on particle type
  // for now, just do this for all star particles

  int it = particle.type_index("star");
  int count = 0;

  // Polynomial coefficients for the fit to the delay time distribution
  double p_delay[6] = {
    3505021.4516666,         16621326.48066255,        4382816.59085307,
    46194173.14420852,      -52836941.28241706,       22967062.02780452,
  };

  // Polynomial coefficients for the fit to the stellar mass distribution
  // as a function of delay time.
  double p_mass[10] = {
    4.42035634891,        -13.2890466089,        26.1103296098,        -30.1876007562,
    21.8976126631,        -10.2544493943,        3.09621304958,        -0.581870413299,
    0.0618795946119,      -0.00284352366041
  };

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


/*
        // negative lifetime are particles that have alreahy gone SN
        // creation time must be > 0
        // only go SN if age >= lifetime
        if ( (plifetime[ipdl] <= 0.0) || (pcreation[ipdc] <= -100.0) ||
             (current_time - pcreation[ipdc]) < plifetime[ipdl]) continue;
        count++;

        // Update particle properties here if needed
        plifetime[ipdl] *= -1.0;                  // set to negative - flag for alreahy gone SNe
        pmass[ipdm] = pmass[ipdm] - total_ejecta_mass_;

*/

// --------------------------------------- start checking if explosion occurs ---

        int explosion_flag = -1, will_explode = -1;
        double soonest_explosion = -1.0, s49_tot = -1.0, td7 = -1.0;
        unsigned long int rand_int;


        if ( (plifetime[ipdl] <= 0.0) ||
             ( (current_time - pcreation[ipdc]) < 3.7e7 * cello::yr_s / enzo_units->time() ) ) {

          soonest_explosion = -1.0;

          // seed random number generator
          // using function to generate unique number:
          //     (a,b) = { a*a + a + b, a >= b;
          //                   a + b*b, a  < b}
          // where b is batch number (keeps vals smaller... since np > nb usually)
          //  WARNGING - THIS DOESN'T ACTUALLY WORK SINCE PARTICLES MOVE AROUND
          srand(( (ip >= ib) ? (ip*ip+ip+ib) : (ip+ib*ib) ));

          explosion_flag = -1;
          will_explode   =  0;

          // figure out how many explosions this particle undergoes.

          // The number of SNe expected from a 10^6 mass cluster (s99)
          const float expected_sn_s99 = 10616.955572;
          // Number of expected SNe for the minimum star particle mass
          const float lambda         = expected_sn_s99 * 1.0E-6 *
                                       enzo_config->method_star_maker_minimum_star_mass;

          // Knuth algorithm for generating a Poisson distribution
          // See http://goo.gl/sgLPcj

          float L = exp(-lambda), p = 1.0;
          int   k = 0;
          while (p>L){
            ++k;
            rand_int = (double(rand())) / (double(RAND_MAX));;
            float u      = float(rand_int%32768) / 32768.0;
            p*u;
          }

          // Now k-1 is approx Poisson(lambda)
          int number_of_sn = k - 1;

          if (number_of_sn == 0){
            explosion_flag = 0;
            continue;
          } else if (number_of_sn < 0){
            std::cout <<"Number of SN in distributed feedback is negative \n";
          }

          // If there are explosions, we need to check the ionizing luminosity
          // of the progenitor and know whether or not it will explode
          // right now, or at a later time

          // Loop over these SN
          for (int kk=0; kk<number_of_sn; kk++){
            // Draw delay time for the event
            rand_int = (double(rand())) / (double(RAND_MAX));
            double  x = float(rand_int%32768) / 32768.0;
            double delay_time = p_delay[0] + p_delay[1]*x + p_delay[2]*x*x +
                               p_delay[3]*x*x*x + p_delay[4]*x*x*x*x + p_delay[5]*x*x*x*x*x; // yr

            td7 = delay_time*1.0E7;
            double progenitor_mass = p_mass[0] + pmass[1]*td7 + pmass[2]*td7*td7 +
                                    p_mass[3]*td7*td7*td7 + p_mass[4]*td7*td7*td7*td7 +
                                    p_mass[5]*td7*td7*td7*td7*td7 +
                                    p_mass[6]*td7*td7*td7*td7*td7*td7 +
                                    p_mass[7]*td7*td7*td7*td7*td7*td7*td7 +
                                    p_mass[8]*td7*td7*td7*td7*td7*td7*td7*td7 +
                                    p_mass[9]*td7*td7*td7*td7*td7*td7*td7*td7*td7;
            progenitor_mass = std::pow(10.0, progenitor_mass); // in Msun

            delay_time   *=  cello::yr_s / enzo_units->time(); // code units
            double relative_time = current_time - pcreation[ipdc]; // code units
            if ((delay_time > relative_time) && (delay_time < relative_time + enzo_block->dt  )) {
              if (explosion_flag == -1){
                explosion_flag = 1;
              } else if (explosion_flag == 0){
                // fail
                ERROR("EnzoMethodDistributedFeedback::compute()",
                      "There was supposed to be no explosions, but this particle is flagged");
              } else if (explosion_flag >0){
                explosion_flag += 1; // add to the counters - this is not implemented behavior
              } else{
                ERROR("EnzoMethodDistributedFeedback::compute()",
                      "Somehow explosion_flag is below -1");
              }

            } // end check delay time

            if (relative_time < delay_time){
              // SN has not yet gone off. Get ionizing luminosity and
              // sum for the particle

              s49_tot += s49Lookup(progenitor_mass); // in 10^49 photons s^-1
              will_explode = 1;
              if ( (soonest_explosion < 0) || (delay_time < soonest_explosion) ){
                soonest_explosion = delay_time;
              }

            } // end check relative_time

//          Behavior not yet possible in Enzo-E
//            if ( relative_time < delay_time + 0.1 * cello::Myr_s / enzo_units->time()){
//              // change type to must refine
//            } else {
//              change back to star
//            }

          } // end loop over SN

          if (explosion_flag == -1) explosion_flag = 0;
        } // end lifetime check



     // ----------------------------------------
     // now compute ejecta mass and stuff from s99
     //    td7 is delay time in units of 10^7 years

        double wind_mass   = s99_wind_mass(td7);
     //          This ejection rate is correct for the 10^6 solar mass 'cluster'
     //           used to compute these rates with starburst 99.  Reduce to account
     //           for the size of the star particle
        wind_mass = wind_mass * enzo_config->method_star_maker_minimum_star_mass * 1.0E-6;
        wind_mass = wind_mass / cello::yr_s; // in Msun / s
        wind_mass = (wind_mass * enzo_units->time()) * enzo_block->dt; // Msun this timestep

        double tsoon7      = soonest_explosion * enzo_units->time() / (1.0E7 * cello::yr_s);
        double wind_energy = s99_wind_energy(td7, tsoon7); // in cm^2/s^2 (i.e. per unit mass in cgs)
        wind_energy        = wind_energy * wind_mass * cello::mass_solar; // now total E in erg

        double sn_mass = 0.0, sn_energy = 0.0;
        if (explosion_flag > 0){
          sn_mass = s99_sn_mass(td7) * explosion_flag; // sn mass in Msun
          sn_energy = 1.0E51 * explosion_flag;         // sn energy in erg
        }

        double m_eject = wind_mass + sn_mass;          // total mass in Msun
        double energy  = wind_energy + sn_energy;      // total energy in erg

        if (m_eject*cello::mass_solar / enzo_units->mass() > pmass[ipdm]){
          std::cout << "WARNING: S99 Distributed Feedback is loosing too much mass\n";
          std::cout << "setting particle mass to zero, but continuing anyway\n";
          std::cout << pmass[ipdm] << m_eject*cello::mass_solar/enzo_units->mass() << wind_mass << sn_mass << "\n";
          std::cout << will_explode << explosion_flag << "\n";
        } // mass will be removed elsewhere

// ---------------------------------------- end determining explosion properties ----


        // get corresponding grid position of particle
        // and shift it if it is too close to the grid boundaries
        double xpos = px[ipdp];
        double ypos = py[ipdp];
        double zpos = pz[ipdp];

        this->add_ionization_feedback(block, xpos, ypos, zpos,
                                      s49_tot, will_explode);

// -----

        this->inject_feedback(block, xpos, ypos, zpos,
                              m_eject, (energy/1.0E51), enzo_config->method_feedback_ke_fraction,
                              pvx[ipdv], pvy[ipdv], pvz[ipdv]);
        // remove mass - error checking on the std::max is handled above with warning
        pmass[ipdm] = std::max( 0.0, pmass[ipdm] - m_eject*cello::mass_solar/enzo_units->mass());
/*
        this->inject_feedback(block, xpos, ypos, zpos,
                              enzo_config->method_feedback_ejecta_mass,
                              enzo_config->method_feedback_supernova_energy,
                              enzo_config->method_feedback_ke_fraction,
                              pvx[ipdv], pvy[ipdv], pvz[ipdv] );
*/

        count++;
      } // end loop over particles
    } // end loop over batches

    if (count > 0){
      std::cout << "Number of feedback particles:   " << count << "\n";
    }

  } // end particle check

  return;
}

void EnzoMethodDistributedFeedback::add_ionization_feedback(
                                                        Block * block,
                                                        double xpos, double ypos, double zpos,
                                                        const double & s49_tot,
                                                        const int & will_explode){

  if (!(this->use_ionization_feedback_)) return;

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  Field field = block->data()->field();
  // Obtain grid sizes and ghost sizes

  enzo_float * d           = (enzo_float *) field.values("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");

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

  double cell_volume = (hx*hy*hz);
  double inv_volume  = 1.0 / cell_volume;

  //  stencil_rad_ is integer separation from cell center
  //  and edge of injection region (i.e. 1 for 3x3 injection grid)
  if (  shift_cell_center_ &&
       ( ((xpos - (stencil_rad_+1)*hx) < xm) || // note: xm/xp is min/max without including ghost
         ((xpos + (stencil_rad_+1)*hx) > xp) || //     +1 b/c of CIC inerpolation onto grid
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
  // this must account for ghost zones
  double xcell = (xpos - xm) / hx + gx - 0.5;
  double ycell = (ypos - ym) / hy + gy - 0.5;
  double zcell = (zpos - zm) / hz + gz - 0.5;

  int ix       = ((int) floor(xcell + 0.5));
  int iy       = ((int) floor(ycell + 0.5));
  int iz       = ((int) floor(zcell + 0.5));

  int index =  INDEX(ix,iy,iz,mx,my);

  // Case B recombination, assuming T = 10^4 K
  const float alpha = 2.60E-13; // cm^3 / s

  // AE: possibly actually compute mu from species fields
  double ndens = d[index] / enzo_config->ppm_mol_weight / cello::mass_hydrogen;

  double stromgren_radius = std::pow( (3.0 * s49_tot * 1.0E49) /
                            (4.0 * cello::pi * alpha * ndens*ndens),1.0/3.0);
  double stromgren_volume = (4.0/3.0)*cello::pi*stromgren_radius*stromgren_radius*stromgren_radius;

  double ionized          = cello::kboltz * 1.0E4 / enzo_config->ppm_mol_weight /
                            cello::mass_hydrogen / (enzo_units->length()*enzo_units->length()) *
                            enzo_units->time() * enzo_units->time();

  if (stromgren_volume <= cell_volume){
    ionized = ionized * stromgren_volume * inv_volume;
  }

  if ((will_explode == 1) && (ge[index] < ionized)){
    double diff = ionized - ge[index];
    ge[index]   = ge[index] + diff;
    te[index]   = te[index] + diff;
  }

  return;
}

void EnzoMethodDistributedFeedback::inject_feedback(
                                          Block * block,
                                          double xpos, double ypos, double zpos,
                                          double m_eject,  // in Msun
                                          double E_51, // in 10^51 erg
                                          double ke_fraction,
                                          enzo_float pvx,    //default -9999
                                          enzo_float pvy,    //default -9999
                                          enzo_float pvz){   //default -9999
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

  double density_per_cell = mass_per_cell * inv_vol;
  double energy_density_per_cell = energy_per_cell * inv_vol;


  //  stencil_rad_ is integer separation from cell center
  //  and edge of injection region (i.e. 1 for 3x3 injection grid)
  if (  shift_cell_center_ &&
       ( ((xpos - (stencil_rad_+1)*hx) < xm) || // note: xm/xp is min/max without including ghost
         ((xpos + (stencil_rad_+1)*hx) > xp) || //     +1 b/c of CIC inerpolation onto grid
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
  // this must account for ghost zones
  double xcell = (xpos - xm) / hx + gx - 0.5;
  double ycell = (ypos - ym) / hy + gy - 0.5;
  double zcell = (zpos - zm) / hz + gz - 0.5;

  int ix       = ((int) floor(xcell + 0.5));
  int iy       = ((int) floor(ycell + 0.5));
  int iz       = ((int) floor(zcell + 0.5));

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

  enzo_float *u_local=0, *v_local=0, *w_local=0, *d_local=0, *ge_local=0, *te_local=0;
  enzo_float *ke_before=0, *metal_local=0;

  // number of feedback cells + 1 cell in each dimension
  int num_loc = (stencil_ + 1) * (stencil_ + 1) * (stencil_ + 1);

  u_local     = new enzo_float[num_loc];
  v_local     = new enzo_float[num_loc];
  w_local     = new enzo_float[num_loc];
  d_local     = new enzo_float[num_loc];
  ge_local    = new enzo_float[num_loc];
  te_local    = new enzo_float[num_loc];
  metal_local = new enzo_float[num_loc];
  ke_before   = new enzo_float[num_loc];

  // assign initial values to these
  for (int i = 0; i < num_loc; i++){
    u_local[i] = 0.0; v_local[i] = 0.0; w_local[i] = 0.0;
    d_local[i] = 0.0; te_local[i] = 0.0; ge_local[i] = 0.0;
    ke_before[i] = 0.0; metal_local[i] = 0.0;
  }

  double ke_f = 0.0;

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
          if ( (index < 0) || (index >= mx*my*mz)) continue;

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

  double E_therm = (1.0 - ke_f) * energy_density_per_cell;

  // compute kinetic energy in the localized region on the grid before
  // the explosion
  int loc_index = 0;
  for (int k = iz - stencil_rad_ ; k <= iz + stencil_rad_ + 1; k++){
    for (int j = iy - stencil_rad_ ; j <= iy + stencil_rad_ + 1; j++){
      for (int i = ix - stencil_rad_ ; i <= ix + stencil_rad_ + 1; i++){

        int index = INDEX(i,j,k,mx,my);

        if ( (index < 0) || (index >= mx*my*mz)){
          loc_index++;
          continue;
        }

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
                         mx, my, mz,
                         ix, iy, iz, 1);

  // compute the total mass and energy in the cells before the explosion
  double sum_mass_init, sum_energy_init, sum_ke_init;
  this->sum_mass_energy(v3[0], v3[1], v3[2], d, ge, te,
                        mx, my, mz, ix, iy, iz,
                        sum_mass_init, sum_energy_init, sum_ke_init);

  // compute the mass and momentum properties of adding feedback to the
  // psuedo-grid centered on the particle position before doing so for
  // the real grid. this is used to compute coefficient properties below
  // and prep for the CIC deposition
  this->add_feedback_to_grid(u_local, v_local, w_local, d_local,
                             ge_local, te_local, metal_local,
                             stencil_+1, stencil_+1, stencil_+1,          // mx,my,mz for local grid
                             stencil_rad_, stencil_rad_, stencil_rad_, // local grid cell center  - should be 1 for 3x3x3 stencil (0,1,2) - 2 for 5x5x5 stencil (0,1, 2, 3,4)
                             dxc, dyc, dzc,
                             density_per_cell, 1.0, 0.0);

  // momenum injection - compute coefficients for injection
  double mom_per_cell = 0.0;
  if (ke_f > 0){
    double A=0.0, B=0.0, C=0.0;
    this->compute_coefficients( v3[0], v3[1], v3[2], d, ge,
                                u_local, v_local, w_local, d_local,
                                mx, my, mz, ix, iy, iz, A, B, C);

    A            = A - (sum_ke_init + ke_f*energy_density_per_cell);
    mom_per_cell = (-B + std::sqrt(B*B - 4.0 * A * C)) / (2.0 * C);
  } // add switch here? to do just momentum injection but no KE

  // Need to deposit metal(s) here

  this->add_feedback_to_grid(v3[0], v3[1], v3[2], d, ge, te, metal,
                             mx, my, mz, ix, iy, iz, dxc, dyc, dzc,
                             density_per_cell, mom_per_cell, E_therm);

  double sum_mass_final, sum_energy_final, sum_ke_final;

  this->sum_mass_energy(v3[0], v3[1], v3[2], d, ge, te,
                        mx, my, mz, ix, iy, iz,
                        sum_mass_final, sum_energy_final, sum_ke_final);

  // AE NOTE:
  //    should do some error checking here
  //    with kinetic energy and momentum

  // Now convert momentum back to velocity
  this->convert_momentum(v3[0], v3[1], v3[2], d,
                         pvx, pvy, pvz,
                         mx, my, mz, ix, iy, iz, 0);

  // Adjust total energy

  double ke_injected = 0.0;
  double delta_ke    = 0.0;
  double ke_after    = 0.0;
  loc_index = 0;
  for (int k = iz - stencil_rad_; k <= iz + stencil_rad_ + 1; k++){
    for (int j = iy - stencil_rad_; j <= iy + stencil_rad_ + 1; j++){
      for (int i = ix - stencil_rad_; i <= ix + stencil_rad_ + 1; i++){

        int index =  INDEX(i,j,k,mx,my);

        if (index < 0 || index >= mx*my*mz){
          loc_index++;
          continue;
        }

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


  delete [] u_local; u_local = NULL;
  delete [] v_local; v_local = NULL;
  delete [] w_local; w_local = NULL;
  delete [] d_local; d_local = NULL;
  delete [] ge_local; ge_local = NULL;
  delete [] te_local; te_local = NULL;
  delete [] metal_local; metal_local = NULL;
  delete [] ke_before; ke_before = NULL;

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

        if ( (index < 0) || (index >= mx*my*mz)) continue;

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

        if ( (index < 0) || (index >= mx*my*mz)) continue;

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

              if ( (index < 0) || (index >= mx*my*mz)) continue;

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
