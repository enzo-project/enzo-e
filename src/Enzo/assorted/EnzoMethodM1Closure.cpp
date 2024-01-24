// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodM1Closure.cpp
/// @author   William Hicks (whicks@ucsd.edu)
/// @date     Mon Aug 16 17:05:23 PDT 2021
/// @brief    Implements the EnzoMethodM1Closure class

#include "cello.hpp"

#include "enzo.hpp"

//#define DEBUG_PRINT_REDUCED_VARIABLES
//#define DEBUG_PRINT_GROUP_PARAMETERS
//#define DEBUG_RECOMBINATION
//#define DEBUG_INJECTION
//#define DEBUG_TRANSPORT
//#define DEBUG_ATTENUATION

//----------------------------------------------------------------------

EnzoMethodM1Closure ::EnzoMethodM1Closure(const int N_groups)
  : Method()
    , N_groups_(N_groups)
    , ir_injection_(-1)
{

  const int rank = cello::rank();

  cello::define_field("photon_density"); // photon number density

  if (rank >= 1) {
    cello::define_field("flux_x");
  }                                         
  if (rank >= 2) {
    cello::define_field("flux_y");
  }
  if (rank >= 3) {
    cello::define_field("flux_z");
  }

  for (int i=0; i<N_groups_; i++) {
    std::string istring = std::to_string(i); 
    cello::define_field("photon_density_" + istring);
    if (rank >= 1) cello::define_field("flux_x_" + istring);
    if (rank >= 2) cello::define_field("flux_y_" + istring);
    if (rank >= 3) cello::define_field("flux_z_" + istring);
  }

  // define other fields accessed by this method
  cello::define_field_in_group ("HI_density",    "color");
  cello::define_field_in_group ("HII_density",   "color");
  cello::define_field_in_group ("HeI_density",   "color");
  cello::define_field_in_group ("HeII_density",  "color");
  cello::define_field_in_group ("HeIII_density", "color");
  cello::define_field_in_group ("e_density",     "color");

  cello::define_field("pressure");
  cello::define_field("temperature"); // needed for recombination rates

  if (rank >= 1) {
    cello::define_field("P00"); // elements of the radiation pressure tensor.
  }
  if (rank >= 2) {
    cello::define_field("P10");
    cello::define_field("P01");
    cello::define_field("P11");
  }
  if (rank >= 3) {
    cello::define_field("P02");
    cello::define_field("P12");
    cello::define_field("P20");
    cello::define_field("P21");
    cello::define_field("P22");
  }

  // fields for refresh+accumulate
  for (int i=0; i<N_groups_; i++) {
    std::string istring = std::to_string(i);
    cello::define_field("photon_density_" + istring + "_deposit");
  } 

  // Initialize default Refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("photon_density");

  if (rank >= 1) refresh->add_field("flux_x");
  if (rank >= 2) refresh->add_field("flux_y");
  if (rank >= 3) refresh->add_field("flux_z");

  for (int i=0; i<N_groups_; i++) {
    std::string istring = std::to_string(i); 
    refresh->add_field("photon_density_" + istring);
    if (rank >= 1) refresh->add_field("flux_x_" + istring);
    if (rank >= 2) refresh->add_field("flux_y_" + istring);
    if (rank >= 3) refresh->add_field("flux_z_" + istring);  
  }

  // Store frequency group attributes as ScalarData variables.
  // Variables with suffix "mL" store the numerators/denominator
  // of eqs. (B6)-(B8). 
  // mL = mass_star * luminosity_star 
  ScalarDescr * scalar_descr = cello::scalar_descr_double();
 
  // only three ionizable species (HI, HeI, HeII)
  // photodissociation cross sections for H2 are added
  // if method_m1_closure_H2_photodissociation = true
  int N_species_ = 3 + enzo::config()->method_m1_closure_H2_photodissociation; 
  for (int i=0; i<N_groups_; i++) {
    scalar_descr->new_value( eps_string(i) );
    scalar_descr->new_value(  mL_string(i) );
    scalar_descr->new_value( eps_string(i) + mL_string(i) );

    for (int j=0; j<N_species_; j++) {
      scalar_descr->new_value( sigN_string(i,j) );
      scalar_descr->new_value( sigE_string(i,j) );

      scalar_descr->new_value( sigN_string(i,j) + mL_string(i) );
      scalar_descr->new_value( sigE_string(i,j) + mL_string(i) );
    }
  }

  // ScalarData variable for storing uniform LWB photon density
  if (enzo::config()->method_m1_closure_lyman_werner_background) {
    scalar_descr->new_value("N_LWB");
  }
 
  // Initialize Refresh object for after injection step 
  ir_injection_ = add_refresh_();

  cello::simulation()->refresh_set_name(ir_injection_, name()+":injection");
  Refresh * refresh_injection = cello::refresh(ir_injection_);

  refresh_injection->set_accumulate(true);

  // don't need to add fluxes to injection refresh 
  // because only the photon_density fields are updated there

  for (int i=0; i<N_groups_; i++) {
    std::string istring = std::to_string(i); 

    refresh_injection->add_field_src_dst
       ("photon_density_"+istring+"_deposit", "photon_density_"+istring); 
  }

  // read in data tables
  M1_tables = new M1Tables();

  refresh_injection->set_callback(CkIndex_EnzoBlock::p_method_m1_closure_solve_transport_eqn()); 
}

M1Tables::M1Tables ()
{
  // Read in any data tables relevant to the M1 closure method
  const EnzoConfig * enzo_config = enzo::config();
  if (enzo_config->method_m1_closure_flux_function == "HLL") {
    read_hll_eigenvalues(enzo_config->method_m1_closure_hll_file);
  } 
}
//----------------------------------------------------------------------

void EnzoMethodM1Closure ::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | N_groups_;
  p | ir_injection_;
}

//----------------------------------------------------------------------

void EnzoMethodM1Closure::compute ( Block * block ) throw()
{

  // need to execute this method on ALL blocks (even non-leaves) because
  // there is a global reduction at the end of call_inject_photons().
  // Charm requires all members of the chare array to participate in 
  // global reductions. If I call compute_done here for non-leaf blocks,
  // they will still participate in the global sum at the end of inject_photons(), 
  // but they will also execute the callback function following the contribute() call.
  // This means the non-leaves would end up calling compute_done() twice which results
  // in skipped cycles
  
  const EnzoConfig * enzo_config = enzo::config();

  Field field = block->data()->field();
  int mx,my,mz;
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);

  const int m = mx*my*mz;

  EnzoBlock * enzo_block = enzo::block(block);

  int N_groups = enzo_config->method_m1_closure_N_groups;
  int N_species = 3 + enzo_config->method_m1_closure_H2_photodissociation;

  if (enzo_config->method_m1_closure_thermochemistry) {
    // compute the temperature
    EnzoComputeTemperature compute_temperature(enzo::fluid_props(),
                                               enzo_config->physics_cosmology);

    compute_temperature.compute(enzo_block);
  }

  Scalar<double> scalar = block->data()->scalar_double();

  // reset "mL" sums to zero
  // TODO: only do this once every N cycles, where N is a parameter
  for (int i=0; i<N_groups; i++) {
    *(scalar.value( scalar.index( eps_string(i) + mL_string(i) ))) = 0.0;
    *(scalar.value( scalar.index( mL_string(i) ) )) = 0.0;
    for (int j=0; j<N_species; j++) {
      *(scalar.value( scalar.index( sigN_string(i,j) + mL_string(i) ))) = 0.0;
      *(scalar.value( scalar.index( sigE_string(i,j) + mL_string(i) ))) = 0.0;
    }
  }

  // initialize deposition fields to zero
  for (int i=0; i<N_groups; i++) {
    std::string istring = std::to_string(i);
    enzo_float * N_i_d = (enzo_float *) field.values("photon_density_" + istring + "_deposit");
    for (int j=0; j<m; j++)
    {
      N_i_d[j] = 0.0; 
    }
  }

  //start photon injection step
  //This function will start the transport step after a refresh
  this->call_inject_photons(enzo_block);

#ifdef DEBUG_PRINT_GROUP_PARAMETERS
  for (int i=0; i<N_groups; i++) {
    for (int j=0; j<N_species; j++) {
      CkPrintf("[i,j] = [%d,%d]; sigN = %1.2e cm^2; sigE = %1.2e cm^2; eps = %1.2e eV\n",i,j,
               *(scalar.value( scalar.index( sigN_string(i,j) ))),
               *(scalar.value( scalar.index( sigE_string(i,j) ))),
               *(scalar.value( scalar.index(  eps_string(i)   ))) / enzo_constants::erg_eV);
    }
  }
#endif


  return;
}

//----------------------------------------------------------------------

double EnzoMethodM1Closure::timestep ( Block * block ) throw()
{
  Data * data = block->data();
  Field field = data->field();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  double h_min = std::numeric_limits<double>::max();
  if (rank >= 1) h_min = std::min(h_min,hx);
  if (rank >= 2) h_min = std::min(h_min,hy);
  if (rank >= 3) h_min = std::min(h_min,hz);

  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  double courant = this->courant();
  double clight_frac = enzo_config->method_m1_closure_clight_frac;
  return courant * h_min / (3.0 * clight_frac*enzo_constants::clight / enzo_units->velocity());
}

//-----------------------------------------------------------------------

double EnzoMethodM1Closure::integrate_simpson(double a, double b,
                    int n, // Number of intervals
                    std::function<double(double,double,double,int)> f, double v1, double v2, int v3) throw()
{
    // solve 1D integral using composite simpson's rule
    double h = (b - a) / n;

    // Internal sample points, there should be n - 1 of them
    double sum_odds = 0.0;
    for (int i = 1; i < n; i += 2)
    {
        sum_odds += f(a + i * h, v1, v2, v3);
    }
    double sum_evens = 0.0;
    for (int i = 2; i < n; i += 2)
    {
        sum_evens += f(a + i * h, v1, v2, v3);
    }

    return (f(a,v1,v2,v3) + f(b,v1,v2,v3) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}

double EnzoMethodM1Closure::planck_function(double nu, double T, double clight, int dependent_variable) throw()
{
  double prefactor = 0.0;
  switch (dependent_variable)
  {
    case 0: //no prefactor
       prefactor = 1.0;
    case 1: //photon density
       prefactor = 8*cello::pi*nu*nu / (clight*clight*clight); 
       break;
    case 2: //energy density
       prefactor = 8*cello::pi*enzo_constants::hplanck*nu*nu*nu / (clight*clight*clight);
       break;
  }
  
  return prefactor / ( std::exp(enzo_constants::hplanck*nu/(enzo_constants::kboltz*T)) - 1 );

}

//-------------------------- INJECTION STEP ----------------------------

double EnzoMethodM1Closure::get_star_temperature(double M) throw()
{
  double L, R;
  // mass-luminosity relations for main sequence stars
  if      (M/enzo_constants::mass_solar < 0.43) L = 0.23*pow(M/enzo_constants::mass_solar,2.3);
  else if (M/enzo_constants::mass_solar < 2   ) L =      pow(M/enzo_constants::mass_solar,4.0);
  else if (M/enzo_constants::mass_solar < 55  ) L =  1.4*pow(M/enzo_constants::mass_solar,3.5);
  else L = 32000*(M/enzo_constants::mass_solar);

  // mass-radius relations (not valid for particles representing stellar populations)
  if (M/enzo_constants::mass_solar < 1) R = pow(M/enzo_constants::mass_solar, 0.8);
  else R = pow(M/enzo_constants::mass_solar, 0.57);
  
  L *= enzo_constants::luminosity_solar; // solar luminosity in erg/s
  R *= 6.957e10; // solar radius in cm

  return pow( L/(4*cello::pi*R*R*5.670e-7), 0.25);
}

// ----------------------

double EnzoMethodM1Closure::get_radiation_custom(EnzoBlock * enzo_block, 
           double energy, double pmass, double plum, 
           double dt, double inv_vol, int igroup) throw()
{
  const EnzoConfig * enzo_config = enzo::config();

  Scalar<double> scalar = enzo_block->data()->scalar_double();

  // if either particle_luminosity parameter is set, give all particles the same luminosity,
  // otherwise just use whatever value was passed into this function ('luminosity' attribute)

  if (enzo_config->method_m1_closure_particle_luminosity >= 0.0) {
    plum = enzo_config->method_m1_closure_particle_luminosity; // erg/s
  }

  // only add energy fraction of radiation into this group according to SED                                           
  double plum_i = plum * enzo_config->method_m1_closure_SED[igroup] / (energy * enzo_constants::erg_eV); // converted to photons/s

  double mL = pmass*plum_i; 
 
  // loop through ionizable species
  int N_species = 3 + enzo_config->method_m1_closure_H2_photodissociation;
  
  for (int j=0; j<N_species; j++) {
    double sigma_j = sigma_vernier(energy,j); // cm^2   

    #ifdef DEBUG_INJECTION
      CkPrintf("MethodM1Closure::get_radiation_custom -- j = %d; energy = %f eV; sigma_j = %1.2e cm^2; mL = %1.2e \n", j, energy, sigma_j, mL);
    #endif
 
    *(scalar.value( scalar.index(sigN_string(igroup, j) + mL_string(igroup) ))) += sigma_j * mL;
    *(scalar.value( scalar.index(sigE_string(igroup, j) + mL_string(igroup) ))) += sigma_j * mL;
  }

  *(scalar.value( scalar.index(mL_string(igroup)) )) += mL;
  *(scalar.value( scalar.index( eps_string(igroup   ) + mL_string(igroup) ))) += energy*enzo_constants::erg_eV * mL;

  #ifdef DEBUG_INJECTION
    CkPrintf("MethodM1Closure::get_radiation_custom -- Ndot = %1.2e photons/s \n", plum_i);
  #endif

  return  plum_i * inv_vol * dt; 
}

// -------

double EnzoMethodM1Closure::get_radiation_blackbody(EnzoBlock * enzo_block,  
                   double pmass, double freq_lower, double freq_upper, double clight, 
                   double dt, double cell_volume, int igroup) throw()
{
  // Does all calculations in CGS
  const EnzoConfig * enzo_config = enzo::config();

  int n = 10; // number of partitions for simpson's method

  //need to use lambda expression to pass in planck_function() function as a parameter
  //because member function in c++ are automatically attached to the `this` pointer,
  //so you need to capture `this` in order for the compiler to know which function
  //to point to
  
  if (freq_lower == 0.0) freq_lower = 1.0; //planck function undefined at zero
                                         //1 Hz is a very small frequency compared to ~1e16 Hz
                            
             
  // Get temperature of star
  double T = enzo_config -> method_m1_closure_temperature_blackbody;
  if (T > 0.0) { 
    T = get_star_temperature(pmass);
  }

  int planck_case_N = 1;
  int planck_case_E = 2;
                                        
  double N_integrated = integrate_simpson(freq_lower,freq_upper,n, 
                 [this](double a, double b, double c, int d) {return planck_function(a,b,c,d);},
                 T,clight,planck_case_N); 
  double E_integrated = integrate_simpson(freq_lower,freq_upper,n, 
                 [this](double a, double b, double c, int d) {return planck_function(a,b,c,d);},
                 T,clight,planck_case_E);

  //----------

  double luminosity = N_integrated * cell_volume/dt; // photons per second
  double mL = pmass*luminosity; // cgs 

  //----------Calculate photon group attributes--------
  Scalar<double> scalar = enzo_block->data()->scalar_double(); 
 
  //eq. B3 ----> eps = int(E_nu dnu) / int(N_nu dnu)
  *(scalar.value( scalar.index( eps_string(igroup) + mL_string(igroup) ) ))
                                    +=
                E_integrated / N_integrated * mL;

  int N_species = 3 + enzo_config->method_m1_closure_H2_photodissociation;
  for (int j=0; j<N_species; j++) { // loop over ionizable species

    // eq. B4 ----> sigmaN = int(sigma_nuj * N_nu dnu)/int(N_nu dnu)
    *(scalar.value( scalar.index(sigN_string(igroup, j) + mL_string(igroup) ) )) 
                                    +=
           integrate_simpson(freq_lower,freq_upper,n, 
                [this,j](double nu, double b, double c, int d){
                  return sigma_vernier(enzo_constants::hplanck*nu / enzo_constants::erg_eV, j)
                                *planck_function(nu,b,c,d);
                },
                T,clight,planck_case_N) / N_integrated * mL;

    // eq. B5 ----> sigmaE = int(sigma_nuj * E_nu dnu)/int(E_nu dnu)
    *(scalar.value( scalar.index(sigE_string(igroup, j) + mL_string(igroup)) ))
                                    +=
           integrate_simpson(freq_lower,freq_upper,n, 
                [this,j](double nu, double b, double c, int d){
                  return sigma_vernier(enzo_constants::hplanck*nu / enzo_constants::erg_eV, j)
                                *planck_function(nu,b,c,d);
                },
                T,clight,planck_case_E) / E_integrated * mL;
  }
 
  *(scalar.value( scalar.index(mL_string(igroup)) )) += mL;

  #ifdef DEBUG_INJECTION
    CkPrintf("MethodM1Closure::get_radiation_blackbody -- [freq_lower, freq_upper] = [%1.2e, %1.2e], N_integrated = %1.2e cm^-3, T = %1.2e K, Ndot = %1.2e photons/s \n", 
                     freq_lower, freq_upper, N_integrated, T, luminosity);
  #endif

  // return photon density update 
  return N_integrated; 
}

// ----

void EnzoMethodM1Closure::inject_photons ( EnzoBlock * enzo_block, int igroup ) throw()
{
  // Solve dN_i/dt = Ndot^*_i
  // routine for identifying star particles and getting their grid position
  // copy/pasted from enzo_EnzoMethodFeedback::compute_()

  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();
  double munit = enzo_units->mass();
  double lunit = enzo_units->length();
  double tunit = enzo_units->time();
  double Nunit = enzo_units->photon_number_density();

  double f_esc = enzo_config->method_m1_closure_photon_escape_fraction;

  Field field = enzo_block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);
  
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);
 
  double xm,ym,zm;
  double xp,yp,zp;
  double hx,hy,hz;
  enzo_block->data()->lower(&xm,&ym,&zm);
  enzo_block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  double cell_volume = hx*hy*hz * enzo_units->volume(); 

  double dt = enzo_block->state()->dt() * tunit;

  // get relevant field variables
  enzo_float * N          = (enzo_float *) field.values(
        "photon_density_"+std::to_string(igroup));
  enzo_float * N_deposit  = (enzo_float *) field.values(
        "photon_density_"+std::to_string(igroup)+"_deposit");

  Particle particle = enzo_block->data()->particle();
  int it = particle.type_index("star");
  
  // if no stars, don't do anything
  if (particle.num_particles(it) == 0) return;

  const int ia_m = particle.attribute_index (it, "mass");
  const int ia_L = particle.attribute_index (it, "luminosity");
  const int ia_x = (rank >= 1) ? particle.attribute_index (it, "x") : -1;
  const int ia_y = (rank >= 2) ? particle.attribute_index (it, "y") : -1;
  const int ia_z = (rank >= 3) ? particle.attribute_index (it, "z") : -1;

  const int dm = particle.stride(it, ia_m);
  const int dp = particle.stride(it, ia_x);
  const int dL = particle.stride(it, ia_L);

  const int nb = particle.num_batches(it);

  // bin energies in eV
  double E_lower = (enzo_config->method_m1_closure_energy_lower)[igroup];
  double E_upper = (enzo_config->method_m1_closure_energy_upper)[igroup];
  double E_mean = (enzo_config->method_m1_closure_energy_mean)[igroup];

  // convert energies to frequency in Hz
  double freq_lower = E_lower * enzo_constants::erg_eV / enzo_constants::hplanck;
  double freq_upper = E_upper * enzo_constants::erg_eV / enzo_constants::hplanck;
  double clight = enzo_config->method_m1_closure_clight_frac * enzo_constants::clight;

  // which type of radiation spectrum to use
  std::string radiation_spectrum = enzo_config->method_m1_closure_radiation_spectrum; 
  for (int ib=0; ib<nb; ib++){
    enzo_float *px=0, *py=0, *pz=0;
    enzo_float *pmass=0, *plum=0;
   
    pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
    plum  = (enzo_float *) particle.attribute_array(it, ia_L, ib);

    px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

    int np = particle.num_particles(it,ib);

    // loop through particles within each batch
    for (int ip=0; ip<np; ip++) {
      int ipdp = ip*dp;
      int ipdm = ip*dm;
      int ipdL = ip*dL;

      double xp = px[ipdp];
      double yp = py[ipdp];
      double zp = pz[ipdp];

      // get 3D grid index for particle - account for ghost zones!!
      double ix = (int) std::floor((xp - xm) / hx) + gx;
      double iy = (int) std::floor((yp - ym) / hy) + gy;
      double iz = (int) std::floor((zp - zm) / hz) + gz;

      // now get index of this cell
      int i = INDEX(ix,iy,iz,mx,my);

      // deposit photons
      double pmass_cgs = pmass[ipdm]*enzo_units->mass();

      double dN;
      if (radiation_spectrum == "blackbody") {
        // This function will take in particle mass as a parameter, and will fit
        // and integrate over the SED to get the total injection rate. The current
        // implementation in RAMSES assumes that SED of gas stays constant.
        dN = get_radiation_blackbody(enzo_block, pmass_cgs, freq_lower, freq_upper, clight,
                                     dt, cell_volume, igroup);
      }
      else if (radiation_spectrum == "custom") {
        // This function samples a user-defined SED to get the injection rate into each group.
        // Uses the `luminosity` particle attribute to calculate the photon injection rate.
        // If `particle_luminosity` is set to something > 0, all particles will be given the same 
        // luminosity specified by that parameter
        double plum_cgs = plum[ipdL] * munit * lunit * lunit/ (tunit * tunit * tunit);
        dN = get_radiation_custom(enzo_block, E_mean, pmass_cgs, plum_cgs, 
                                  dt, 1/cell_volume, igroup);
      }
      dN *= f_esc / Nunit; // put back into code units
    
      // CIC deposit with cloud radius of 1 cell width
      double wx, wy, wz;
      for (int iz_ = iz-1; iz_ <= iz+1; iz_++) {

        double zcell = zm + (iz_+0.5 - gz)*hz;
        wz = std::max(1.0 - std::abs(zp - zcell) / hz, 0.0);

        for (int iy_ = iy-1; iy_ <= iy+1; iy_++) {

          double ycell = ym + (iy_+0.5 - gy)*hy;
          wy = std::max(1.0 - std::abs(yp - ycell) / hy, 0.0);

          for (int ix_ = ix-1; ix_ <= ix+1; ix_++) {
            // cell positions
            double xcell = xm + (ix_+0.5 - gx)*hx;
            wx  = std::max(1.0 - std::abs(xp - xcell) / hx, 0.0);

            int i_ = INDEX(ix_,iy_,iz_,mx,my);
            double dN_cic = wx*wy*wz*dN;
            N_deposit[i_] += dN_cic;
            // recall that the refresh machinery uses exchanged boundary data
            // from the field held by N_deposit to DIRECTLY update the field
            // held by N (i.e. so we may as well update the rest of N now)
            N[i_] += dN_cic;
          }
        }
      }
      
    } // end loop over particles
  }  // end loop over batches

} // end function

//---------------------------------------------------------------------

void M1Tables::read_hll_eigenvalues(std::string hll_file) throw()
{
    std::fstream inFile;
    inFile.open(hll_file, std::ios::in);

    ASSERT("M1Tables::read_hll_eigenvalues()", "hll_file failed to open!",
           inFile.is_open());

    // store table in vectors
    int line_count = 10201; // 101*101
    
    hll_table_f_.resize(line_count);
    hll_table_theta_.resize(line_count);
    hll_table_lambda_min_.resize(line_count);
    hll_table_lambda_max_.resize(line_count);
    hll_table_col3_.resize(line_count);
    hll_table_col4_.resize(line_count);

    int i = 0;
    while(inFile >> this->hll_table_f_[i] >> this->hll_table_theta_[i] >> 
                    this->hll_table_lambda_min_[i] >> 
                    this->hll_table_col3_[i] >> this->hll_table_col4_[i] >> 
                    this->hll_table_lambda_max_[i]) i++;

    inFile.close();
}



void EnzoMethodM1Closure::compute_hll_eigenvalues (double f, double theta, double * lmin, double * lmax, 
                                                    double clight) throw() 
{
  double lf = f*100;
  double lt = theta/cello::pi * 100;

  int i = std::min(int(lf),99);
  int j = std::min(int(lt),99);

  double dd1 = lf - i;
  double dd2 = lt - j;
  double de1 = 1 - dd1; 
  double de2 = 1 - dd2;

  *lmin = 0.0;
  *lmin += de1*de2*M1_tables->hll_table_lambda_min(i  ,j  );
  *lmin += dd1*de2*M1_tables->hll_table_lambda_min(i+1,j  );
  *lmin += de1*dd2*M1_tables->hll_table_lambda_min(i  ,j+1);
  *lmin += dd1*dd2*M1_tables->hll_table_lambda_min(i+1,j+1);

  *lmax = 0.0;
  *lmax += de1*de2*M1_tables->hll_table_lambda_max(i  ,j  );
  *lmax += dd1*de2*M1_tables->hll_table_lambda_max(i+1,j  );
  *lmax += de1*dd2*M1_tables->hll_table_lambda_max(i  ,j+1);
  *lmax += dd1*dd2*M1_tables->hll_table_lambda_max(i+1,j+1);
}

//----------------------------------------------------------------------

double EnzoMethodM1Closure::flux_function (double U_l, double U_lplus1,
					    double Q_l, double Q_lplus1, double clight,
                                            double lmin, double lmax,  
					    std::string type) throw()
{
  // returns face-flux of a cell at index idx
  if (type == "GLF") {
    return 0.5*(  Q_l+Q_lplus1 - clight*(U_lplus1-U_l) ); 
  }

  else if (type == "HLL") {
    return (lmax*Q_l - lmin*Q_lplus1 + lmax*lmin*clight*(U_lplus1-U_l)) / (lmax - lmin);
  }

  else {
    ERROR("EnzoMethodM1Closure::flux_function",
	     "flux_function type not recognized");
    return 0.0; 
  }
}

//--------------------------------------------------------------------------

double EnzoMethodM1Closure::deltaQ_faces (double U_l, double U_lplus1, double U_lminus1,
                                                  double Q_l, double Q_lplus1, double Q_lminus1,
                                                  double clight, double lmin, double lmax, 
                                                  std::string flux_type) throw()
{
  // calls flux_function(), and calculates Q_{i-1/2} - Q_{i+1/2}
  
  return flux_function(U_lminus1, U_l     , Q_lminus1, Q_l     , clight, lmin, lmax, flux_type) - 
         flux_function(U_l      , U_lplus1, Q_l      , Q_lplus1, clight, lmin, lmax, flux_type); 
}

//--------------------------------------------------------------------------

void EnzoMethodM1Closure::get_reduced_variables (double * chi, double (*n)[3], int i, double clight,
                               enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz) throw()
{
  double Fnorm = sqrt(Fx[i]*Fx[i] + Fy[i]*Fy[i] + Fz[i]*Fz[i]);
  double f = N[i] > 0 ? std::min(Fnorm / (clight*N[i] ), 1.0) : 0.0; // reduced flux ( 0 < f < 1)
  *chi = (3 + 4*f*f) / (5 + 2*sqrt(4-3*f*f)); // isotropy measure (1/3 < chi < 1)
  if (Fnorm > 0.0) {
    (*n)[0] = Fx[i]/Fnorm;
    (*n)[1] = Fy[i]/Fnorm; 
    (*n)[2] = Fz[i]/Fnorm;
  }

  #ifdef DEBUG_PRINT_REDUCED_VARIABLES
    CkPrintf("i = %d; N = %1.2e; Fx = %1.2e; Fy = %1.2e; Fz = %1.2e; f = %1.2e; chi = %1.2e; n=[%1.2e,%1.2e,%1.2e]\n", i, N[i], Fx[i], Fy[i], Fz[i], f, *chi, (*n)[0], (*n)[1], (*n)[2]); 
  #endif
}

//--------------------------------------------------------------------------

void EnzoMethodM1Closure::get_pressure_tensor (EnzoBlock * enzo_block, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz, double clight) 
                       throw()
{
  Field field = enzo_block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);

  // if rank >= 1
  enzo_float * P00 = (enzo_float *) field.values("P00");
  // if rank >= 2
  enzo_float * P10 = (enzo_float *) field.values("P10");
  enzo_float * P01 = (enzo_float *) field.values("P01");
  enzo_float * P11 = (enzo_float *) field.values("P11");
  // if rank >= 3
  enzo_float * P02 = (enzo_float *) field.values("P02");
  enzo_float * P12 = (enzo_float *) field.values("P12");
  enzo_float * P20 = (enzo_float *) field.values("P20");
  enzo_float * P21 = (enzo_float *) field.values("P21"); 
  enzo_float * P22 = (enzo_float *) field.values("P22"); 

  // Need to directly calculate pressure tensor elements 
  // one layer deep into the ghost zones because active cells
  // need information about their neighbors, and there's no 
  // guarantee that a neighboring block will have updated
  // its pressure tensor by the time this block starts
  //
  // Note that we're actually storing c^P, since that's the actual
  // value that's being converted to a flux 
  for (int iz=gz-1; iz<mz-gz+1; iz++) { 
   for (int iy=gy-1; iy<my-gy+1; iy++) {
    for (int ix=gx-1; ix<mx-gx+1; ix++) {
      int i = INDEX(ix,iy,iz,mx,my); //index of current cell
      double chi;
      double n[3] = {0.0, 0.0, 0.0}; 
      get_reduced_variables( &chi, &n, i, clight, 
                            N, Fx, Fy, Fz);
      double iterm = 0.5*(1.0-chi);   // identity term
      double oterm = 0.5*(3.0*chi-1); // outer product term
      double cc = clight * clight;
      P00[i] = cc * N[i] * (oterm *n[0]*n[0] + iterm );
      P10[i] = cc * N[i] *  oterm *n[1]*n[0];
      P01[i] = cc * N[i] *  oterm *n[0]*n[1];
      P11[i] = cc * N[i] * (oterm *n[1]*n[1] + iterm );
      P02[i] = cc * N[i] *  oterm *n[0]*n[2];
      P12[i] = cc * N[i] *  oterm *n[1]*n[2];
      P20[i] = cc * N[i] *  oterm *n[2]*n[0];
      P21[i] = cc * N[i] *  oterm *n[2]*n[1];
      P22[i] = cc * N[i] * (oterm *n[2]*n[2] + iterm );
    }
   }
  }

}

//--------------------------------------------------------------------------


void EnzoMethodM1Closure::get_U_update (EnzoBlock * enzo_block, double * N_update, 
                       double * Fx_update, double * Fy_update, double * Fz_update, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz,
                       double hx, double hy, double hz, double dt, double clight, 
                       int i, int idx, int idy, int idz) throw()
{
  Field field = enzo_block->data()->field();
  // if rank >= 1
  enzo_float * P00 = (enzo_float *) field.values("P00");
  // if rank >= 2
  enzo_float * P10 = (enzo_float *) field.values("P10");
  enzo_float * P01 = (enzo_float *) field.values("P01");
  enzo_float * P11 = (enzo_float *) field.values("P11");
  // if rank >= 3
  enzo_float * P02 = (enzo_float *) field.values("P02");
  enzo_float * P12 = (enzo_float *) field.values("P12");
  enzo_float * P20 = (enzo_float *) field.values("P20");
  enzo_float * P21 = (enzo_float *) field.values("P21"); 
  enzo_float * P22 = (enzo_float *) field.values("P22"); 

  // if using HLL flux function, compute eigenvalues here
  std::string flux_type = enzo::config()->method_m1_closure_flux_function;

  // HLL min and max eigenvalues
  // +/- clight corresponds to GLF flux function
  double lmin_x = -1.0, lmin_y = -1.0, lmin_z = -1.0;
  double lmax_x =  1.0, lmax_y =  1.0, lmax_z =  1.0;

  if (flux_type == "HLL") {
    double Fnorm = sqrt(Fx[i]*Fx[i] + Fy[i]*Fy[i] + Fz[i]*Fz[i]);
    double f = std::min(Fnorm / (N[i]*clight), 1.0);

    double theta_x = acos(std::min(Fx[i] / Fnorm, -1.0));
    double theta_y = acos(std::min(Fy[i] / Fnorm, -1.0));
    double theta_z = acos(std::min(Fz[i] / Fnorm, -1.0));

    compute_hll_eigenvalues(f, theta_x, &lmin_x, &lmax_x, clight);
    compute_hll_eigenvalues(f, theta_y, &lmin_y, &lmax_y, clight);
    compute_hll_eigenvalues(f, theta_z, &lmin_z, &lmax_z, clight);
  } 
  
  // if rank >= 1
  *N_update += dt/hx * deltaQ_faces( N[i],  N[i+idx],  N[i-idx],
                                   Fx[i], Fx[i+idx], Fx[i-idx], clight, lmin_x, lmax_x, flux_type );
    
  *Fx_update += dt/hx * deltaQ_faces(Fx[i], Fx[i+idx], Fx[i-idx],
                                   P00[i],
                                   P00[i+idx],
                                   P00[i-idx], clight, lmin_x, lmax_x, flux_type );

  // if rank >= 2
  *N_update += dt/hy * deltaQ_faces( N[i],  N[i+idy],  N[i-idy],
                                   Fy[i], Fy[i+idy], Fy[i-idy], clight, lmin_y, lmax_y, flux_type );

  *Fx_update += dt/hy * deltaQ_faces(Fx[i], Fx[i+idy], Fx[i-idy],
                                   P10[i],
                                   P10[i+idy],
                                   P10[i-idy], clight, lmin_y, lmax_y, flux_type );

  *Fy_update += dt/hx * deltaQ_faces(Fy[i], Fy[i+idx], Fy[i-idx],
                                   P01[i],
                                   P01[i+idx],
                                   P01[i-idx], clight, lmin_x, lmax_x, flux_type );

  *Fy_update += dt/hy * deltaQ_faces(Fy[i], Fy[i+idy], Fy[i-idy],
                                   P11[i],
                                   P11[i+idy],
                                   P11[i-idy], clight, lmin_y, lmax_y, flux_type );


   // if rank >= 3
  *N_update += dt/hz * deltaQ_faces( N[i],  N[i+idz],  N[i-idz],
                                   Fz[i], Fz[i+idz], Fz[i-idz], clight, lmin_z, lmax_z, flux_type );

  *Fx_update += dt/hz * deltaQ_faces(Fx[i], Fx[i+idz], Fx[i-idz],
                                   P20[i],
                                   P20[i+idz],
                                   P20[i-idz], clight, lmin_z, lmax_z, flux_type );

  *Fy_update += dt/hz * deltaQ_faces(Fy[i], Fy[i+idz], Fy[i-idz],
                                   P21[i],
                                   P21[i+idz],
                                   P21[i-idz], clight, lmin_z, lmax_z, flux_type);

  *Fz_update += dt/hx * deltaQ_faces(Fz[i], Fz[i+idx], Fz[i-idx],
                                   P02[i],
                                   P02[i+idx],
                                   P02[i-idx], clight, lmin_x, lmax_x, flux_type);

  *Fz_update += dt/hy * deltaQ_faces(Fz[i], Fz[i+idy], Fz[i-idy],
                                   P12[i],
                                   P12[i+idy],
                                   P12[i-idy], clight, lmin_y, lmax_y, flux_type);

  *Fz_update += dt/hz * deltaQ_faces(Fz[i], Fz[i+idz], Fz[i-idz],
                                   P22[i],
                                   P22[i+idz],
                                   P22[i-idz], clight, lmin_z, lmax_z, flux_type);

}

//----------------------------------

double EnzoMethodM1Closure::sigma_vernier (double energy, int type) throw()
{
  // copy FindCrossSection.C from Enzo
  // Uses fits from Vernier et al. (1996) to calculate photoionization cross-section 
  // between photons of energy E and gas of a given species. Then need to average this
  // value over  
 
  double sigma;
  double e_th, e_max, e0, sigma0, ya, P, yw, y0, y1;

  switch (type) { 
    // HI
  case 0:
    e_th = 13.6;
    e_max = 5e4;
    e0 = 4.298e-1;
    sigma0 = 5.475e4;
    ya = 32.88;
    P = 2.963;
    yw = y0 = y1 = 0.0;
    break;

    // HeI
  case 1:
    e_th = 24.59;
    e_max = 5e4;
    e0 = 13.61;
    sigma0 = 9.492e2;
    ya = 1.469;
    P = 3.188;
    yw = 2.039;
    y0 = 0.4434;
    y1 = 2.136;
    break;
  
    // HeII
  case 2:
    e_th = 54.42;
    e_max = 5e4;
    e0 = 1.720;
    sigma0 = 1.369e4;
    ya = 32.88;
    P = 2.963;
    yw = y0 = y1 = 0.0;
    break;
  
    // H2I (Lyman-Werner band)
  case 3:
    e_th = 11.18;
    e_max = 13.60;
    sigma0 = 3.71;
    break; 
  }

  // return 0 if below ionization threshold or above max energy
  if ((energy < e_th) || (energy > e_max)) return 0.0;

  double fy;
  if (type != 3) {
    double x, y;
    x = energy/e0 - y0;
    y = sqrt(x*x + y1*y1);
    fy = ((x-1.0)*(x-1.0) + yw*yw) * pow(y, 0.5*P-5.5) * 
      pow((1.0 + sqrt(y/ya)), -P);
  }
  else {
    fy = 1.0;
  }

  sigma = sigma0 * fy * 1e-18;

  return sigma;

}

//---------------------------------

void EnzoMethodM1Closure::get_photoionization_and_heating_rates (EnzoBlock * enzo_block, double clight) throw() 
{
  // Calculates photoionization and heating rates in each cell according to RAMSES-RT prescription
  // See pg. 14 of https://grackle.readthedocs.io/_/downloads/en/latest/pdf/ for relavent Grackle
  // parameters. 
  // ionization -- first term of eq. A21 -- sum_i(sigmaN*clight*Ni), where i iterates over frequency groups
  // ionization rates should be in code_time^-1
  // heating rates should be in erg s^-1 cm^-3 / nHI 
  // NOTE: Rates are zero if mean energy in a bin is less than the ionization threshold.

  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();
  double Nunit = enzo_units->photon_number_density();

  Field field = enzo_block->data()->field();
  Scalar<double> scalar = enzo_block->data()->scalar_double(); 

  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);
   
  double xm,ym,zm;
  double xp,yp,zp;
  enzo_block->lower(&xm,&ym,&zm);
  enzo_block->upper(&xp,&yp,&zp);

  enzo_float * HI_density    = (enzo_float *) field.values("HI_density");
  enzo_float * HeI_density   = (enzo_float *) field.values("HeI_density");
  enzo_float * HeII_density  = (enzo_float *) field.values("HeII_density");
  
  enzo_float * RT_HI_ionization_rate   = (enzo_float *) field.values("RT_HI_ionization_rate");
  enzo_float * RT_HeI_ionization_rate  = (enzo_float *) field.values("RT_HeI_ionization_rate");
  enzo_float * RT_HeII_ionization_rate = (enzo_float *) field.values("RT_HeII_ionization_rate");

  enzo_float * RT_heating_rate = (enzo_float *) field.values("RT_heating_rate");

  std::vector<enzo_float*> chemistry_fields = {HI_density, HeI_density, HeII_density};

  std::vector<enzo_float*> ionization_rate_fields = {RT_HI_ionization_rate, 
                               RT_HeI_ionization_rate, RT_HeII_ionization_rate};

  double tunit = enzo_units->time();
  double rhounit = enzo_units->density();

  std::vector<double> Eion = {13.6*enzo_constants::erg_eV, 24.59*enzo_constants::erg_eV, 54.42*enzo_constants::erg_eV};
  double mH = enzo_constants::mass_hydrogen;
  std::vector<double> masses = {mH,4*mH,4*mH};

  std::vector<enzo_float*> photon_densities = {};
  for (int igroup=0; igroup<enzo_config->method_m1_closure_N_groups; igroup++) { 
    photon_densities.push_back( (enzo_float *) field.values("photon_density_" + std::to_string(igroup))) ;
  }

  int N_species = 3; // HI, HeI, HeII
  // loop through cells
  for (int i=0; i<mx*my*mz; i++) {
    double nHI = std::max(HI_density[i] * rhounit / mH, 1e-20); // cgs 
    double heating_rate = 0.0; 
    for (int j=0; j<N_species; j++) { //loop over species
      double ionization_rate = 0.0;
      for (int igroup=0; igroup<enzo_config->method_m1_closure_N_groups; igroup++) { //loop over groups
        double sigmaN = *(scalar.value( scalar.index( sigN_string(igroup,j) ))); // cm^2 
        double sigmaE = *(scalar.value( scalar.index( sigE_string(igroup,j) ))); // cm^2
        double eps    = *(scalar.value( scalar.index(  eps_string(igroup  ) ))); // erg 

        double N_i = (photon_densities[igroup])[i] * Nunit; // cm^-3
        double n_j = (chemistry_fields[j])[i] * rhounit / masses[j]; //number density of species j
           
        double ediff = eps*sigmaE - Eion[j]*sigmaN;
            
        ionization_rate += sigmaN*clight*N_i;
        heating_rate    += std::max( n_j*clight*N_i*ediff, 0.0 ); // Equation A16

      }
      (ionization_rate_fields[j])[i] = ionization_rate * tunit; //update fields with new value, put ionization rates in 1/time_units
    }
        
  RT_heating_rate[i] = heating_rate / nHI; // units of erg/s/cm^3/nHI
  }

  // H2 photodissociation from LW radiation
  if (enzo_config->method_m1_closure_H2_photodissociation) {
    enzo_float * RT_H2_photodissociation_rate = (enzo_float *) field.values("RT_H2_dissociation_rate");
    for (int i=0; i<mx*my*mz; i++) {
      double N = (photon_densities[0])[i] * Nunit; // LW-group assumed to be group 0
      double sigmaN = *(scalar.value( scalar.index( sigN_string(0,3) ))); // cm^2 
      RT_H2_photodissociation_rate[i] = sigmaN*clight*N * tunit;
    }
  }
 
}

//---------------------------------

double EnzoMethodM1Closure::get_alpha (double T, int species, char rec_case) throw()
{
  // Return recombination rate coefficients according to appendix E2
  // Takes temperature, species, and recombination type (case A or B) 
  // and spits out the rate coefficient
  
  double lambda, lambda_0;
  double a,b,c,d;
  
  switch (rec_case) {
    case 'A': // case A recombination
      switch (species) {
        case 0: // HII + e -> HI + photon
          lambda = 315614 / T;
          lambda_0 = 0.522;
          a = 1.269e-13;
          b = 1.503;
          c = 0.47;
          d = 1.923;
          break;
        case 1: // HeII + e -> HeI + photon
          lambda = 570670 / T;
          return 3e-14 * pow(lambda,0.654);
        case 2: // HeIII + e -> HeII + photon
          lambda = 1263030 / T;
          lambda_0 = 0.522;
          a = 2.538e-13;
          b = 1.503;
          c = 0.47;
          d = 1.923;
          break;
      }    

      break;

    case 'B': // case B recombination
      switch (species) {
        case 0: // HII + e -> HI + photon
          lambda = 315614 / T;
          lambda_0 = 2.74;
          a = 2.753e-14;
          b = 1.5;
          c = 0.407;
          d = 2.242;
          break;
        case 1: // HeII + e -> HeI + photon
          lambda = 570670 / T;
          return 1.26e-14 * pow(lambda,0.75);
        case 2: // HeIII + e -> HeII + photon
          lambda = 1263030 / T;
          lambda_0 = 2.74;
          a = 5.506e-14;
          b = 1.5;
          c = 0.407;
          d = 2.242;
          break;
      }    

      break;

  } 
  return a * pow(lambda,b) * pow( 1+pow(lambda/lambda_0,c), -d); 
}

//---------------------------------

int EnzoMethodM1Closure::get_b_boolean (double E_lower, double E_upper, int species) throw()
{
  // boolean 1 or 0 which specifies whether or not photon from given recombination
  // lies within given energy range

  switch (species) { // species after recombination
    case 0: // HI
      if ( (E_lower <= 13.6 ) && (13.6  < E_upper)) return 1;
      else break;
    case 1: // HeI
      if ( (E_lower <= 24.59) && (24.59 < E_upper)) return 1;
      else break;
    case 2: // HeII
      if ( (E_lower <= 54.42) && (54.42 < E_upper)) return 1;
      else break; 
  }
  
  // if energy not in this group, return 0
  return 0;
}

//---------------------------------

double EnzoMethodM1Closure::C_add_recombination (EnzoBlock * enzo_block,
                                                enzo_float * T, int i, int igroup,
                                                double E_lower, double E_upper) throw()
{
  // update photon_density to account for recombination
  // 2nd half of eq 25, using backwards-in-time quantities for all variables.
  // this is called once for each group.

  Field field = enzo_block->data()->field();
  if (! field.is_field("density")) return 0.0;

  EnzoUnits * enzo_units = enzo::units();
  double rhounit = enzo_units->density();
  double Cunit = enzo_units->photon_number_density() / enzo_units->time(); 

  enzo_float * e_density = (enzo_float *) field.values("e_density");

  std::vector<std::string> chemistry_fields = {"HI_density", 
                                               "HeI_density", "HeII_density"};
  double mH  = enzo_constants::mass_hydrogen; // cgs

  std::vector<double> masses = {mH,4*mH, 4*mH};

  double C = 0.0;
  for (std::size_t j=0; j<chemistry_fields.size(); j++) {  
    enzo_float * density_j = (enzo_float *) field.values(chemistry_fields[j]);
     
    int b = get_b_boolean(E_lower, E_upper, j);

    double alpha_A = get_alpha(T[i], j, 'A');  // cgs
    double alpha_B = get_alpha(T[i], j, 'B');

    double n_j = density_j[i]*rhounit/masses[j];
    double n_e = e_density[i]*rhounit/mH; // electrons have same mass as protons in code units

    C += b*(alpha_A-alpha_B) * n_j*n_e / Cunit;

#ifdef DEBUG_RECOMBINATION
    CkPrintf("MethodM1Closure::C_add_recombination -- j=%d; alpha_A = %1.3e; alpha_B = %1.3e; n_j = %1.3e; n_e = %1.3e; b_boolean = %d\n", j, alpha_A, alpha_B, n_j, n_e, b);
#endif
  }

#ifdef DEBUG_RECOMBINATION
  CkPrintf("MethodM1Closure::C_add_recombination -- [E_lower, E_upper] = [%.2f, %.2f]; dN_dt[i] = %1.3e; dt = %1.3e\n", E_lower, E_upper, C, enzo_block->dt);
#endif 

  return C;
}

//---------------------------------

double EnzoMethodM1Closure::D_add_attenuation ( EnzoBlock * enzo_block, 
                                             double clight, int i, int igroup) throw()
{
  // Attenuate radiation

  EnzoUnits * enzo_units = enzo::units();
  double rhounit = enzo_units->density();
  double tunit = enzo_units->time();

  Field field = enzo_block->data()->field();

  if (! field.is_field("density")) return 0.0;
 
  std::vector<std::string> chemistry_fields = {"HI_density", 
                                               "HeI_density", "HeII_density"};

  double mH = enzo_constants::mass_hydrogen;
  std::vector<double> masses = {mH,4*mH, 4*mH};
 
  Scalar<double> scalar = enzo_block->data()->scalar_double();
  double D = 0.0;
  for (std::size_t j=0; j<chemistry_fields.size(); j++) {  
    enzo_float * density_j = (enzo_float *) field.values(chemistry_fields[j]);
    double n_j = density_j[i]*rhounit / masses[j];     
    double sigN_ij = *(scalar.value( scalar.index( sigN_string(igroup, j) )));
    
    D += n_j * clight*sigN_ij * tunit; // code_time^-1

    #ifdef DEBUG_ATTENUATION
      CkPrintf("[i,j]=[%d,%d]; sigN_ij=%1.2e; n_j=%1.2e; clight=%1.2e\n", igroup, j, sigN_ij, n_j, clight);
    #endif
  }

#ifdef DEBUG_ATTENUATION
    CkPrintf("i=%d; D=%e s^-1; dt = %e\n",i, D / tunit, enzo_block->dt*enzo_units->time());
#endif 
  return D; 
}

//----------------------

void EnzoMethodM1Closure::solve_transport_eqn ( EnzoBlock * enzo_block, int igroup ) throw()
{
  // Solve dU/dt + del[F(U)] = 0; F(U) = { (Fx,Fy,Fz), c^2 P }
  //                                U  = { N, (Fx,Fy,Fz) }
  // M1 closure: P_i = D_i * N_i, where D_i is the Eddington tensor for 
  // photon group i
 
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  Field field = enzo_block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);
   
  double xm,ym,zm;
  double xp,yp,zp;
  enzo_block->lower(&xm,&ym,&zm);
  enzo_block->upper(&xp,&yp,&zp);

  // array incremenent (because 3D array of field values are flattened to 1D)
  const int idx = 1;
  const int idy = mx;
  const int idz = mx*my; 

  // energy bounds for this group (leave in eV)
  double E_lower = enzo_config->method_m1_closure_energy_lower[igroup]; 
  double E_upper = enzo_config->method_m1_closure_energy_upper[igroup]; 
  
  std::string istring = std::to_string(igroup);
  enzo_float * N  = (enzo_float *) field.values("photon_density_" + istring);
  enzo_float * Fx = (enzo_float *) field.values("flux_x_" + istring);
  enzo_float * Fy = (enzo_float *) field.values("flux_y_" + istring);
  enzo_float * Fz = (enzo_float *) field.values("flux_z_" + istring);

  enzo_float * T = (enzo_float *) field.values("temperature");

  const int m = mx*my*mz;

  // extra copy of fields needed to store
  // the evolved values until the end
  enzo_float * Nnew  = new enzo_float[m]; 
  enzo_float * Fxnew = new enzo_float[m];
  enzo_float * Fynew = new enzo_float[m];
  enzo_float * Fznew = new enzo_float[m];

  double lunit = enzo_units->length();
  double tunit = enzo_units->time();
  double Nunit = enzo_units->photon_number_density();

  double dt = enzo_block->state()->dt();
  double hx = (xp-xm)/(mx-2*gx);
  double hy = (yp-ym)/(my-2*gy);
  double hz = (zp-zm)/(mz-2*gz);
  double clight_cgs = enzo_config->method_m1_closure_clight_frac*enzo_constants::clight; 
  double clight_code = clight_cgs * tunit/lunit;
  
  for (int i=0; i<m; i++)
  {
    Nnew [i] = N [i];
    Fxnew[i] = Fx[i];
    Fynew[i] = Fy[i];
    Fznew[i] = Fz[i];
  }


  //calculate the radiation pressure tensor
  get_pressure_tensor(enzo_block, N, Fx, Fy, Fz, clight_code);
  
  double Nmin = enzo_config->method_m1_closure_min_photon_density / Nunit;

  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
        int i = INDEX(ix,iy,iz,mx,my); //index of current cell
        double N_update=0, Fx_update=0, Fy_update=0, Fz_update=0;
      
        get_U_update( enzo_block, &N_update, &Fx_update, &Fy_update, &Fz_update,
                             N, Fx, Fy, Fz, hx, hy, hz, dt, clight_code,
                             i, idx, idy, idz ); 
        
        // get updated fluxes
        Fxnew[i] += Fx_update;
        Fynew[i] += Fy_update;
        Fznew[i] += Fz_update;

        // now get updated photon densities
        Nnew[i] = std::max(Nnew[i] + N_update, Nmin);


      #ifdef DEBUG_TRANSPORT
        CkPrintf("i = %d; N_update = %f; Fx_update = %f; Nnew[i] = %f; hx = %f; dt = %f \n", i, N_update, Fx_update, Nnew[i], hx, dt);
      #endif


        // add interactions with matter 

        double C = 0.0; // photon creation term
        double D = 0.0; // photon destruction term

        if (enzo_config->method_m1_closure_attenuation) {
          D = D_add_attenuation(enzo_block, clight_cgs, i, igroup);
        }

        if (enzo_config->method_m1_closure_recombination_radiation) {
          // update photon density due to recombinations
          // Grackle does recombination chemistry, but doesn't
          // do anything about the radiation that comes out of recombination
          C = C_add_recombination(enzo_block, T, i, igroup, E_lower, E_upper);
        }
        
        // update radiation fields due to thermochemistry (see appendix A)
        double mult = 1.0/(1+dt*D);
        Nnew [i] = std::max((Nnew [i] + dt*C) * mult, Nmin);
        Fxnew[i] = Fxnew[i] * mult;
        Fynew[i] = Fynew[i] * mult;
        Fznew[i] = Fznew[i] * mult;
      }
    } 
  } 
  
  // now copy values over  
  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
        int i = INDEX(ix,iy,iz,mx,my); //index of current cell
        N [i] = Nnew [i];
        Fx[i] = Fxnew[i];
        Fy[i] = Fynew[i];
        Fz[i] = Fznew[i];
  
        if ( isnan(N[i]) ) {
          ERROR("EnzoMethodM1Closure::solve_transport_eqn()", 
                "N[i] is NaN!\n");
        }
      }
    }
  } 

  delete [] Nnew;
  delete [] Fxnew;
  delete [] Fynew;
  delete [] Fznew;
}

//----------------------------------------------------------------------

void EnzoMethodM1Closure::add_LWB(EnzoBlock * enzo_block, double J21) 
{

  // Adds cosmological Lyman-Werner background. If J21 > 0, LWB intensity will be 
  // constant with a value specified by J21. 
  // Otherwise, J21 will be calculated using the redshift-dependent polynomial fit from 
  // Wise et al. (2012)
  // The LW radiation field is assumed to be group 0
  const EnzoConfig * enzo_config = enzo::config();
  double energy = enzo_config->method_m1_closure_energy_mean[0];
  ASSERT("EnzoMEthodM1Closure::add_LWB", 
         "Adding LWB, but photon group 0 is not in the lyman-werner band",
         (11.18 < energy) && (energy < 13.60));
 
  Field field = enzo_block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones

  double hx, hy, hz;
  enzo_block->cell_width(&hx,&hy,&hz);

  EnzoUnits * enzo_units = enzo::units();
  double Nunit = enzo_units->photon_number_density();
    
  double JLW;
  if (J21 >= 0.0) {
    JLW = J21 * 1e-21;
  }

  else {
    ASSERT("EnzoMethodM1Closure::add_LWB", 
           "Must be running a cosmology simulation to use a time-dependent Lyman-Werner background", 
           enzo::cosmology());

    double A = -2.567;
    double B = 0.4562;
    double C = -0.02680;
    double D = 5.882e-4;
    double E = -5.056e-6;

    double z = enzo_block->state()->redshift();
    if (z > 30) { // function valid for z < 30 (not many Pop III stars at z > 30)
      return;
    }
    double log_J21 = A + B*z + C*z*z + D*z*z*z + E*z*z*z*z;
    JLW = pow(10, log_J21) * 1e-21; // erg s^-1 cm^-2 Hz^-1 sr^-1
  }

  enzo_float * N = (enzo_float *) field.values("photon_density_0");

  double energy_cgs = energy * enzo_constants::erg_eV;
  double dnu = (enzo_config->method_m1_closure_energy_upper[0] - 
                enzo_config->method_m1_closure_energy_lower[0]) * 
                enzo_constants::erg_eV / enzo_constants::hplanck; // frequency in Hz

  double Nbackground = 4*cello::pi * JLW/(energy_cgs*enzo_constants::clight) * dnu / Nunit;

  Scalar<double> scalar = enzo_block->data()->scalar_double();
  
  bool is_first_cycle =
    (enzo_block->state()->cycle() == enzo_config->initial_cycle);
  double Nbackground_previous =
    is_first_cycle ? 0.0 : *(scalar.value(scalar.index("N_LWB")));
 
  for (int i=0; i<mx*my*mz; i++) {
    // subtract the background from the previous timestep and add the new
    N[i] += Nbackground - Nbackground_previous;
  }

  *(scalar.value( scalar.index("N_LWB") )) = Nbackground;

}

//----------------------------------------------------------------------

void EnzoMethodM1Closure::call_inject_photons(EnzoBlock * enzo_block) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  
  const int N_groups = enzo_config->method_m1_closure_N_groups;
  const int N_species = 3 + enzo_config->method_m1_closure_H2_photodissociation; 
  
  if (enzo_block->is_leaf()) { // only inject photons for leaf blocks

    if (enzo_config->method_m1_closure_lyman_werner_background) {
      add_LWB(enzo_block, enzo_config->method_m1_closure_LWB_J21);
    }

    for (int igroup=0; igroup<N_groups; igroup++) {
      this->inject_photons(enzo_block, igroup);
    }
  }
 
  // set group mean cross sections and energies
  //
  // "vernier_average" -- calculates cross section from sigma_vernier() function,
  //    then averages that value over all star particles in the simulation, weighted by 
  //    mass * luminosity
  //  "vernier" -- just sets cross sections equal to values from sigma_vernier()
  //  "custom" -- sets cross sections to user-specified values in the parameter file
  Scalar<double> scalar = enzo_block->data()->scalar_double();

  for (int i=0; i<N_groups; i++) {
    *(scalar.value( scalar.index( eps_string(i) ))) = 0.0;
    for (int j=0; j<N_species; j++) 
    {
     *(scalar.value( scalar.index( sigN_string(i,j) ))) = 0.0;
     *(scalar.value( scalar.index( sigE_string(i,j) ))) = 0.0;
    }
  }  
    
  if (enzo_config->method_m1_closure_cross_section_calculator == "vernier_average") {
    // do global reduction of sigE, sigN, and eps over star particles
    // then do refresh -> solve_transport_eqn()

    // flattened array of sigN, sigE, and eps "mL" variables
    // N_groups number of mL_i, eps_i variables
    // N_groups*N_species number of sigN and sigE variables
    //
    // Gives 2*N_groups+2*N_groups*N_species as the total number of variables

    std::vector<double> temp(2*N_groups+2*N_groups*N_species);

    // fill temp with ScalarData "mL" quantities
    int index = 0;
    for (int i=0; i<N_groups; i++) {
      index = i;
      temp[index] = *(scalar.value( scalar.index(mL_string(i)) ));
    }
    for (int i=0; i<N_groups; i++) {
      index = N_groups + i;
      temp[index] = *(scalar.value( scalar.index( eps_string(i) + mL_string(i) )));
    }
    for (int i=0; i<N_groups; i++) {
      for (int j=0; j<N_species; j++) {
        index = 2*N_groups + i*N_species + j;
        temp[index] = *(scalar.value( scalar.index( sigN_string(i,j) + mL_string(i) )));
      }
    }
    for (int i=0; i<N_groups; i++) {
      for (int j=0; j<N_species; j++) {
        index = 2*N_groups + N_groups*N_species + i*N_species + j;
        temp[index] = *(scalar.value( scalar.index( sigE_string(i,j) + mL_string(i) )));
      }
    }

    CkCallback callback (CkIndex_EnzoBlock::p_method_m1_closure_set_global_averages(NULL),
             enzo_block->proxy_array());

    enzo_block->contribute(temp, CkReduction::sum_double, callback);
  } else { // just set sigmaN = sigmaE = either sigma_vernier or custom value, and eps = mean(energy)

    if (! enzo_block->is_leaf() ) {
      enzo_block->compute_done();
      return;
    }

    for (int i=0; i<N_groups; i++) {
      double E_lower = (enzo_config->method_m1_closure_energy_lower)[i];
      double E_upper = (enzo_config->method_m1_closure_energy_upper)[i];
      double energy = (enzo_config->method_m1_closure_energy_mean)[i]; // eV
      *(scalar.value( scalar.index( eps_string(i) ))) = energy*enzo_constants::erg_eV; // erg
      
      if (enzo_config->method_m1_closure_cross_section_calculator == "vernier") {
        // set sigmaN = sigmaE = sigma_vernier
        for (int j=0; j<N_species; j++) {
          double sigma_j = sigma_vernier(energy,j); // cm^2
          *(scalar.value( scalar.index( sigN_string(i,j) ))) = sigma_j;
          *(scalar.value( scalar.index( sigE_string(i,j) ))) = sigma_j;
        }
      } else if (enzo_config->method_m1_closure_cross_section_calculator == "custom") {
        // set sigmaN = sigmaE = custom values
        for (int j=0; j<N_species; j++) {
          int sig_index = i*N_species + j;
          double sigmaN_ij = enzo_config->method_m1_closure_sigmaN[sig_index]; // cm^2
          double sigmaE_ij = enzo_config->method_m1_closure_sigmaE[sig_index];
          *(scalar.value( scalar.index( sigN_string(i,j) ))) = sigmaN_ij;
          *(scalar.value( scalar.index( sigE_string(i,j) ))) = sigmaE_ij;
        }
      }
 
    }

    cello::refresh(ir_injection_)->set_active(enzo_block->is_leaf()); 
    enzo_block->refresh_start(ir_injection_, CkIndex_EnzoBlock::p_method_m1_closure_solve_transport_eqn());   
  }
  
}

//-----------------------------------

void EnzoBlock::p_method_m1_closure_set_global_averages(CkReductionMsg * msg)
{
  EnzoMethodM1Closure * method = static_cast<EnzoMethodM1Closure*> (this->method()); 
  method->set_global_averages(this, msg);
}

//----

void EnzoMethodM1Closure::set_global_averages(EnzoBlock * enzo_block, CkReductionMsg * msg) throw()
{
  // contribute does global reduction over ALL blocks by default (not just leaves)
  // call compute_done here for non leaves so that we don't waste time 
  // pushing these blocks through solve_transport_eqn().

  const EnzoConfig * enzo_config = enzo::config();
  const int N_groups = enzo_config->method_m1_closure_N_groups;
  const int N_species = 3 + enzo_config->method_m1_closure_H2_photodissociation;
  
  if (! enzo_block->is_leaf()) {
    enzo_block->compute_done(); 
    return;  
  }

  double * temp = (double *)msg->getData(); // pointer to vector containing numerator/denominators
                                            // of eqs. (B6)-(B8)
                                            // temp[0:N_groups] hold the denominators -> sum(m*L_i)
                                            // temp[N_groups:] hold the numerators -> sum(<eps/sigN/sigE>_ij * m*L_i)

  Scalar<double> scalar = enzo_block->data()->scalar_double();

  int index=0;
  double mult = 0.0;
  for (int i=0; i<N_groups; i++) {
    mult = (temp[i] == 0) ? 0.0 : 1.0/temp[i];
    index = N_groups + i;
    // eq. B6 --> sum(eps*m*L) / sum(m*L)
    *(scalar.value( scalar.index( eps_string(i) ))) = mult*temp[index];
  }

  for (int i=0; i<N_groups; i++) {
    for (int j=0; j<N_species; j++) {
      mult = (temp[i] == 0) ? 0.0 : 1.0/temp[i];
      index = 2*N_groups + i*N_species + j;
      // eq. B7 --> sum(sigN*m*L) / sum(m*L)
      *(scalar.value( scalar.index( sigN_string(i,j) ))) = mult*temp[index];
    }
  }

  for (int i=0; i<N_groups; i++) {
    for (int j=0; j<N_species; j++) {
      mult = (temp[i] == 0) ? 0.0 : 1.0/temp[i];
      index = 2*N_groups + N_groups*N_species + i*N_species + j;
      // eq. B8 --> sum(sigE*m*L) / sum(m*L)
      *(scalar.value( scalar.index( sigE_string(i,j) ))) = mult*temp[index];
    }
  }
  
  cello::refresh(ir_injection_)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_start(ir_injection_, CkIndex_EnzoBlock::p_method_m1_closure_solve_transport_eqn());
}

//-----------------------------------

void EnzoMethodM1Closure::call_solve_transport_eqn(EnzoBlock * enzo_block) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  int N_groups = enzo_config->method_m1_closure_N_groups;
  double clight = enzo_config->method_m1_closure_clight_frac * enzo_constants::clight;

  // loop through groups and solve transport equation for each group
  for (int igroup=0; igroup<N_groups; igroup++) {
    this->solve_transport_eqn(enzo_block, igroup);
  }

  if (enzo_config->method_m1_closure_thermochemistry) {
    // Calculate photoheating and photoionization rates.
    // Sums over frequency groups
    get_photoionization_and_heating_rates(enzo_block, clight);
  }

}

//-------------------------------

void EnzoBlock::p_method_m1_closure_solve_transport_eqn()
{
  EnzoMethodM1Closure * method = static_cast<EnzoMethodM1Closure*> (this->method());

  // solve transport eqn for each group
  method->call_solve_transport_eqn(this);

  // sum group fields
  method->sum_group_fields(this);  
  compute_done();
  return; 
}

//-----------------------------------------

void EnzoMethodM1Closure::sum_group_fields(EnzoBlock * enzo_block) throw()
{
  // this function does two things sums group fields together and stores them in integrated fields
  Field field = enzo_block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);

  const int m = mx*my*mz;
  const EnzoConfig * enzo_config = enzo::config();
        EnzoUnits  * enzo_units  = enzo::units();
  
  enzo_float * N  = (enzo_float *) field.values("photon_density");
  enzo_float * Fx = (enzo_float *) field.values("flux_x");
  enzo_float * Fy = (enzo_float *) field.values("flux_y");
  enzo_float * Fz = (enzo_float *) field.values("flux_z");

  for (int j=0; j<m; j++)
  {
    N [j] = 0.0;
    Fx[j] = 0.0;
    Fy[j] = 0.0;
    Fz[j] = 0.0; 
  }

  int N_groups = enzo_config->method_m1_closure_N_groups;
  for (int i=0; i < N_groups; i++) {
    std::string istring = std::to_string(i);
    enzo_float *  N_i = (enzo_float *) field.values("photon_density_" + istring);
    enzo_float * Fx_i = (enzo_float *) field.values("flux_x_" + istring);
    enzo_float * Fy_i = (enzo_float *) field.values("flux_y_" + istring);
    enzo_float * Fz_i = (enzo_float *) field.values("flux_z_" + istring);
    for (int j=0; j<m; j++)
    {
      N [j] +=  N_i[j];
      Fx[j] += Fx_i[j];
      Fy[j] += Fy_i[j];
      Fz[j] += Fz_i[j];
    }    
  }
}
