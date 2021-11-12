// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRamsesRT.cpp
/// @author   William Hicks (whicks@ucsd.edu)
/// @date     Mon Aug 16 17:05:23 PDT 2021
/// @brief    Implements the EnzoMethodRamsesRT  class

#include "cello.hpp"

#include "enzo.hpp"

//#define DEBUG_PRINT_GROUP_PARAMETERS
//----------------------------------------------------------------------

EnzoMethodRamsesRT ::EnzoMethodRamsesRT(const int N_groups, const double clight)
  : Method()
    , N_groups_(N_groups)
    , clight_(clight)
    , eps_()
    , sigN_()
    , sigE_()
    , gfracN_()
    , gfracF_()
    , ir_injection_(-1)
    , ir_transport_(-1)
    //, igroup_()
{

  const int rank = cello::rank();

  this->required_fields_ = std::vector<std::string> {"photon_density"}; //number density

  if (rank >= 1) {
    this->required_fields_.push_back("flux_x");
    this->required_fields_.push_back("P00"); // elements of the pressure tensor
  }                                         
  if (rank >= 2) {
    this->required_fields_.push_back("flux_y");
    this->required_fields_.push_back("P10");
    this->required_fields_.push_back("P01");
    this->required_fields_.push_back("P11");
  }
  if (rank >= 3) {
    this->required_fields_.push_back("flux_z");
    this->required_fields_.push_back("P02");
    this->required_fields_.push_back("P12");
    this->required_fields_.push_back("P20");
    this->required_fields_.push_back("P21");
    this->required_fields_.push_back("P22");
  }


  for (int i=0; i<N_groups_; i++) {
    std::string istring = std::to_string(i); 
    this->required_fields_.push_back("photon_density_" + istring);
    if (rank >= 1) this->required_fields_.push_back("flux_x_" + istring);
    if (rank >= 2) this->required_fields_.push_back("flux_y_" + istring);
    if (rank >= 3) this->required_fields_.push_back("flux_z_" + istring);   
  }

  // define required fields if they do not exist
  this->define_fields();

  // Initialize default Refresh object

  cello::simulation()->new_refresh_set_name(ir_post_,name());
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

  // Initialize Refresh object for after injection step 
  ir_injection_ = add_new_refresh_();

  cello::simulation()->new_refresh_set_name(ir_injection_, name()+":injection");
  Refresh * refresh_injection = cello::refresh(ir_injection_);
  refresh_injection->add_field("photon_density");

  if (rank >= 1) refresh_injection->add_field("flux_x");
  if (rank >= 2) refresh_injection->add_field("flux_y");
  if (rank >= 3) refresh_injection->add_field("flux_z");
  
  for (int i=0; i<N_groups_; i++) {
    std::string istring = std::to_string(i); 
    refresh_injection->add_field("photon_density_" + istring);
    if (rank >= 1) refresh_injection->add_field("flux_x_" + istring);
    if (rank >= 2) refresh_injection->add_field("flux_y_" + istring);
    if (rank >= 3) refresh_injection->add_field("flux_z_" + istring);   
  }
 
  refresh_injection->set_callback(CkIndex_EnzoBlock::p_method_ramses_rt_solve_transport_eqn()); 
  
  // store frequency group attributes as ScalarData variables
  ScalarDescr * scalar_descr = cello::scalar_descr_double();
  
  int N_species_ = 3; //only three ionizable species (HI, HeI, HeII)
  for (int i=0; i<N_groups_; i++) {
    std::string istring = std::to_string(i);
    eps_ .push_back( scalar_descr->new_value("eps_" +istring ));
    for (int j=0; j<N_species_; j++) {
      std::string jstring = std::to_string(j);
      sigN_.push_back( scalar_descr->new_value("sigN_"+istring+jstring ));
      sigE_.push_back( scalar_descr->new_value("sigE_"+istring+jstring ));
    }
  }

}

//----------------------------------------------------------------------

void EnzoMethodRamsesRT ::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | N_groups_;
  p | clight_;
  p | eps_;
  p | sigN_;
  p | sigE_;
  p | gfracN_;
  p | gfracF_;
  p | ir_injection_;
  p | ir_transport_;
  //p | igroup_;
}

//----------------------------------------------------------------------

void EnzoMethodRamsesRT::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    compute_ (block);
  }

  //block->compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodRamsesRT::timestep ( Block * block ) const throw()
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
  //EnzoUnits * enzo_units = enzo::units();

  //enzo_block = enzo::block(block);
  return h_min / (3*enzo_config->method_ramses_rt_clight);//clight_);
}

//-----------------------------------------------------------------------

double EnzoMethodRamsesRT::integrate_simpson(double a, double b,
                    int n, // Number of intervals
                    std::function<double(double,double,double,int)> f, double v1, double v2, int v3) throw()
{
    //TODO (?): The elegant thing to do here would be to somehow let the parameter f be a 
    //          variadic function, and take va_list as a parameter instead of (v1,v2,v3,...). 
    //          If I'm only using this to integrate planck functions though, I probably don't need to 
    //          waste time making this as general as possible
              
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

double EnzoMethodRamsesRT::planck_function(double nu, double T, double clight, int dependent_variable) throw()
{
  double prefactor = 0.0;
  switch (dependent_variable)
  {
    case 0: //photon density
       prefactor = 8*cello::pi*nu*nu / (clight*clight*clight); 
       break;
    case 1: //energy density
       prefactor = 8*cello::pi*cello::hplanck*nu*nu*nu / (clight*clight*clight); 
       break;
  }
  
  return prefactor / ( pow(cello::e, cello::hplanck*nu/(cello::kboltz*T)) -1 );

}

//-------------------------- INJECTION STEP ----------------------------
double EnzoMethodRamsesRT::get_star_temperature(double M) throw()
{
  double L, R;
  // mass-luminosity relations for main sequence stars
  if      (M/cello::mass_solar < 0.43) L = 0.23*pow(M/cello::mass_solar,2.3);
  else if (M/cello::mass_solar < 2   ) L =      pow(M/cello::mass_solar,4.0);
  else if (M/cello::mass_solar < 55  ) L =  1.4*pow(M/cello::mass_solar,3.5);
  else L = 32000*(M/cello::mass_solar);

  // mass-radius relations (need to find more accurate version for large masses?)
  if (M/cello::mass_solar < 1) R = pow(M/cello::mass_solar, 0.8);
  else R = pow(M/cello::mass_solar, 0.57);
  
  L *= cello::luminosity_solar;
  R *= cello::radius_solar;

  return pow( L/(4*cello::pi*R*R*cello::sigma_SF), 0.25);
}

void EnzoMethodRamsesRT::get_radiation_blackbody(EnzoBlock * enzo_block, enzo_float * N, int i, double T, 
                   double freq_lower, double freq_upper, double dt, double clight, double f_esc) throw()
{
  int n = 10; // number of partitions for simpson's method

  //need to use lambda expression to pass in planck_function() function as a parameter
  //because member function in c++ are automatically attached to the `this` pointer,
  //so you need to capture `this` in order for the compiler to know which function
  //to point to
  
  if (freq_lower == 0.0) freq_lower = 1; //planck function undefined at zero
                                         //1 Hz is a very small frequency compared to ~1e16 Hz
                                         
  double N_integrated = integrate_simpson(freq_lower,freq_upper,n, 
                 [this](double a, double b, double c, int d) {return planck_function(a,b,c,d);},
                 T,clight,0); 
  double E_integrated = integrate_simpson(freq_lower,freq_upper,n, 
                 [this](double a, double b, double c, int d) {return planck_function(a,b,c,d);},
                 T,clight,1);
 
  N[i] = f_esc * N_integrated; 
 
  // "integrate" by just drawing a rectangle through the midpoint
  //double freq_center = freq_lower + 0.5 * (freq_upper - freq_lower);
  //N[i] = f_esc*(freq_upper-freq_lower) * 8*cello::pi*freq_center*freq_center / (clight*clight*clight) 
  //      / ( pow(cello::e, cello::hplanck*freq_center/(cello::kboltz*T)) -1 );

   
   //for testing
   //N[i] = 1e16;
   std::vector<std::string> chemistry_fields = {"HI_density", 
                                               "HeI_density", "HeII_density"};
   std::vector<double> masses = {cello::mass_hydrogen,
                      4*cello::mass_hydrogen, 4*cello::mass_hydrogen};

   //Calculate photon group attributes----
#ifdef DEBUG_PRINT_GROUP_PARAMETERS
   if (enzo_block->method_ramses_rt_igroup > 0) std::cout << "~~~~~~" << std::endl;
#endif
   Scalar<double> scalar = enzo_block->data()->scalar_double(); 
 
   //eq. B3 ----> int(E_nu dnu) / int(N_nu dnu)
   *(scalar.value( scalar.index("eps_" + 
                  std::to_string(enzo_block->method_ramses_rt_igroup)) ))
                              =
                 E_integrated / N_integrated; 

   for (int j=0; j<chemistry_fields.size(); j++) {

     //TODO: Need to average sigma_ij over all stars in the simulation. See eq. B6-B8
  
     // eq. B4 ----> int(sigma_nuj * N_nu dnu)/int(N_nu dnu)
     *(scalar.value( scalar.index("sigN_" + std::to_string(enzo_block->method_ramses_rt_igroup) 
                               + std::to_string(j)) )) 
                                     =
            integrate_simpson(freq_lower,freq_upper,n, 
                 [this,j](double nu, double b, double c, int d){
                   return sigma_vernier(cello::hplanck*nu / cello::erg_eV, j)*planck_function(nu,b,c,d);
                 },
                 T,clight,0) / N_integrated;

     // eq. B5 ----> int(sigma_nuj * E_nu dnu)/int(E_nu dnu)
     *(scalar.value( scalar.index("sigE_" + std::to_string(enzo_block->method_ramses_rt_igroup) 
                              + std::to_string(j)) ))
                                     =
            integrate_simpson(freq_lower,freq_upper,n, 
                 [this,j](double nu, double b, double c, int d){
                   return sigma_vernier(cello::hplanck*nu / cello::erg_eV, j)*planck_function(nu,b,c,d);
                 },
                 T,clight,1) / E_integrated;
#ifdef DEBUG_PRINT_GROUP_PARAMETERS

     const EnzoConfig * enzo_config = enzo::config();
     
     std::cout<<'['  << enzo_config->method_ramses_rt_bin_lower[enzo_block->method_ramses_rt_igroup] 
              << ", " << enzo_config->method_ramses_rt_bin_upper[enzo_block->method_ramses_rt_igroup]
              << "] " << j
              << ' '  << *(scalar.value( scalar.index("sigN_" + 
                          std::to_string(enzo_block->method_ramses_rt_igroup) + std::to_string(j)) ))
              << ' '  << *(scalar.value( scalar.index("sigE_" + 
                          std::to_string(enzo_block->method_ramses_rt_igroup) + std::to_string(j)) )) 
              << ' '  << *(scalar.value( scalar.index("eps_" + 
                          std::to_string(enzo_block->method_ramses_rt_igroup) ))) / cello::erg_eV
              << std::endl;  
#endif         


   // for testing
   //*(scalar.value( scalar.index("sigN_" + std::to_string(enzo_block->method_ramses_rt_igroup) 
   //                           + std::to_string(j)) )) = enzo_block->method_ramses_rt_igroup * 2e-14; 

  }
  //----------------------------------
  

}

void EnzoMethodRamsesRT::inject_photons ( EnzoBlock * enzo_block ) throw()
{
  // Solve dN_i/dt = Ndot^*_i
  // routine for identifying star particles and getting their grid position
  // copy/pasted from enzo_EnzoMethodFeedback::compute_()

  const EnzoConfig * enzo_config = enzo::config();
  double f_esc = 1.0; // enzo_config->method_radiation_injection_escape_fraction;

  Field field = enzo_block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);
  
  int idx = 1;
  int idy = mx;
  int idz = mx*my;
 
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);
 
  double xm,ym,zm;
  double xp,yp,zp;
  double hx,hy,hz;
  enzo_block->data()->lower(&xm,&ym,&zm);
  enzo_block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  // if (cosmology) {...}
 
  double inv_vol = 1.0/(hx*hy*hz); 

  double dt = enzo_block->dt;

  // get relevant field variables
  enzo_float * N  = (enzo_float *) field.values(
        "photon_density_"+std::to_string(enzo_block->method_ramses_rt_igroup));
 
  std::vector<std::string> chemistry_fields = {"HI_density", 
                                               "HeI_density", "HeII_density"};
  std::vector<double> masses = {cello::mass_hydrogen,
                      4*cello::mass_hydrogen, 4*cello::mass_hydrogen};


  /*
  enzo_float * Fx  = (enzo_float *) field.values(
        "flux_x_"+std::to_string(enzo_block->method_ramses_rt_igroup));
  enzo_float * Fy  = (enzo_float *) field.values(
        "flux_y_"+std::to_string(enzo_block->method_ramses_rt_igroup));
  enzo_float * Fz  = (enzo_float *) field.values(
        "flux_z_"+std::to_string(enzo_block->method_ramses_rt_igroup));
  */

  Particle particle = enzo_block->data()->particle();
  int it = particle.type_index("star");
  
  // if no stars, don't do anything
  if (particle.num_particles(it) == 0) return;
  
  const int ia_m = particle.attribute_index (it, "mass");

  const int ia_x = (rank >= 1) ? particle.attribute_index (it, "x") : -1;
  const int ia_y = (rank >= 2) ? particle.attribute_index (it, "y") : -1;
  const int ia_z = (rank >= 3) ? particle.attribute_index (it, "z") : -1;

  const int dm = particle.stride(it, ia_m);
  const int dp = particle.stride(it, ia_x);

  const int nb = particle.num_batches(it);

  double freq_lower = (enzo_config->method_ramses_rt_bin_lower)[enzo_block->method_ramses_rt_igroup] * cello::erg_eV / cello::hplanck;
  double freq_upper = (enzo_config->method_ramses_rt_bin_upper)[enzo_block->method_ramses_rt_igroup] * cello::erg_eV / cello::hplanck;

  double clight = enzo_config->method_ramses_rt_clight;   

  for (int ib=0; ib<nb; ib++){
    enzo_float *px=0, *py=0, *pz=0;
    enzo_float *pmass=0;
   
    pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);

    px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

    int np = particle.num_particles(it,ib);

    // loop through particles within each batch
    for (int ip=0; ip<np; ip++) {
      int ipdp = ip*dp;
      int ipdm = ip*dm;

      //TODO: Add check for whether star particle is still active (has it blown up yet?)
        
      // get corresponding grid position
      double x_part = (px[ipdp] - xm) / hx;
      double y_part = (py[ipdp] - ym) / hy;
      double z_part = (pz[ipdp] - zm) / hz;

      // get 3D grid index for particle - account for ghost zones!!
      int ix = ((int) std::floor(x_part))  + gx;
      int iy = ((int) std::floor(y_part))  + gy;
      int iz = ((int) std::floor(z_part))  + gz;

      // now get index of this cell
      int i = INDEX(ix,iy,iz,mx,my);
      
      // deposit photons
       
      // This function will take in particle mass as a parameter, and will fit
      // and integrate over the SED to get the total injection rate. The current
      // implementation in RAMSES assumes that SED of gas stays constant.
      // How do I keep track of the gas SED if different mass stars radiate differently?
      
      // Get temperature of star 

      double T = get_star_temperature(pmass[ipdm]);

      get_radiation_blackbody(enzo_block, N, i, T, freq_lower, freq_upper, dt, clight,f_esc);
      //TODO: Only call get_cross_section here every N cycles, where N is an input parameter.
      //      Also only call if igroup=0 so that it only gets called once
     

      /*
      // initialize fluxes within a radius of one cell
      for (int k_=iz-1; k_<=iz+1;k_++) {
        for (int j_=iy-1; j_<=iy+1;j_++) {
          for (int i_=ix-1; i_<=ix+1;i_++) {
            int index = INDEX(i_,j_,k_,mx,my);
            if (index == i) continue;
            double distsqr = (i_-ix)*(i_-ix)*hx*hx + (j_-iy)*(i_-iy)*hy*hy + (k_-iz)*(k_-iz)*hz*hz;
                        
            N[index]  = N[i]/27;//distsqr;
            
            Fx[index] = (i_-ix); //i_-ix is either +- 1, so this gives direction of flux
            Fy[index] = (j_-iy);
            Fz[index] = (k_-iz);
          
            // normalize fluxes 
            double Ftot = Fx[index]*Fx[index] + Fy[index]*Fy[index] + Fz[index]*Fz[index];                      
            Fx[index] *= clight*N[index]/Ftot; 
            Fy[index] *= clight*N[index]/Ftot;
            Fz[index] *= clight*N[index]/Ftot; 
         }
        }
      }
      */
//      N[i] += f_esc * injection_rate * dt * inv_vol; // photons/s * timestep / volume

 
      // I don't have to directly alter the fluxes here because that naturally gets
      // taken care of during the transport step 
      
    }
  } 
}

//---------------------------------------------------------------------

//----------------------------------------------------------------------

double EnzoMethodRamsesRT::flux_function (double U_l, double U_lplus1,
					   double Q_l, double Q_lplus1,  
					   std::string type) 
					   throw()
{
  // returns face-flux of a cell at index idx
  const EnzoConfig * enzo_config = enzo::config();
  if (type == "GLF") {
    return 0.5*(  Q_l+Q_lplus1 - enzo_config->method_ramses_rt_clight*(U_lplus1-U_l) ); 
  }

  else {
    ERROR("EnzoMethodRamsesRT::flux_function",
	     "flux_function type not recognized");
    return 0.0; 
  }
}

double EnzoMethodRamsesRT::deltaQ_faces (double U_l, double U_lplus1, double U_lminus1,
                                                  double Q_l, double Q_lplus1, double Q_lminus1)
                                                  throw()
{
  // calls flux_function(), and calculates Q_{i-1/2} - Q_{i+1/2}
  
  return flux_function(U_lminus1, U_l, Q_lminus1, Q_l, "GLF") - flux_function(U_l, U_lplus1, Q_l, Q_lplus1, "GLF"); 
}



void EnzoMethodRamsesRT::get_reduced_variables (double * chi, double (*n)[3], int i, double clight,
                               enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz) throw()
{
        double Fnorm = sqrt(Fx[i]*Fx[i] + Fy[i]*Fy[i] + Fz[i]*Fz[i]);
        double f = Fnorm / (clight*N[i]); // reduced flux
        *chi = (3 + 4*f*f) / (5 + 2*sqrt(4-3*f*f));
        (*n)[0] = Fx[i]/Fnorm;
        (*n)[1] = Fy[i]/Fnorm; 
        (*n)[2] = Fz[i]/Fnorm;
}

void EnzoMethodRamsesRT::get_pressure_tensor (EnzoBlock * enzo_block, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz) 
                       throw()
{
  Field field = enzo_block->data()->field();
  const EnzoConfig * enzo_config = enzo::config();
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
  
  for (int iz=gz-1; iz<mz-gz+1; iz++) { 
   for (int iy=gy-1; iy<my-gy+1; iy++) {
    for (int ix=gx-1; ix<mx-gx+1; ix++) {
      int i = INDEX(ix,iy,iz,mx,my); //index of current cell
      double chi, n[3]; 
      get_reduced_variables( &chi, &n, i, enzo_config->method_ramses_rt_clight, 
                            N, Fx, Fy, Fz);
      P00[i] = N[i] * ( 0.5*(1.0-chi) + 0.5*(3.0*chi - 1)*n[0]*n[0] );
      P10[i] = N[i] * 0.5*(3.0*chi - 1)*n[1]*n[0];
      P01[i] = N[i] * 0.5*(3.0*chi - 1)*n[0]*n[1];
      P11[i] = N[i] * ( 0.5*(1.0-chi) + 0.5*(3.0*chi - 1)*n[1]*n[1] );
      P02[i] = N[i] * 0.5*(3.0*chi - 1)*n[0]*n[2];
      P12[i] = N[i] * 0.5*(3.0*chi - 1)*n[1]*n[2];
      P20[i] = N[i] * 0.5*(3.0*chi - 1)*n[2]*n[0];
      P21[i] = N[i] * 0.5*(3.0*chi - 1)*n[2]*n[1];
      P22[i] = N[i] * ( 0.5*(1.0-chi) + 0.5*(3.0*chi - 1)*n[2]*n[2] );
    }
   }
  }

}

void EnzoMethodRamsesRT::get_U_update (EnzoBlock * enzo_block, double * N_update, 
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

  // if rank >= 1
  *N_update += dt/hx * deltaQ_faces( N[i],  N[i+idx],  N[i-idx],
                                   Fx[i], Fx[i+idx], Fx[i-idx] );
    
  *Fx_update += dt/hx * deltaQ_faces(Fx[i], Fx[i+idx], Fx[i-idx],
                                   clight*clight*P00[i],
                                   clight*clight*P00[i+idx],
                                   clight*clight*P00[i-idx] );

  // if rank >= 2
  *N_update += dt/hy * deltaQ_faces( N[i],  N[i+idy],  N[i-idy],
                                   Fy[i], Fy[i+idy], Fy[i-idy] ); 
    
  *Fx_update += dt/hy * deltaQ_faces(Fx[i], Fx[i+idy], Fx[i-idy],
                                   clight*clight*P10[i],
                                   clight*clight*P10[i+idy],
                                   clight*clight*P10[i-idy] );

  *Fy_update += dt/hx * deltaQ_faces(Fy[i], Fy[i+idx], Fy[i-idx],
                                   clight*clight*P01[i],
                                   clight*clight*P01[i+idx],
                                   clight*clight*P01[i-idx] );

  *Fy_update += dt/hy * deltaQ_faces(Fy[i], Fy[i+idy], Fy[i-idy],
                                   clight*clight*P11[i],
                                   clight*clight*P11[i+idy],
                                   clight*clight*P11[i-idy] );


   // if rank >= 3
  *N_update += dt/hz * deltaQ_faces( N[i],  N[i+idz],  N[i-idz],
                                   Fz[i], Fz[i+idz], Fz[i-idz] );

  *Fx_update += dt/hz * deltaQ_faces(Fx[i], Fx[i+idz], Fx[i-idz],
                                   clight*clight*P20[i],
                                   clight*clight*P20[i+idz],
                                   clight*clight*P20[i-idz] );

  *Fy_update += dt/hz * deltaQ_faces(Fy[i], Fy[i+idz], Fy[i-idz],
                                   clight*clight*P21[i],
                                   clight*clight*P21[i+idz],
                                   clight*clight*P21[i-idz]);

  *Fz_update += dt/hx * deltaQ_faces(Fz[i], Fz[i+idx], Fz[i-idx],
                                   clight*clight*P02[i],
                                   clight*clight*P02[i+idx],
                                   clight*clight*P02[i-idx]);

  *Fz_update += dt/hy * deltaQ_faces(Fz[i], Fz[i+idy], Fz[i-idy],
                                   clight*clight*P12[i],
                                   clight*clight*P12[i+idy],
                                   clight*clight*P12[i-idy]);

  *Fz_update += dt/hz * deltaQ_faces(Fz[i], Fz[i+idz], Fz[i-idz],
                                   clight*clight*P22[i],
                                   clight*clight*P22[i+idz],
                                   clight*clight*P22[i-idz]);


}

//----------------------------------
double EnzoMethodRamsesRT::sigma_vernier (double energy, int type) throw()
{
  // copy FindCrossSection.C from Enzo
  // Uses fits from Vernier et al. (1996) to calculate photoionization cross-section 
  // between photons of energy E and gas of a given species. Then need to average this
  // value over  
 
  float sigma;
  float e_th, e_max, e0, sigma0, ya, P, yw, y0, y1;

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
  }

  float x, y, fy;
  x = energy/e0 - y0;
  y = sqrt(x*x + y1*y1);
  fy = ((x-1.0)*(x-1.0) + yw*yw) * pow(y, 0.5*P-5.5) * 
    pow((1.0 + sqrt(y/ya)), -P);

  sigma = sigma0 * fy * 1e-18;

  return sigma;

}

//---------------------------------

void EnzoMethodRamsesRT::add_attenuation ( EnzoBlock * enzo_block, enzo_float * N, 
                     enzo_float * Fx, enzo_float * Fy, enzo_float * Fz, double clight, int i) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  
  Field field = enzo_block->data()->field();
  int igroup = enzo_block->method_ramses_rt_igroup;

  if (!field.is_field("density")) return;
  
  //enzo_float * density = (enzo_float *) field.values("density");

  //TODO: Only call get_cross_section here every N cycles, where N is an input parameter.
  //      Also only call if igroup=0 so that it only gets called once
  std::vector<std::string> chemistry_fields = {"HI_density", 
                                               "HeI_density", "HeII_density"};
  std::vector<double> masses = {cello::mass_hydrogen,
                      4*cello::mass_hydrogen, 4*cello::mass_hydrogen};

  //----------------------------------
 
 

  //it's okay to use the same cross section for both attenuation (affects N) and 
  //radiation pressure (affects F) because F and N have approximately the 
  //same spectral shape.

  Scalar<double> scalar = enzo_block->data()->scalar_double();
  double d_dt = 0.0;
  for (int j=0; j<chemistry_fields.size(); j++) {  
    enzo_float * density_j = (enzo_float *) field.values(chemistry_fields[j]);     
    double sigN_ij = *(scalar.value( scalar.index("sigN_" + std::to_string(igroup) + std::to_string(j)) ));
    
    d_dt += density_j[i] / masses[j] * clight*sigN_ij;
  }
  

  N [i] -= d_dt*N [i] * enzo_block->dt;
  Fx[i] -= d_dt*Fx[i] * enzo_block->dt;
  Fy[i] -= d_dt*Fy[i] * enzo_block->dt;
  Fz[i] -= d_dt*Fz[i] * enzo_block->dt; 
  
}

//----------------------

void EnzoMethodRamsesRT::solve_transport_eqn ( EnzoBlock * enzo_block ) throw()
{

  // TODO: Adapt for 2D and 3D cases
  
  // Solve dU/dt + del[F(U)] = 0; F(U) = { (Fx,Fy,Fz), c^2 P }
  //                                U  = { N, (Fx,Fy,Fz) }
  // M1 closure: P_i = D_i * N_i, where D_i is the Eddington tensor for 
  // photon group i
 
  const EnzoConfig * enzo_config = enzo::config();

  Field field = enzo_block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);
   
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);
 
  double xm,ym,zm;
  double xp,yp,zp;
  enzo_block->lower(&xm,&ym,&zm);
  enzo_block->upper(&xp,&yp,&zp);

  // array incremenent (because 3D array of field values are flattened to 1D)
  const int idx = 1;
  const int idy = mx;
  const int idz = mx*my;  

  std::string istring = std::to_string(enzo_block->method_ramses_rt_igroup);
  enzo_float * N  = (enzo_float *) field.values("photon_density_" + istring);
  enzo_float * Fx = (enzo_float *) field.values("flux_x_" + istring);
  enzo_float * Fy = (enzo_float *) field.values("flux_y_" + istring);
  enzo_float * Fz = (enzo_float *) field.values("flux_z_" + istring);

  const int m = mx*my*mz;
  // extra copy of fields needed to store
  // the evolved values until the end
  enzo_float * Nnew  = new enzo_float[m]; 
  enzo_float * Fxnew = new enzo_float[m];
  enzo_float * Fynew = new enzo_float[m];
  enzo_float * Fznew = new enzo_float[m];
 
  for (int i=0; i<m; i++)
  {
    if (N [i] < 1e-16) N [i] = 1e-16;
 
    Nnew [i] = N [i];
    Fxnew[i] = Fx[i];
    Fynew[i] = Fy[i];
    Fznew[i] = Fz[i];
  }

  double dt = enzo_block->dt;
  double hx = (xp-xm)/(mx-2*gx);
  double hy = (yp-ym)/(my-2*gy);
  double hz = (zp-zm)/(mz-2*gz);

  // Adjust values for cosmological expansion here? 
  // if (cosmology) {
  // }

  // if 2D, gz = 0 by default (I think) and mz = 1, so the outermost loop
  // would just be for (int iz=0, iz<1, iz++)

  //calculate the pressure tensor
  get_pressure_tensor(enzo_block, N, Fx, Fy, Fz);

  //if (true) get_group_attributes(); // if (enzo_block->cycle % attribute_cycle_step == 0)

  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
        int i = INDEX(ix,iy,iz,mx,my); //index of current cell
        double N_update=0, Fx_update=0, Fy_update=0, Fz_update=0;
      
        get_U_update( enzo_block, &N_update, &Fx_update, &Fy_update, &Fz_update,
                             N, Fx, Fy, Fz, hx, hy, hz, dt, enzo_config->method_ramses_rt_clight,
                             i, idx, idy, idz ); 

        // get updated fluxes

        Fxnew[i] += Fx_update;
        Fynew[i] += Fy_update;
        Fznew[i] += Fz_update;
         
        // now get updated photon densities
        Nnew[i] += N_update;
 
        // add interactions with matter, ignoring chemistry for now. Chemistry due to RT can be handled by Grackle 
        add_attenuation(enzo_block, Nnew, Fxnew, Fynew, Fznew, enzo_config->method_ramses_rt_clight, i);        
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
      }
    }
  } 

  delete [] Nnew;
  delete [] Fxnew;
  delete [] Fynew;
  delete [] Fznew;
}


//----------------------------------------------------------------------
void EnzoMethodRamsesRT::call_inject_photons(EnzoBlock * enzo_block) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  enzo_block->method_ramses_rt_igroup = 0;

#ifdef DEBUG_PRINT_GROUP_PARAMETERS
  //for printing out table of group parameters
  std::cout<<"[E_lower, E_upper] (eV); j; sigN_ij; sigE_ij; eps_i" << std::endl;
  std::cout<<"---------------------------------------------------" << std::endl;  
#endif 
 
  for (int i=0; i<enzo_config->method_ramses_rt_N_groups; i++) {
    this->inject_photons(enzo_block);
    enzo_block->method_ramses_rt_igroup += 1;
  }
  // do refresh, start transport step once that's done
 
  cello::refresh(ir_injection_)->set_active(enzo_block->is_leaf()); 
  enzo_block->new_refresh_start(ir_injection_, CkIndex_EnzoBlock::p_method_ramses_rt_solve_transport_eqn());
}

void EnzoMethodRamsesRT::call_solve_transport_eqn(EnzoBlock * enzo_block) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  enzo_block->method_ramses_rt_igroup = 0;
  for (int i=0; i<enzo_config->method_ramses_rt_N_groups;i++) {
    this->solve_transport_eqn(enzo_block);
    enzo_block->method_ramses_rt_igroup += 1;
  }
}

void EnzoBlock::p_method_ramses_rt_solve_transport_eqn()
{
  EnzoMethodRamsesRT * method = static_cast<EnzoMethodRamsesRT*> (this->method());
  method->call_solve_transport_eqn(this);

  // sum group fields and end compute()
  // TODO: Make tracking integrated group fields optional?
  // Move compute_done() call to after chemistry step once that's started 
  method->sum_group_fields(this); 
  this->compute_done(); 
}

//-----------------------------------------
void EnzoMethodRamsesRT::sum_group_fields(EnzoBlock * enzo_block) throw()
{
  Field field = enzo_block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);

  const int m = mx*my*mz;
  const EnzoConfig * enzo_config = enzo::config();

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

  for (int i=0; i < enzo_config->method_ramses_rt_N_groups; i++) {
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


//======================================================================

void EnzoMethodRamsesRT::compute_ (Block * block) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
 
  // TODO: lump the rest of the radiation that doesn't fall into these bins into one "other" group
  // and evolve that one too.
  //
  //
  // QUESTION: How do I get N_nu from the SED?
  // 
  //
  // LOOK AT DAN REYNOLDS' BRANCH OF ENZO-BW TO SEE HOW THIS IS IMPLEMENTED
  // See ComputeRadiationIntegrals() function and EnforceRadiationBounds()?
  //
  //
  // To test that multigroup is working, sample from a blackbody spectrum and plot the intensity in each group. Should trace out blackbody.

  Field field = block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);

  const int m = mx*my*mz;

  if (block->cycle() == 0)
  {
    for (int i=0; i<enzo_config->method_ramses_rt_N_groups; i++) {
      std::string istring = std::to_string(i);
      enzo_float *  N_i = (enzo_float *) field.values("photon_density_" + istring);
      enzo_float * Fx_i = (enzo_float *) field.values("flux_x_" + istring);
      enzo_float * Fy_i = (enzo_float *) field.values("flux_y_" + istring);
      enzo_float * Fz_i = (enzo_float *) field.values("flux_z_" + istring);
      for (int j=0; j<m; j++)
      {
        // initialize fields to zero 
        // TODO: check whether value parameter is set to choose whether or not to do this,
        // since we could also initialize the fields with the value parameter. 
        N_i [j] = 1e-16;
        Fx_i[j] = 1e-16;
        Fy_i[j] = 1e-16;
        Fz_i[j] = 1e-16; 
      }    
    }
  }

  //start photon injection step
  //This function will start the transport step after a refresh
  EnzoBlock * enzo_block = enzo::block(block);

  this->call_inject_photons(enzo_block);
}
