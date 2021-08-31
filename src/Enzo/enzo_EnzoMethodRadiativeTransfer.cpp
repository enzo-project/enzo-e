// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRadiativeTransfer.cpp
/// @author   William Hicks (whicks@ucsd.edu)
/// @date     Mon Aug 16 17:05:23 PDT 2021
/// @brief    Implements the EnzoMethodRadiativeTransfer  class

#include "cello.hpp"

#include "enzo.hpp"

//#define DEBUG_CHECK_TIMESTEP
//#define DEBUG_RT_TRANSPORT_STEP
//#define DEBUG_PRINT_REDUCED_VARIABLES
//#define DEBUG_CHECK_NAN
//#define DEBUG_CHECK_NAN1
//----------------------------------------------------------------------

EnzoMethodRadiativeTransfer ::EnzoMethodRadiativeTransfer()
  : Method()
{

  const int rank = cello::rank();

  this->required_fields_ = std::vector<std::string> {"photon_density"}; //number density

  if (rank >= 1) this->required_fields_.push_back("flux_x");
  if (rank >= 2) this->required_fields_.push_back("flux_y");
  if (rank >= 3) this->required_fields_.push_back("flux_z");

  // define required fields if they do not exist
  this->define_fields();

  // Initialize default Refresh object

  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("photon_density");

  if (rank >= 1) refresh->add_field("flux_x");
  if (rank >= 2) refresh->add_field("flux_y");
  if (rank >= 3) refresh->add_field("flux_z");

}

//----------------------------------------------------------------------

void EnzoMethodRadiativeTransfer ::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
}

//----------------------------------------------------------------------

void EnzoMethodRadiativeTransfer::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    compute_ (block);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodRadiativeTransfer::timestep ( Block * block ) const throw()
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
  const double clight = enzo_config->method_radiative_transfer_clight;


  return h_min / (3*clight);
}


//----------------------------------------------------------------------
void EnzoMethodRadiativeTransfer::inject_photons ( Block * block ) throw()
{
  // Solve dN_i/dt = Ndot^*_i 
}

//----------------------------------------------------------------------

double EnzoMethodRadiativeTransfer::flux_function (double U_l, double U_lplus1,
					   double Q_l, double Q_lplus1,  
					   std::string type,
					   double clight) 
					   throw()
{
  // returns face-flux of a cell at index idx
  if (type == "GLF") {
    return 0.5*(  Q_l+Q_lplus1 - clight*(U_lplus1-U_l) ); 
  }

  else {
    ERROR("EnzoMethodRadiativeTransfer::flux_function",
	     "flux_function type not recognized");
    return 0.0; 
  }
}

double EnzoMethodRadiativeTransfer::deltaQ_faces (double U_l, double U_lplus1, double U_lminus1,
                                                  double Q_l, double Q_lplus1, double Q_lminus1,
                                                  double clight) throw()
{
  // calls flux_function(), and calculates Q_{i-1/2} - Q_{i+1/2}
  
  return flux_function(U_lminus1, U_l, Q_lminus1, Q_l, "GLF", clight) - flux_function(U_l, U_lplus1, Q_l, Q_lplus1, "GLF", clight); 
}



void EnzoMethodRadiativeTransfer::get_reduced_variables (double * chi, double (*n)[3], int i, double clight,
                                                         enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz) throw()
{
        double Fnorm = sqrt(Fx[i]*Fx[i] + Fy[i]*Fy[i] + Fz[i]*Fz[i]);
        double f = Fnorm / (clight*N[i]); // reduced flux
        *chi = (3 + 4*f*f) / (5 + 2*sqrt(4-3*f*f));

#ifdef DEBUG_CHECK_NAN
        // We need 0 < f < 1
        if (isnan(*chi)) std::cout << " f: " << f << " N[i]: " << N[i] << std::endl;   
#endif
        (*n)[0] = Fx[i]/Fnorm;
        (*n)[1] = Fy[i]/Fnorm; 
        (*n)[2] = Fz[i]/Fnorm;
}

void EnzoMethodRadiativeTransfer::update_fluxes_1D (double * N_update, double * F_update, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz, 
                       double h, double dt, double clight, int i, int increment, int axis) throw()
{
  // calculate evolved photon density N

  double N_l = N[i], N_lplus1 = N[i+increment], N_lminus1 = N[i-increment];
  double F_l, F_lplus1, F_lminus1;
  if      (axis == 0) F_l = Fx[i], F_lplus1 = Fx[i+increment], F_lminus1 = Fx[i-increment];
  else if (axis == 1) F_l = Fy[i], F_lplus1 = Fy[i+increment], F_lminus1 = Fy[i-increment];
  else                F_l = Fz[i], F_lplus1 = Fz[i+increment], F_lminus1 = Fz[i-increment];

  *N_update += dt/h * deltaQ_faces( N_l, N_lplus1, N_lminus1,
                                    F_l, F_lplus1, F_lminus1, 
                                    clight );
                                                        
  // Calculate column `k = axis` of pressure tensor and update flux along that axis

  double P_l_k[3], P_lplus1_k[3], P_lminus1_k[3]; // TODO: Account for < 3 dimensions 
  double chi_l, chi_lplus1, chi_lminus1; // need to keep these for diagonal terms. Free to overwrite n though
  double n[3];
  int k = axis;
  
  //construct column `axis` of pressure tensor for cell `i` and its neighbors along direction of `axis`
  for (int j=0; j<3; j++) {
    get_reduced_variables(&chi_l, &n, i, clight, N, Fx, Fy, Fz);
    P_l_k[j] = 0.5*(3*chi_l-1) * n[j]*n[k] * N_l;
   
    get_reduced_variables(&chi_lplus1, &n, i+increment, clight, N, Fx, Fy, Fz);
    P_lplus1_k[j] = 0.5*(3*chi_lplus1-1) * n[j]*n[k] * N_lplus1;

    get_reduced_variables(&chi_lminus1, &n, i-increment, clight, N, Fx, Fy, Fz);
    P_lminus1_k[j] = 0.5*(3*chi_lminus1-1) * n[j]*n[k] * N_lminus1;   
  }  

#ifdef DEBUG_PRINT_REDUCED_VARIABLES             
    std::cout << " i: "           << i
              << " chi_l: "       << chi_l
              << " chi_lplus1: "  << chi_lplus1
              << " chi_lminus1: " << chi_lminus1 << std::endl;
#endif
 
  // diagonals += (1-chi)/2 
  P_l_k[k]       += 0.5*(1-chi_l) * N_l;
  P_lplus1_k[k]  += 0.5*(1-chi_lplus1) * N_lplus1;
  P_lminus1_k[k] += 0.5*(1-chi_lminus1) * N_lminus1;

  // update F (absorb into previous loop?)
  // F^{n+1}_i = F^n_i + dt/dx *c^2*(P_k^n_{i+1/2}-P_k^n_{i-1/2}), F and P_k are vectors here
  for (int j=0; j<3; j++) {
    *F_update += dt/h * deltaQ_faces(F_l, F_lplus1, F_lminus1,
                                     clight*clight * P_l_k[j], 
                                     clight*clight * P_lplus1_k[j], 
                                     clight*clight * P_lminus1_k[j], clight);
  }
  
}

void EnzoMethodRadiativeTransfer::transport_photons ( Block * block, double clight ) throw()
{
  // TODO: Make Nnew, Fxnew, Fynew, Fznew input parameters so that we don't need
  // the loop at the end that copies the values over to the original fields 
  // (Like in EnzoMethodHeat::compute_)

  // TODO: Adapt for 2D and 3D cases
  
  // TODO: Store elements of pressure tensor as field variable and access to avoid unnecessary recalculation

  // Solve dU/dt + del[F(U)] = 0; F(U) = { (Fx,Fy,Fz), c^2 P }
  //                                U  = { N, (Fx,Fy,Fz) }
  // M1 closure: P_i = D_i * N_i, where D_i is the Eddington tensor for 
  // photon group i
 

  Field field = block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);
   
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);
 
  double xm,ym,zm;
  double xp,yp,zp;
  block->lower(&xm,&ym,&zm);
  block->upper(&xp,&yp,&zp);

  // array incremenent (because 3D array of field values are flattened to 1D)
  const int idx = 1;
  const int idy = mx;
  const int idz = mx*my;

  enzo_float * N  = (enzo_float *) field.values("photon_density");
  enzo_float * Fx = (enzo_float *) field.values("flux_x");
  enzo_float * Fy = (enzo_float *) field.values("flux_y");
  enzo_float * Fz = (enzo_float *) field.values("flux_z");

  //~~~~~~~~ Make copy of fields to be updated ~~~~~~~~~~~~~//
  const int m = mx*my*mz;
  enzo_float * Nnew  = new enzo_float[m];
  enzo_float * Fxnew = new enzo_float[m];
  enzo_float * Fynew = new enzo_float[m];
  enzo_float * Fznew = new enzo_float[m]; 

  for (int i=0; i<m; i++)
  {
    Nnew [i] = N [i];
    Fxnew[i] = Fx[i];
    Fynew[i] = Fy[i];
    Fznew[i] = Fz[i];
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  double dt = timestep(block);
  double hx = (xp-xm)/(mx-2*gx);
  double hy = (yp-ym)/(my-2*gy);
  double hz = (zp-zm)/(mz-2*gz);

  // need to get a copy of each flux field  and density fields so we can update appropriately
  

  // Adjust values for cosmological expansion here? 
  // if (cosmology) {
  // }

  // loop through frequency bins
  
  // loop through cells
 
  // if 2D, gz = 0 by default (I think) and mz = 1, so the outermost loop
  // would just be for (int iz=0, iz<1, iz++)
  //

  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<my-gx; ix++) {
        int i = INDEX(ix,iy,iz,mx,my); //index of current cell
    
        double N_update=0, Fx_update=0, Fy_update=0, Fz_update=0;
        
        // if rank >= 1
        update_fluxes_1D( &N_update, &Fx_update, N, Fx, Fy, Fz,
                          hx, dt, clight, i, idx, 0);
        // if rank >= 2
        update_fluxes_1D( &N_update, &Fy_update, N, Fx, Fy, Fz,
                          hy, dt, clight, i, idy, 1);
        // if rank >= 3
        update_fluxes_1D( &N_update, &Fz_update, N, Fx, Fy, Fz,
                          hz, dt, clight, i, idz, 2);

        // get updated fluxes
        Fxnew[i] += Fx_update;
        Fynew[i] += Fy_update;
        Fznew[i] += Fz_update;
         

        // now get updated photon densities
        Nnew[i] += N_update;
       
      }
    }   
  } 
  
  // now copy values over (can eliminate this loop by passing *new variables as input parameters)
  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<my-gx; ix++) {
        int i = INDEX(ix,iy,iz,mx,my); //index of current cell
        N [i] = Nnew [i];
        Fx[i] = Fxnew[i];
        Fy[i] = Fynew[i];
        Fz[i] = Fznew[i];
     
#ifdef DEBUG_RT_TRANSPORT_STEP
        std::cout << " i: "     << i     
                  << " N [i]: " << N [i]
                  << " Fx[i]: " << Fx[i]
                  << " Fy[i]: " << Fy[i]
                  << " Fz[i]: " << Fz[i] << std::endl;     
#endif
      }
    }
  }
}


//----------------------------------------------------------------------


void EnzoMethodRadiativeTransfer::thermochemistry ( Block * block ) throw()
{

}

//======================================================================

void EnzoMethodRadiativeTransfer::compute_ (Block * block) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  
  const int N_groups = enzo_config->method_radiative_transfer_N_groups;
  const double min_freq = enzo_config->method_radiative_transfer_min_freq;
  const double max_freq = enzo_config->method_radiative_transfer_max_freq;
  const std::string flux_function = enzo_config->method_radiative_transfer_flux_function;
  const double clight = enzo_config->method_radiative_transfer_clight;


  // implement for only one frequency bin for now
  // Will need to call these functions in a loop over all the frequency groups
  // Need to get a frequency spectrum for the photons somehow (which SED???)
 
#ifdef DEBUG_CHECK_TIMESTEP
  std::cout << "calculated dt: " << this->timestep(block) << " actual dt: " << block->dt() << std::endl;
#endif 
  this->inject_photons(block);
  this->transport_photons(block, clight);
  this->thermochemistry(block); 

}
