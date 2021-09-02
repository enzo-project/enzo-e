// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRadiativeTransfer.cpp
/// @author   William Hicks (whicks@ucsd.edu)
/// @date     Mon Aug 16 17:05:23 PDT 2021
/// @brief    Implements the EnzoMethodRadiativeTransfer  class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodRadiativeTransfer ::EnzoMethodRadiativeTransfer()
  : Method()
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

  // define required fields if they do not exist
  this->define_fields();

  // Initialize default Refresh object

  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("photon_density");

  if (rank >= 1) {
    refresh->add_field("flux_x");
    refresh->add_field("P00"); // elements of the pressure tensor
  }
  if (rank >= 2) {
    refresh->add_field("flux_y");
    refresh->add_field("P10");
    refresh->add_field("P01");
    refresh->add_field("P11");
  }
  if (rank >= 3) {
    refresh->add_field("flux_z");
    refresh->add_field("P02");
    refresh->add_field("P12");
    refresh->add_field("P20");
    refresh->add_field("P21");
    refresh->add_field("P22");
  }
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

        (*n)[0] = Fx[i]/Fnorm;
        (*n)[1] = Fy[i]/Fnorm; 
        (*n)[2] = Fz[i]/Fnorm;
}

void EnzoMethodRadiativeTransfer::get_pressure_tensor (Block * block, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz, 
                       double clight) throw()
{
  Field field = block->data()->field();
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
      get_reduced_variables(&chi, &n, i, clight, N, Fx, Fy, Fz);
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

void EnzoMethodRadiativeTransfer::get_U_update (Block * block, double * N_update, 
                       double * Fx_update, double * Fy_update, double * Fz_update, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz, 
                       double hx, double hy, double hz, double dt, double clight, 
                       int i, int idx, int idy, int idz) throw()
{
  Field field = block->data()->field();
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
                                   Fx[i], Fx[i+idx], Fx[i-idx], 
                                   clight );
    
  *Fx_update += dt/hx * deltaQ_faces(Fx[i], Fx[i+idx], Fx[i-idx],
                                   clight*clight*P00[i],
                                   clight*clight*P00[i+idx],
                                   clight*clight*P00[i-idx], clight);

  // if rank >= 2
  *N_update += dt/hy * deltaQ_faces( N[i],  N[i+idy],  N[i-idy],
                                   Fy[i], Fy[i+idy], Fy[i-idy], 
                                   clight );
    
  *Fx_update += dt/hy * deltaQ_faces(Fx[i], Fx[i+idy], Fx[i-idy],
                                   clight*clight*P10[i],
                                   clight*clight*P10[i+idy],
                                   clight*clight*P10[i-idy], clight);

  *Fy_update += dt/hx * deltaQ_faces(Fy[i], Fy[i+idx], Fy[i-idx],
                                   clight*clight*P01[i],
                                   clight*clight*P01[i+idx],
                                   clight*clight*P01[i-idx], clight);

  *Fy_update += dt/hy * deltaQ_faces(Fy[i], Fy[i+idy], Fy[i-idy],
                                   clight*clight*P11[i],
                                   clight*clight*P11[i+idy],
                                   clight*clight*P11[i-idy], clight);


   // if rank >= 3
  *N_update += dt/hz * deltaQ_faces( N[i],  N[i+idz],  N[i-idz],
                                   Fz[i], Fz[i+idz], Fz[i-idz], 
                                   clight );  

  *Fx_update += dt/hz * deltaQ_faces(Fx[i], Fx[i+idz], Fx[i-idz],
                                   clight*clight*P20[i],
                                   clight*clight*P20[i+idz],
                                   clight*clight*P20[i-idz], clight);

  *Fy_update += dt/hz * deltaQ_faces(Fy[i], Fy[i+idz], Fy[i-idz],
                                   clight*clight*P21[i],
                                   clight*clight*P21[i+idz],
                                   clight*clight*P21[i-idz], clight);

  *Fz_update += dt/hx * deltaQ_faces(Fz[i], Fz[i+idx], Fz[i-idx],
                                   clight*clight*P02[i],
                                   clight*clight*P02[i+idx],
                                   clight*clight*P02[i-idx], clight);

  *Fz_update += dt/hy * deltaQ_faces(Fz[i], Fz[i+idy], Fz[i-idy],
                                   clight*clight*P12[i],
                                   clight*clight*P12[i+idy],
                                   clight*clight*P12[i-idy], clight);

  *Fz_update += dt/hz * deltaQ_faces(Fz[i], Fz[i+idz], Fz[i-idz],
                                   clight*clight*P22[i],
                                   clight*clight*P22[i+idz],
                                   clight*clight*P22[i-idz], clight);


}

void EnzoMethodRadiativeTransfer::transport_photons ( Block * block, double clight ) throw()
{
  // TODO: Make Nnew, Fxnew, Fynew, Fznew input parameters so that we don't need
  // the loop at the end that copies the values over to the original fields 
  // (Like in EnzoMethodHeat::compute_)

  // TODO: Adapt for 2D and 3D cases
  
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

  //calculate the pressure tensor
  get_pressure_tensor(block, N, Fx, Fy, Fz, clight);
 
  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
        int i = INDEX(ix,iy,iz,mx,my); //index of current cell
        double N_update=0, Fx_update=0, Fy_update=0, Fz_update=0;
      
        // TODO: Just pass in *new variable instead of update variable so you don't
        // have to do the += steps after
        get_U_update( block, &N_update, &Fx_update, &Fy_update, &Fz_update,
                             N, Fx, Fy, Fz, hx, hy, hz, dt, clight,
                             i, idx, idy, idz ); 
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
 
  this->inject_photons(block);
  this->transport_photons(block, clight);
  this->thermochemistry(block); 

}
