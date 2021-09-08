// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRadiationInjection.cpp
/// @author   William Hicks (whicks@ucsd.edu)
/// @date     Mon Aug 16 17:05:23 PDT 2021
/// @brief    Implements the EnzoMethodRadiationInjection  class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodRadiationInjection ::EnzoMethodRadiationInjection()
  : Method()
{

  this->required_fields_ = std::vector<std::string> {"photon_density"}; //number density

  // define required fields if they do not exist
  this->define_fields();

  // Initialize default Refresh object

  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("photon_density");

}

//----------------------------------------------------------------------

void EnzoMethodRadiationInjection ::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
}

//----------------------------------------------------------------------

void EnzoMethodRadiationInjection::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    compute_ (block);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodRadiationInjection::timestep ( Block * block ) const throw()
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
  const double clight = enzo_config->method_radiation_injection_clight;


  return h_min / (3*clight);
}


//----------------------------------------------------------------------
void EnzoMethodRadiationInjection::compute_ ( Block * block ) throw()
{
  // Solve dN_i/dt = Ndot^*_i
  // routine for identifying star particles and getting their grid position
  // copy/pasted from enzo_EnzoMethodFeedback::compute_()

  const EnzoConfig * enzo_config = enzo::config();
  double injection_rate = 1e10; // enzo_config->method_radiation_injection_injection_rate;
  double f_esc = 1.0; // enzo_config->method_radiation_injection_escape_fraction;

  Field field = block->data()->field();
  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);
   
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);
 
  double xm,ym,zm;
  double xp,yp,zp;
  double hx,hy,hz;
  block->lower(&xm,&ym,&zm);
  block->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  // if (cosmology) {...}
 
  double inv_vol = 1.0/(hx*hy*hz); 

  double dt = timestep(block);

  // get relevant field variables
  enzo_float * N  = (enzo_float *) field.values("photon_density");

  Particle particle = block->data()->particle();
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

  // loop through batches
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
      // TODO: Write function that references an SED to calculate injection_rate
      // call that function here.
      // This function will take in particle mass as a parameter, and will fit
      // and integrate over the SED to get the total injection rate. The current
      // implementation in RAMSES assumes that SED of gas stays constant.
      // How do I keep track of the gas SED if different mass stars radiate differently?
      N[i] += f_esc * injection_rate * dt * inv_vol; // photons/s * timestep / volume 
      // I don't have to directly alter the fluxes here because that naturally gets
      // taken care of during the transport step 
    }
  } 
}
//======================================================================
