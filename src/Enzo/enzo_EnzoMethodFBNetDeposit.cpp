// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodFBNetDeposit.hpp
/// @author   William Hicks (whicks@ucsd.edu) 
/// @date     Wed Jun 7 16:14:38 PDT 2023
/// @brief    [\ref Enzo] Implements of EnzoMethodFBNetDeposit
///           method.


#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodFBNetDeposit::EnzoMethodFBNetDeposit ()
  : Method()
{
  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);

  // define the metal fields
  refresh->add_field("metal_density");
  refresh->add_field("PopIII_metal_density");
  refresh->add_field("PopIII_SNe_metal_density");
  refresh->add_field("PopIII_HNe_metal_density");
  refresh->add_field("PopIII_PISNe_metal_density");

  refresh->add_particle(cello::particle_descr()->type_index("popIIIremnant"));
}

//----------------------------------------------------------------------

void EnzoMethodFBNetDeposit::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
}

//----------------------------------------------------------------------

void EnzoMethodFBNetDeposit::compute ( Block * block ) throw()
{
  EnzoBlock * enzo_block = enzo::block(block);
  
  double sphere_x = 0.0, sphere_y = 0.0, sphere_z = 0.0;
  double sphere_r = 0.0;
  double sphere_mSNe = 0.0, sphere_mHNe = 0.0, sphere_mPISNe = 0.0;
  if (block->is_leaf()) {

    Field field = enzo_block->data()->field();

    int mx,my,mz;
    field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
    int gx,gy,gz;
    field.ghost_depth(0,&gx, &gy, &gz);

    double xm,ym,zm;
    double xp,yp,zp;
    double hx,hy,hz;
    enzo_block->data()->lower(&xm,&ym,&zm);
    enzo_block->data()->upper(&xp,&yp,&zp);
    field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

    Particle particle = enzo_block->data()->particle();

    const int it = particle.type_index("popIIIremnant");

    const int ia_x = particle.attribute_index (it, "x");
    const int ia_y = particle.attribute_index (it, "y");
    const int ia_z = particle.attribute_index (it, "z");
    const int ia_r = particle.attribute_index (it, "r");
    const int ia_vx = particle.attribute_index (it, "vx");
    const int ia_vy = particle.attribute_index (it, "vy");
    const int ia_vz = particle.attribute_index (it, "vz");
  
    const int ia_mSNe = particle.attribute_index(it, "yield_SNe");
    const int ia_mHNe = particle.attribute_index(it, "yield_HNe");
    const int ia_mPISNe = particle.attribute_index(it, "yield_PISNe");

    const int ia_t = particle.attribute_index(it, "creation_time");

    const int dp = particle.stride(it, ia_x);
    const int dr = particle.stride(it, ia_r);
    const int dv = particle.stride(it, ia_vx);
    const int dm_SNe = particle.stride(it, ia_mSNe);
    const int dm_HNe = particle.stride(it, ia_mHNe);
    const int dm_PISNe = particle.stride(it, ia_mPISNe);
    const int dt = particle.stride(it, ia_t);

    const int nb = particle.num_batches(it);

    // loop through batches
    for (int ib=0; ib<nb; ib++) {
      enzo_float * px   = (enzo_float *) particle.attribute_array(it, ia_x, ib);
      enzo_float * py   = (enzo_float *) particle.attribute_array(it, ia_y, ib);
      enzo_float * pz   = (enzo_float *) particle.attribute_array(it, ia_z, ib);
      enzo_float * pr   = (enzo_float *) particle.attribute_array(it, ia_r, ib);
      enzo_float * pvx  = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
      enzo_float * pvy  = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
      enzo_float * pvz  = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
      enzo_float * pmSNe = (enzo_float *) particle.attribute_array(it, ia_mSNe, ib);
      enzo_float * pmHNe = (enzo_float *) particle.attribute_array(it, ia_mHNe, ib);
      enzo_float * pmPISNe = (enzo_float *) particle.attribute_array(it, ia_mPISNe, ib);
      enzo_float * pform  = (enzo_float *) particle.attribute_array(it, ia_t, ib);

      int np = particle.num_particles(it,ib);

      // loop through particles
      for (int ip=0; ip<np; ip++) {
        int ipdp = ip*dp;
        int ipdr = ip*dr;
        int ipdv = ip*dv;
        int ipdm_SNe = ip*dm_SNe;
        int ipdm_HNe = ip*dm_HNe;
        int ipdm_PISNe = ip*dm_PISNe;
        int ipdt = ip*dt;


        CkPrintf("EnzoMethodFBDeposit:: creation_time = %f; block time = %f; block_cycle = %d\n", pform[ipdt], block->time(), block->cycle());
        if (pform[ipdt] == block->time()) {
          sphere_x = px[ipdp]; 
          sphere_y = py[ipdp]; 
          sphere_z = pz[ipdp];
          sphere_r = pr[ipdp];
          sphere_mSNe = pmSNe[ipdm_SNe];
          sphere_mHNe = pmHNe[ipdm_HNe];
          sphere_mPISNe = pmPISNe[ipdm_PISNe]; 
        }
      }
    }

  }

  CkCallback callback (CkIndex_EnzoBlock::p_method_fbnet_update_mesh(NULL),
                       enzo::block(block)->proxy_array());

  // note that StarNet will spawn at most one sphere per inference block
  // in a timestep

  double center[3] = {sphere_x, sphere_y, sphere_z};

  EnzoObjectFeedbackSphere sphere(center, sphere_r, sphere_mSNe, sphere_mHNe, sphere_mPISNe);
  EnzoObjectFeedbackSphere sphere_[1] = {sphere};

  int n = 0;
  SIZE_ARRAY_TYPE(n,EnzoObjectFeedbackSphere,sphere_, 1);

  enzo_block->contribute(n, sphere_, CkReduction::set, callback);
}

//------------------------------------

void EnzoBlock::p_method_fbnet_update_mesh(CkReductionMsg * msg) 
{

  // put spheres down on the mesh (if inference_method == "starnet")
  if (this->is_leaf()) {
    CkReduction::setElement *current = (CkReduction::setElement*) msg->getData();
    int i=0;
    while(current != NULL) {
      EnzoObjectFeedbackSphere * sphere = (EnzoObjectFeedbackSphere *) &current->data;
      FBNet::update_mesh(this, *sphere);
      i += 1;
      current = current->next();
    }
  }

  compute_done();
}

//-----------------------------------

double EnzoMethodFBNetDeposit::timestep( Block * block ) throw()
{
  return std::numeric_limits<double>::max();
}
