// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-07
/// @brief    Implements the EnzoMethodGravity class


#include "cello.hpp"
#include "enzo.hpp"

#include "enzo.decl.h"

// #define DEBUG_COPY_B
// #define DEBUG_COPY_POTENTIAL

//----------------------------------------------------------------------

EnzoMethodGravity::EnzoMethodGravity
(int index_solver,
 double grav_const,
 int order,
 bool accumulate)
  : Method(),
    index_solver_(index_solver),
    grav_const_(grav_const),
    order_(order),
    ir_exit_(-1)
{


  // Change this if fields used in this routine change
  // declare required fields
  this->required_fields_ = std::vector<std::string>
                           {"density","density_total","B","potential",
                            "acceleration_x","acceleration_y","acceleration_z"};
#ifdef DEBUG_FIELD_FACE
  this->required_fields_.insert(this->required_fields_.end(),
                                {"debug_1","debug_2"});
#endif
#ifdef DEBUG_COPY_B
  this->required_fields_.push_back("B_copy");
#endif
#ifdef DEBUG_COPY_POTENTIAL
  this->required_fields_.push_back("potential_copy");
#endif
#ifdef DEBUG_COPY_DENSITY
  this->required_fields_.push_back("density_total_copy");
#endif
#ifdef READ_ENZO_POTENTIAL
  this->required_fields_.push_back({"potential_enzo","potential_dff"});
#endif

  if (accumulate){
    this->required_fields_.insert(this->required_fields_.end(),
                                  {"density_particle","density_particle_accumulate"});
  }

  // now define fields if they do not exist
  this->define_fields();


  // Refresh adds density_total field faces and one layer of ghost
  // zones to "B" field

  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("acceleration_x");
  refresh->add_field("acceleration_y");
  refresh->add_field("acceleration_z");
  refresh->add_field("density");
  // Accumulate is used when particles are deposited into density_total

  if (accumulate) {
    refresh->set_accumulate(true);
    refresh->add_field_src_dst
      ("density_particle","density_particle_accumulate");
    refresh->add_field_src_dst("density_total","B");
  }

  ir_exit_ = add_new_refresh_();
  cello::simulation()->new_refresh_set_name(ir_post_,name()+":exit");
  Refresh * refresh_exit = cello::refresh(ir_exit_);

  refresh_exit->add_field("potential");

  refresh_exit->set_callback(CkIndex_EnzoBlock::p_method_gravity_end());
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{
  // Initialize the linear system

  Field field = block->data()->field();

  /// access problem-defining fields for eventual RHS and solution
  const int ib  = field.field_id ("B");
  const int id  = field.field_id("density");
  const int idt = field.field_id("density_total");
  const int idensity = (idt != -1) ? idt : id;

  // Solve the linear system
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  const int m = mx*my*mz;

  enzo_float * B = (enzo_float*) field.values (ib);
#ifdef DEBUG_COPY_B
  const int ib_copy = field.field_id ("B_copy");
  enzo_float * B_copy = (enzo_float*) field.values (ib_copy);
#endif
  enzo_float * D = (enzo_float*) field.values (idensity);

  for (int i=0; i<m; i++) D[i] += B[i];

  // Add density_particle values to density_particle_accumulate ghosts

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (block->is_leaf()) {
    if (cosmology) {

      int gx,gy,gz;
      field.ghost_depth(0,&gx,&gy,&gz);
      gx=gy=gz=0;
      for (int iz=gz; iz<mz-gz; iz++) {
	for (int iy=gy; iy<my-gy; iy++) {
	  for (int ix=gx; ix<mx-gx; ix++) {
	    int i = ix + mx*(iy + my*iz);
	    D[i]=-(D[i]-1.0);
	    B[i]  = D[i];
	  }
	}
      }
    } else {
      field.scale(ib, -4.0 * (cello::pi) * grav_const_, idensity);
    }

  } else {

    for (int i=0; i<mx*my*mz; i++) B[i] = 0.0;

  }

#ifdef DEBUG_COPY_B
  for (int i=0; i<m; i++) B_copy[i] = B[i];
#endif

  Solver * solver = enzo::problem()->solver(index_solver_);

  // May exit before solve is done...
  solver->set_callback (CkIndex_EnzoBlock::p_method_gravity_continue());

  const int ix = field.field_id ("potential");

  std::shared_ptr<Matrix> A (std::make_shared<EnzoMatrixLaplace>(order_));

  solver->set_field_x(ix);
  solver->set_field_b(ib);

  solver->apply (A, block);
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_gravity_continue()
{
  // So do refresh with barrier synch (note barrier instead of
  // neighbor synchronization otherwise will conflict with Method
  // refresh ("Charm++ fatal error: mis-matched client callbacks in
  // reduction messages")

  EnzoMethodGravity * method = static_cast<EnzoMethodGravity*> (this->method());
  method->refresh_potential(this);

}

//----------------------------------------------------------------------
void EnzoMethodGravity::refresh_potential (EnzoBlock * enzo_block) throw()
{
  cello::refresh(ir_exit_)->set_active(enzo_block->is_leaf());
  enzo_block->new_refresh_start
    (ir_exit_, CkIndex_EnzoBlock::p_method_gravity_end());
}
//----------------------------------------------------------------------

void EnzoBlock::p_method_gravity_end()
{
  EnzoMethodGravity * method = static_cast<EnzoMethodGravity*> (this->method());
  method->compute_accelerations(this);
  // wait for all Blocks before continuing
  compute_done();
}

void EnzoMethodGravity::compute_accelerations (EnzoBlock * enzo_block) throw()
{

  Field field = enzo_block->data()->field();
  int gx,gy,gz;
  int mx,my,mz;
  field.ghost_depth(0,&gx,&gy,&gz);
  field.dimensions (0,&mx,&my,&mz);
  const int m = mx*my*mz;
  enzo_float * potential = (enzo_float*) field.values ("potential");

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (cosmology) {

    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt   = enzo_block->timestep();
    double time = enzo_block->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    //    cosmology-> compute_expansion_factor (&a,&dadt,time);

    for (int i=0; i<m; i++) potential[i] /= cosmo_a;
  }

  /// compute acceleration fields from potential

  EnzoComputeAcceleration compute_acceleration(cello::rank(), order_);

  compute_acceleration.compute(enzo_block);

  // Clear "B" and "density_total" fields for next call
  // Note density_total may not be defined

  enzo_float * B = (enzo_float*) field.values("B");
  for (int i=0; i<m; i++) B[i] = 0.0;

  enzo_float * de_t = (enzo_float*) field.values("density_total");
  if (de_t) for (int i=0; i<m; i++) de_t[i] = 0.0;

  if (potential) {
    for (int i=0; i<m; i++) potential[i] = 0.0;
  }

}

//----------------------------------------------------------------------

double EnzoMethodGravity::timestep (Block * block) const throw()
{
  return timestep_(block);
}

//----------------------------------------------------------------------

double EnzoMethodGravity::timestep_ (Block * block) const throw()
{
  Field field = block->data()->field();

  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");

  enzo_float dt = std::numeric_limits<enzo_float>::max();

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (cosmology) {
    const int rank = cello::rank();
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt   = block->dt();
    double time = block->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    if (rank >= 1) hx*=cosmo_a;
    if (rank >= 2) hy*=cosmo_a;
    if (rank >= 3) hz*=cosmo_a;
  }

  if (ax) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hx/(fabs(ax[i]+1e-20)))));
	}
      }
    }
  }
  if (ay) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hy/(fabs(ay[i]+1e-20)))));
	}
      }
    }
  }
  if (az) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hz/(fabs(az[i]+1e-20)))));
	}
      }
    }
  }

  return 0.5*dt;
}
