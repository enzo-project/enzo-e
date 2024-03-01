// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-07
/// @brief    Implements the EnzoMethodGravity class


#define TRACE_SUPERCYCLE

#include "cello.hpp"
#include "enzo.hpp"

#include "charm_enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodGravity::EnzoMethodGravity
(int index_solver,
 double grav_const,
 int order,
 bool accumulate,
 int index_prolong,
 double dt_max,
 std::string type_super)
  : Method(),
    index_solver_(index_solver),
    grav_const_(grav_const),
    order_(order),
    ir_exit_(-1),
    index_prolong_(index_prolong),
    dt_max_(dt_max),
    type_super_()
{
  ASSERT1("EnzoMethodGravity::EnzoMethodGravity()",
          "Unknown type_super %s (must be \"potential\" or \"accelerations\")",
          type_super.c_str(),
          ((type_super == "potential") || (type_super == "accelerations")));
  
  if (type_super == "potential") {
    type_super_ = super_type_potential;
  } else if (type_super=="accelerations") {
    type_super_ = super_type_accelerations;
  }
  // Change this if fields used in this routine change
  // declare required fields
  cello::define_field ("density");
  cello::define_field ("density_total");
  cello::define_field ("B");
  cello::define_field ("potential");
  const int rank = cello::rank();
  if (rank >= 1) cello::define_field ("acceleration_x");
  if (rank >= 2) cello::define_field ("acceleration_y");
  if (rank >= 3) cello::define_field ("acceleration_z");

  if (max_supercycle_ > 1) {
    cello::define_field ("potential_curr");
    cello::define_field ("potential_prev");
    if (rank >= 1) {
      cello::define_field ("acceleration_x_curr");
      cello::define_field ("acceleration_x_prev");
    }
    if (rank >= 2) {
      cello::define_field ("acceleration_y_curr");
      cello::define_field ("acceleration_y_prev");
    }
    if (rank >= 3) {
      cello::define_field ("acceleration_z_curr");
      cello::define_field ("acceleration_z_prev");
    }
  }

  if (accumulate){
    cello::define_field ("density_particle");
    cello::define_field ("density_particle_accumulate");
  }

  // Define Block scalars

  ScalarDescr * scalar_descr_double = cello::scalar_descr_double();
  is_time_curr_ = scalar_descr_double->new_value(name()+"_time_curr");
  is_time_prev_ = scalar_descr_double->new_value(name()+"_time_prev");

  // Refresh adds density_total field faces and one layer of ghost
  // zones to "B" field

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->set_prolong(index_prolong_);

  refresh_add_potentials_   (refresh);
  refresh_add_accelerations_(refresh);
  // Accumulate is used when particles are deposited into density_total

  if (accumulate) {

    // EnzoProlong does not work with accumulate!
    std::string prolong_name =
      cello::problem()->prolong(index_prolong_)->name();

    ASSERT1("EnzoMethodGravity::EnzoMethodGravity()",
            "Requesting accumulated particle mass refresh: "
            "rerun with parameter Method : %s : prolong = \"linear\"",
            name().c_str(),
            (prolong_name != "enzo"));

    refresh->set_accumulate(true);
    refresh->add_field_src_dst
      ("density_particle","density_particle_accumulate");
    refresh->add_field_src_dst("density_total","B");
  }

  ir_exit_ = add_refresh_();
  cello::simulation()->refresh_set_name(ir_exit_,name()+":exit");
  Refresh * refresh_exit = cello::refresh(ir_exit_);
  refresh_exit->set_prolong(index_prolong_);

  refresh_add_potentials_   (refresh_exit);
  refresh_add_accelerations_(refresh_exit);

  refresh_exit->set_callback(CkIndex_EnzoBlock::p_method_gravity_end());

  ASSERT1 ("EnzoMethodGravity::EnzoMethodGravity()",
           "type_super parameter must be either 0 (extrapolate potentials) or "
           "1 (extrapolate accelerations", type_super_,
           (type_super_ == super_type_potential) ||
           (type_super_ == super_type_accelerations));
}

//----------------------------------------------------------------------

void EnzoMethodGravity::refresh_add_potentials_(Refresh * refresh)
{
  refresh->add_field("potential");
  if (max_supercycle_ > 1) {
    refresh->add_field("potential_curr");
    refresh->add_field("potential_prev");
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::refresh_add_accelerations_(Refresh * refresh)
{
  const int rank = cello::rank();
  if (rank >= 1) {
    refresh->add_field("acceleration_x");
    if (max_supercycle_ > 1) {
      refresh->add_field("acceleration_x_curr");
      refresh->add_field("acceleration_x_prev");
    }
  }
  if (rank >= 2) {
    refresh->add_field("acceleration_y");
    if (max_supercycle_ > 1) {
      refresh->add_field("acceleration_y_curr");
      refresh->add_field("acceleration_y_prev");
    }
  }
  if (rank >= 3) {
    refresh->add_field("acceleration_z");
    if (max_supercycle_ > 1) {
      refresh->add_field("acceleration_z_curr");
      refresh->add_field("acceleration_z_prev");
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{
    EnzoBlock * enzo_block = static_cast<EnzoBlock*>(block);

    if (enzo::simulation()->state()->cycle() == enzo::config()->initial_cycle) {
    // Check if the pm_deposit method is being used and precedes the
    // gravity method.
    ASSERT("EnzoMethodGravity",
           "Error: pm_deposit method must precede gravity method.",
           enzo::problem()->method_precedes("pm_deposit", "gravity"));
  }
  // Initialize the linear system

#ifdef TRACE_SUPERCYCLE
  if (block->index().is_root()) {
    const auto & method_state = block->state()->method(index());
    CkPrintf ("TRACE_SUPERCYCLE cycle %d step %d num_steps %d\n",
              block->state()->cycle(),
              method_state.step(),
              method_state.num_steps());
  }
#endif

  // Solve only if this method's current substep is 0
  bool solve_this_step = (block->state()->method(index()).step() == 0);

  if (solve_this_step) {

#ifdef TRACE_SUPERCYCLE
    if (block->index().is_root()) {
      CkPrintf ("TRACE_SUPERCYCLE cycle %d: solve\n",block->state()->cycle());
    }
#endif

    Field field = block->data()->field();
    /// access problem-defining fields for eventual RHS and solution
    const int ib  = field.field_id ("B");
    const int id  = field.field_id("density");
    const int idt = field.field_id("density_total");
    ASSERT ("EnzoMethodGravity::compute",
            "missing required field density_total",
            idt != -1);
    // Solve the linear system
    int mx,my,mz;
    int gx,gy,gz;
    field.dimensions (0,&mx,&my,&mz);
    field.ghost_depth(0,&gx,&gy,&gz);

    const int m = mx*my*mz;
    enzo_float * B = (enzo_float*) field.values (ib);
    enzo_float * DT = (enzo_float*) field.values (idt);

    for (int i=0; i<m; i++) {
      DT[i] += B[i];
    }

    // Add density_particle values to density_particle_accumulate ghosts

    if (block->is_leaf()) {
      EnzoPhysicsCosmology * cosmology = enzo::cosmology();
      if (cosmology) {
        int gx,gy,gz;
        field.ghost_depth(0,&gx,&gy,&gz);
        gx=gy=gz=0;
        for (int iz=gz; iz<mz-gz; iz++) {
          for (int iy=gy; iy<my-gy; iy++) {
            for (int ix=gx; ix<mx-gx; ix++) {
              int i = ix + mx*(iy + my*iz);
              // In cosmological simulations, density units are defined such that `rho_bar_m` is
              // 1.0, and time units are defined such that `4 * pi * G * rho_bar_m` is 1.0, where
              // `G` is the gravitational constant, and `rho_bar_m` is the mean matter density
              // of the universe. These choices of units result in Poisson's equation having a
              // much simplified form.
              DT[i]=-(DT[i]-1.0);
              B[i]  = DT[i];
            }
          }
        }
      } else {
        const double scale = -4.0 * (cello::pi) * grav_const_;
        for (int iz=0; iz<mz; iz++) {
          for (int iy=0; iy<my; iy++) {
            for (int ix=0; ix<mx; ix++) {
              const int i = ix + mx*(iy + my*iz);
              B[i] = scale * DT[i];
            }
          }
        }
      }

    } else {

      for (int i=0; i<mx*my*mz; i++) B[i] = 0.0;

    }

    Solver * solver = enzo::problem()->solver(index_solver_);

    // May exit before solve is done...
    solver->set_callback (CkIndex_EnzoBlock::p_method_gravity_continue());

    // Save previous potential

    if ((max_supercycle_ > 1) &&
        (type_super_ == super_type_potential)) {
      super_save_potential_(enzo_block);
    }

    int ix = (max_supercycle_ > 1) ?
      field.field_id("potential_curr") :
      field.field_id("potential");

    ASSERT("EnzoMethodGravity::compute()",
           "max_supercycle > 1 but potential_curr field not defined",
           (ix >= 0));

    std::shared_ptr<Matrix> A (std::make_shared<EnzoMatrixLaplace>(order_));
    solver->set_field_x(ix);
    solver->set_field_b(ib);
    solver->apply (A, block);

  } else { // ( ! solve_this_step )

#ifdef TRACE_SUPERCYCLE
    if (block->index().is_root()) {
      CkPrintf ("TRACE_SUPERCYCLE cycle %d: extrapolate\n",block->state()->cycle());
    }
#endif

    // Skip solve; compute accelerations using extrapolated potential
    if (enzo_block && enzo_block->is_leaf()) {
      compute_accelerations(enzo_block);
    }

    enzo_block->compute_done();
  }
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
  enzo_block->refresh_start
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

//----------------------------------------------------------------------

void EnzoMethodGravity::super_save_potential_(EnzoBlock * enzo_block)
{
  Field field = enzo_block->data()->field();
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int m = mx*my*mz;
  enzo_float * potential_curr = (enzo_float*) field.values ("potential_curr");
  enzo_float * potential_prev = (enzo_float*) field.values ("potential_prev");

  s_time_prev_(enzo_block) = s_time_curr_(enzo_block);

  if (potential_curr && potential_prev) {
    for (int i=0; i<m; i++) { potential_prev[i] = potential_curr[i]; }
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::super_extrapolate_potential_(EnzoBlock * enzo_block)
{
  // extrapolate step: extropolate potential from potential_prev and potential_curr
  double tp = s_time_prev_(enzo_block);
  double tc = s_time_curr_(enzo_block);
  double t = enzo_block->state()->time();

  // compute extrapolation coefficients
  const double cp = (tc - t) / (tc - tp);
  const double cc = (t - tp) / (tc - tp);

  // perform potential extrapolation
  Field field = enzo_block->data()->field();
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int m = mx*my*mz;
  enzo_float * potential      = (enzo_float*) field.values ("potential");
  enzo_float * potential_curr = (enzo_float*) field.values ("potential_curr");
  enzo_float * potential_prev = (enzo_float*) field.values ("potential_prev");
  ASSERT1("EnzoMethodGravity::compute_accelerations()",
          "max_supercycle = %d != 1 but potential_[curr|prev] fields not defined",
          max_supercycle_,
          (potential_curr && potential_prev));
  for (int i=0; i<m; i++) {
    potential[i] = cc*potential_curr[i] + cp*potential_prev[i];
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::super_extrapolate_accelerations_(EnzoBlock * enzo_block)
{
  // extrapolate step: extropolate acceleration_[xyz] from
  // potential_[xyz]_prev and potential_[xyz]_curr
  double tp = s_time_prev_(enzo_block);
  double tc = s_time_curr_(enzo_block);
  double t = enzo_block->state()->time();

  // compute extrapolation coefficients
  const double cp = (tc - t) / (tc - tp);
  const double cc = (t - tp) / (tc - tp);


  // perform acceleration_x extrapolation
  Field field = enzo_block->data()->field();
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int m = mx*my*mz;

  enzo_float * ax      = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ax_curr = (enzo_float*) field.values ("acceleration_x_curr");
  enzo_float * ax_prev = (enzo_float*) field.values ("acceleration_x_prev");
  ASSERT1("EnzoMethodGravity::super_extrapolate_accelerations()",
          "max_supercycle = %d != 1 but acceleration_x_[curr|prev] fields not defined",
          max_supercycle_,
          (ax_curr && ax_prev));
  for (int i=0; i<m; i++) {
    ax[i] = cc*ax_curr[i] + cp*ax_prev[i];
  }

  const int rank = cello::rank();
  enzo_float * ay      = (enzo_float*) field.values ("acceleration_y");
  enzo_float * ay_curr = (enzo_float*) field.values ("acceleration_y_curr");
  enzo_float * ay_prev = (enzo_float*) field.values ("acceleration_y_prev");
  ASSERT1("EnzoMethodGravity::super_extrapolate_accelerations()",
          "max_supercycle = %d != 1 but acceleration_y_[curr|prev] fields not defined",
          max_supercycle_,
          rank < 2 || (ay_curr && ay_prev));
  if (rank >= 2) {
    for (int i=0; i<m; i++) {
      ay[i] = cc*ay_curr[i] + cp*ay_prev[i];
    }
  }

  enzo_float * az      = (enzo_float*) field.values ("acceleration_z");
  enzo_float * az_curr = (enzo_float*) field.values ("acceleration_z_curr");
  enzo_float * az_prev = (enzo_float*) field.values ("acceleration_z_prev");
  ASSERT1("EnzoMethodGravity::super_extrapolate_accelerations()",
          "max_supercycle = %d != 1 but acceleration_z_[curr|prev] fields not defined",
          max_supercycle_,
          rank < 3 || (az_curr && az_prev));
  if (rank >= 3) {
    for (int i=0; i<m; i++) {
      az[i] = cc*az_curr[i] + cp*az_prev[i];
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::super_save_accelerations_(EnzoBlock * enzo_block)
{
  Field field = enzo_block->data()->field();
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int m = mx*my*mz;
  enzo_float * ax_curr = (enzo_float*) field.values ("acceleration_x_curr");
  enzo_float * ax_prev = (enzo_float*) field.values ("acceleration_x_prev");
  enzo_float * ay_curr = (enzo_float*) field.values ("acceleration_y_curr");
  enzo_float * ay_prev = (enzo_float*) field.values ("acceleration_y_prev");
  enzo_float * az_curr = (enzo_float*) field.values ("acceleration_z_curr");
  enzo_float * az_prev = (enzo_float*) field.values ("acceleration_z_prev");

  s_time_prev_(enzo_block) = s_time_curr_(enzo_block);

  if (ax_curr && ax_prev) {
    for (int i=0; i<m; i++) { ax_prev[i] = ax_curr[i]; }
  }
  if (ay_curr && ay_prev) {
    for (int i=0; i<m; i++) { ay_prev[i] = ay_curr[i]; }
  }
  if (az_curr && az_prev) {
    for (int i=0; i<m; i++) { az_prev[i] = az_curr[i]; }
  }
}
//----------------------------------------------------------------------

void EnzoMethodGravity::compute_accelerations (EnzoBlock * enzo_block) throw()
{

  Field field = enzo_block->data()->field();
  enzo_float * potential      = (enzo_float*) field.values ("potential");
  enzo_float * potential_prev = (enzo_float*) field.values ("potential_prev");
  enzo_float * potential_curr = (enzo_float*) field.values ("potential_curr");

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int m = mx*my*mz;

  // Extrapolate potential from "potential_curr" and "potential_prev"
  // if supercycling and not a solve step

  const bool is_root = enzo_block->index().is_root();

  if (max_supercycle_ > 1) {
    for (int i=0; i<m; i++) {
      potential[i] = potential_curr[i];
    }
  }

  if (type_super_ == super_type_potential) {

    if (max_supercycle_ > 1) {

      const auto & method_state = enzo_block->state()->method(index());
      const double step = method_state.step();

      if (step == 0) {

        // solve step: copy computed potential_curr to potential
        for (int i=0; i<m; i++) {
          potential[i] = potential_curr[i];
        }
        s_time_curr_(enzo_block) = enzo_block->state()->time();

      } else {

        super_extrapolate_potential_(enzo_block);

      }
    }
  }

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (cosmology) {
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    auto dt   = enzo_block->state()->dt();
    auto time = enzo_block->state()->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    for (int i=0; i<m; i++) potential[i] /= cosmo_a;
  }

  /// compute acceleration fields from potential

  EnzoComputeAcceleration compute_acceleration(cello::rank(), order_);

  if ((max_supercycle_ > 1) &&
      (type_super_ == super_type_accelerations)) {
    const auto & method_state = enzo_block->state()->method(index());
    const double step = method_state.step();
    enzo_float * ax =      (enzo_float*) field.values ("acceleration_x");
    enzo_float * ay =      (enzo_float*) field.values ("acceleration_y");
    enzo_float * az =      (enzo_float*) field.values ("acceleration_z");
    enzo_float * ax_curr = (enzo_float*) field.values ("acceleration_x_curr");
    enzo_float * ay_curr = (enzo_float*) field.values ("acceleration_y_curr");
    enzo_float * az_curr = (enzo_float*) field.values ("acceleration_z_curr");

    if (step == 0) {

      super_save_accelerations_(enzo_block);

      compute_acceleration.compute(enzo_block);

      s_time_curr_(enzo_block) = enzo_block->state()->time();
      if (ax && ax_curr) for (int i=0; i<m; i++) ax_curr[i] = ax[i];
      if (ay && ay_curr) for (int i=0; i<m; i++) ay_curr[i] = ay[i];
      if (az && az_curr) for (int i=0; i<m; i++) az_curr[i] = az[i];

    } else { // not a solve step

      super_extrapolate_accelerations_(enzo_block);
    }

  } else {

    compute_acceleration.compute(enzo_block);

  }

  // Clear "B" and "density_total" fields for next call
  // Note density_total may not be defined

  enzo_float * B = (enzo_float*) field.values("B");
  for (int i=0; i<m; i++) B[i] = 0.0;

  enzo_float * de_t = (enzo_float*) field.values("density_total");
  if (de_t) for (int i=0; i<m; i++) de_t[i] = 0.0;

}

//----------------------------------------------------------------------

double EnzoMethodGravity::timestep (Block * block) throw()
{
  return timestep_(block);
}

//----------------------------------------------------------------------

double EnzoMethodGravity::timestep_ (Block * block) throw()
{
  Field field = block->data()->field();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");

  const int rank = cello::rank();
  
  enzo_float dt = std::numeric_limits<enzo_float>::max();

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);
  
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (cosmology) {
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt = block->state()->dt();
    double time = block->state()->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    if (rank >= 1) hx*=cosmo_a;
    if (rank >= 2) hy*=cosmo_a;
    if (rank >= 3) hz*=cosmo_a;
  }
  
  double mean_cell_width;

  if (rank == 1) mean_cell_width = hx;
  if (rank == 2) mean_cell_width = sqrt(hx*hy);
  if (rank == 3) mean_cell_width = cbrt(hx*hy*hz);
  
  // Timestep is sqrt(mean_cell_width / (a_mag_max + epsilon)),
  // where a_mag_max is the maximum acceleration magnitude
  // across all cells in the block, and epsilon defined as
  // mean_cell_width / dt_max_^2. This means that when acceleration
  // is zero everywhere, the timestep is equal to dt_max_
  
  const double epsilon = mean_cell_width / (dt_max_ * dt_max_);

  // Find th maximum of the square of the magnitude of acceleration
  // across all active cells, then get the square root of this value
  
  double a_mag_2_max = 0.0;
  double a_mag_2;

  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
	int i=ix + mx*(iy + iz*my);
	if (rank == 1) a_mag_2 = ax[i] * ax[i];
	if (rank == 2) a_mag_2 = ax[i] * ax[i] + ay[i] * ay[i];
	if (rank == 3) a_mag_2 = ax[i] * ax[i] + ay[i] * ay[i]
			       + az[i] * az[i];
	a_mag_2_max = std::max(a_mag_2_max,a_mag_2);
      }
    }
  }

  const double a_mag_max = sqrt(a_mag_2_max);
  dt = sqrt(mean_cell_width / (a_mag_max + epsilon)) ;

  return 0.5*dt;
}
