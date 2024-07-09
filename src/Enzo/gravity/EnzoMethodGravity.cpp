// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-07
/// @brief    Implements the EnzoMethodGravity class


#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/gravity/gravity.hpp"

// #define TRACE_SUPERCYCLE

//----------------------------------------------------------------------

EnzoMethodGravity::EnzoMethodGravity(ParameterGroup p, int index_solver,
                                     int index_prolong,
                                     int max_super,
                                     std::string type_super)
  : Method(),
    index_solver_(index_solver),
    order_(p.value_integer("order",4)),
    ir_exit_(-1),
    index_prolong_(index_prolong),
    type_super_(),
    dt_max_(p.value_float("dt_max",1.0e10))
{
  max_supercycle_ = max_super;
  // Initialize type_super_ super-cycling type
  if (type_super == "potential") {
    type_super_ = super_type_potential;
  } else if (type_super=="accelerations") {
    type_super_ = super_type_accelerations;
  } else {
    ERROR1("EnzoMethodGravity::EnzoMethodGravity()",
          "Unknown type_super %s (must be \"potential\" or \"accelerations\")",
           type_super.c_str());
  }

  // Declare required fields
  const bool accumulate = p.value_logical("accumulate",true);

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

  if (is_supercycle()) {

    if (type_super_ == super_type_potential) {

      super_define_fields_("potential","potential_curr", "potential_prev");

    } else if (type_super_ == super_type_accelerations) {

      if (rank >= 1) {
        super_define_fields_("acceleration_x",
                             "acceleration_x_curr",
                             "acceleration_x_prev");
      }
      if (rank >= 2) {
        super_define_fields_("acceleration_y",
                             "acceleration_y_curr",
                             "acceleration_y_prev");
      }
      if (rank >= 3) {
        super_define_fields_("acceleration_z",
                             "acceleration_z_curr",
                             "acceleration_z_prev");
      }
    }
  }

  if (accumulate){
    cello::define_field ("density_particle");
    cello::define_field ("density_particle_accumulate");
  }

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
  if ( is_supercycle_potential() ) {
    refresh->add_field("potential_curr");
    refresh->add_field("potential_prev");
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::refresh_add_accelerations_(Refresh * refresh)
{
  const int rank = cello::rank();
  if (rank >= 1) refresh->add_field("acceleration_x");
  if (rank >= 2) refresh->add_field("acceleration_y");
  if (rank >= 3) refresh->add_field("acceleration_z");
  if ( is_supercycle_accelerations() ) {
    if (rank >= 1) {
      refresh->add_field("acceleration_x_curr");
      refresh->add_field("acceleration_x_prev");
    }
    if (rank >= 2) {
      refresh->add_field("acceleration_y_curr");
      refresh->add_field("acceleration_y_prev");
    }
    if (rank >= 3) {
      refresh->add_field("acceleration_z_curr");
      refresh->add_field("acceleration_z_prev");
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*>(block);

  if (cello::is_initial_cycle(InitCycleKind::fresh_or_noncharm_restart)) {
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

    // Prepare right-hand side vector B

    for (int i=0; i<m; i++) B[i] = 0.0;

    if (block->is_leaf()) {

      EnzoPhysicsCosmology * cosmology = enzo::cosmology();

      if (cosmology) {
        for (int i=0; i<m; i++) {
          // In cosmological simulations, density units are defined
          // such that `rho_bar_m` is 1.0, and time units are defined
          // such that `4 * pi * G * rho_bar_m` is 1.0, where `G` is
          // the gravitational constant, and `rho_bar_m` is the mean
          // matter density of the universe. These choices of units
          // result in Poisson's equation having a much simplified
          // form.
          DT[i] = - (DT[i] - 1.0);
          B[i]  = DT[i];
        }

      } else { // ! cosmology

        const double scale = -4.0 * (cello::pi) * (enzo::grav_constant_codeU());
        for (int i=0; i<m; i++) {
          B[i] = scale * DT[i];
        }

      }

    }

    Solver * solver = enzo::problem()->solver(index_solver_);

    // May exit before solve is done...
    solver->set_callback (CkIndex_EnzoBlock::p_method_gravity_continue());

    // Save previous potential

    if ( is_supercycle_potential() ) {
      super_shift_fields_(enzo_block);
    }

    int ix = field.field_id("potential");

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
    if (enzo_block->is_leaf()) {
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

void EnzoMethodGravity::compute_accelerations (EnzoBlock * enzo_block) throw()
{

  // Extrapolate potential from "potential_curr" and "potential_prev"
  // if supercycling and not a solve step

  const bool is_root = enzo_block->index().is_root();

  if (is_supercycle_potential() && is_solve_cycle_(enzo_block)) {

    super_save_fields_(enzo_block);

  }

  if ( is_supercycle_potential() ) {

    if (is_solve_cycle_(enzo_block)) {

      super_update_time_(enzo_block, enzo_block->state()->time());

    } else {

      super_extrapolate_fields_(enzo_block, enzo_block->state()->time());

    }
  }


  Field field = enzo_block->data()->field();
  const int m = field.dimensions (0);

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (cosmology) {
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    enzo_float * potential = (enzo_float*) field.values ("potential");
    auto dt   = enzo_block->state()->dt();
    auto time = enzo_block->state()->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    for (int i=0; i<m; i++) potential[i] /= cosmo_a;
  }

  /// compute acceleration fields from potential

  EnzoComputeAcceleration compute_acceleration(cello::rank(), order_);

  if (is_supercycle_accelerations() ) {
    enzo_float * ax =      (enzo_float*) field.values ("acceleration_x");
    enzo_float * ay =      (enzo_float*) field.values ("acceleration_y");
    enzo_float * az =      (enzo_float*) field.values ("acceleration_z");
    enzo_float * ax_curr = (enzo_float*) field.values ("acceleration_x_curr");
    enzo_float * ay_curr = (enzo_float*) field.values ("acceleration_y_curr");
    enzo_float * az_curr = (enzo_float*) field.values ("acceleration_z_curr");

    if (is_solve_cycle_(enzo_block)) {

      super_shift_fields_(enzo_block);

      compute_acceleration.compute(enzo_block);

      super_save_fields_(enzo_block);

    } else { // not a solve step

      super_extrapolate_fields_(enzo_block, enzo_block->state()->time());
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

  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");

  const int rank = cello::rank();

  enzo_float dt = std::numeric_limits<enzo_float>::max();

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  double mean_cell_width;

  if (rank == 1) mean_cell_width = hx;
  if (rank == 2) mean_cell_width = sqrt(hx*hy);
  if (rank == 3) mean_cell_width = cbrt(hx*hy*hz);

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  if (cosmology) {
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt = block->state()->dt();
    double time = block->state()->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    mean_cell_width *= cosmo_a;
  }

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

  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
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
