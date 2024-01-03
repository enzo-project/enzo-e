// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-07
/// @brief    Implements the EnzoMethodGravity class

// #define TEST_FIELD_VALUES_TIME

#include "cello.hpp"
#include "enzo.hpp"

#include "charm_enzo.hpp"

// #define DEBUG_COPY_B
// #define DEBUG_COPY_DENSITIES
// #define DEBUG_COPY_POTENTIAL
// #define DEBUG_COPY_ACCELERATION

//----------------------------------------------------------------------

EnzoMethodGravity::EnzoMethodGravity
(int index_solver,
 double grav_const,
 int order,
 bool accumulate,
 int index_prolong,
 double dt_max)
  : Method(),
    index_solver_(index_solver),
    grav_const_(grav_const),
    order_(order),
    ir_exit_(-1),
    index_prolong_(index_prolong),
    dt_max_(dt_max)
{
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

#ifdef DEBUG_FIELD_FACE
  cello::define_field ("debug_1");
  cello::define_field ("debug_2");
#endif
#ifdef DEBUG_COPY_B
  cello::define_field ("B_copy");
#endif
#ifdef DEBUG_COPY_POTENTIAL
  cello::define_field ("potential_copy");
#endif
#ifdef DEBUG_COPY_DENSITIES
  cello::define_field ("density_total_copy");
#endif
#ifdef READ_ENZO_POTENTIAL
  cello::define_field ("potential_enzo");
  cello::define_field ("potential_dff");
#endif

  if (accumulate){
    cello::define_field ("density_particle");
    cello::define_field ("density_particle_accumulate");
  }

  // Refresh adds density_total field faces and one layer of ghost
  // zones to "B" field

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->set_prolong(index_prolong_);

  refresh->add_field("acceleration_x");
  refresh->add_field("acceleration_y");
  refresh->add_field("acceleration_z");
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
  refresh_exit->add_field("potential");

  refresh_exit->set_callback(CkIndex_EnzoBlock::p_method_gravity_end());
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{
  if (enzo::simulation()->state().cycle() == enzo::config()->initial_cycle) {
    // Check if the pm_deposit method is being used and precedes the
    // gravity method.
    ASSERT("EnzoMethodGravity",
           "Error: pm_deposit method must precede gravity method.",
           enzo::problem()->method_precedes("pm_deposit", "gravity"));
  }
  // Initialize the linear system

  Field field = block->data()->field();
  /// access problem-defining fields for eventual RHS and solution
  const int ib  = field.field_id ("B");
  const int id  = field.field_id("density");
  const int idt = field.field_id("density_total");
  const int idensity = (idt != -1) ? idt : id;
  ASSERT ("EnzoMethodGravity::compute",
          "modifying density in EnzoMethodGravity?",
          idensity != id);
  
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
#ifdef DEBUG_COPY_DENSITIES  
  const int id_copy = field.field_id ("D_copy");
  const int idt_copy = field.field_id ("DT_copy");
  enzo_float * D_copy = (enzo_float*) field.values (id_copy);

  enzo_float * DT_copy = (enzo_float*) field.values (idt_copy);
#endif  

  double time = block->state().time();
  enzo_float * D;
  // default current time
  D = (enzo_float*) field.values (idensity);
#ifdef TEST_FIELD_VALUES_TIME
  if (block->state().cycle() > 0) {
    // if history available use given time
    const double t0 = field.history_time(0);
    const double t1 = field.history_time(1);
    enzo_float * d0 = (enzo_float *) field.values("density",0);
    enzo_float * d1 = (enzo_float *) field.values("density",1);
    double a0e=-0.25;
    double ai=0.25;
    double a1e=1.25;
    const double ti = (ai*t0 + (1.0-ai)*t1);
    const double t0e = (a0e*t0 + (1.0-a0e)*t1);
    const double t1e = (a1e*t0 + (1.0-a1e)*t1);
    auto D0e = field.values_at<enzo_float>("density",t0e);
    auto D0 = field.values_at<enzo_float>("density",t0);
    auto DI = field.values_at<enzo_float>("density",ti);
    auto D1 = field.values_at<enzo_float>("density",t1);
    auto D1e = field.values_at<enzo_float>("density",t1e);
    double s0e=0, s0ec=0;
    double s0=0, s0c=0;
    double si=0, sic=0;
    double s1=0, s1c=0;
    double s1e=0, s1ec=0;
    for (int i=0; i<m; i++) {
      s0e += a0e*d0[i] + (1.0-a0e)*d1[i];
      s0ec += D0e[i];
      s0 += d0[i];
      s0c += D0[i];
      si += ai*d0[i] + (1.0-ai)*d1[i];
      sic += DI[i];
      s1 += d1[i];
      s1c += D1[i];
      s1e += a1e*d0[i] + (1.0-a1e)*d1[i];
      s1ec += D1e[i];
    }
    CkPrintf ("DEBUG_FIELD time %g %g\n",t0,t1);
    CkPrintf ("DEBUG_FIELD %s sum  %g %g\n",block->name().c_str(),s0,s1);
    CkPrintf ("DEBUG_FIELD %s D0e %20.16g %20.16g\n",
              block->name().c_str(),s0e,s0ec);
    CkPrintf ("DEBUG_FIELD %s D0 %20.16g %20.16g\n",
              block->name().c_str(),s0,s0c);
    CkPrintf ("DEBUG_FIELD %s DI %20.16g %20.16g\n",
              block->name().c_str(),si,sic);
    CkPrintf ("DEBUG_FIELD %s D1 %20.16g %20.16g\n",
              block->name().c_str(),s1,s1c);
    CkPrintf ("DEBUG_FIELD %s D1e %20.16g %20.16g\n",
              block->name().c_str(),s1e,s1ec);
    CkPrintf ("DEBUG_FIELD field %p %p\n",d0,d1);
  }
#endif

  double ds=0,bs=0;
  for (int i=0; i<m; i++) {
    D[i] += B[i];
    ds+=D[i];
    bs+=B[i];
  }

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
	    // In cosmological simulations, density units are defined such that `rho_bar_m` is
	    // 1.0, and time units are defined such that `4 * pi * G * rho_bar_m` is 1.0, where
	    // `G` is the gravitational constant, and `rho_bar_m` is the mean matter density
	    // of the universe. These choices of units result in Poisson's equation having a
	    // much simplified form.
	    D[i]=-(D[i]-1.0);
	    B[i]  = D[i];
	  }
	}
      }
    } else {
      const double scale = -4.0 * (cello::pi) * grav_const_;
      for (int iz=0; iz<mz; iz++) {
        for (int iy=0; iy<my; iy++) {
          for (int ix=0; ix<mx; ix++) {
            const int i = ix + mx*(iy + my*iz);
            B[i] = scale * D[i];
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

  const int ix = field.field_id ("potential");
  std::shared_ptr<Matrix> A (std::make_shared<EnzoMatrixLaplace>(order_));
  solver->set_field_x(ix);
  solver->set_field_b(ib);
#ifdef DEBUG_COPY_B
  if (B_copy) for (int i=0; i<m; i++) B_copy[i] = B[i];
#endif	
#ifdef DEBUG_COPY_DENSITIES
  enzo_float * DT = (enzo_float*) field.values (idt);
  if (DT_copy) for (int i=0; i<m; i++) DT_copy[i] = DT[i];
  if (D_copy) for (int i=0; i<m; i++) D_copy[i] = D[i];
#endif	
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
    double time = enzo_block->state().time();
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

#ifdef DEBUG_COPY_POTENTIAL
  enzo_float * potential_copy = (enzo_float*) field.values ("potential_copy");
  if (potential_copy) {
    for (int i=0; i<m; i++) potential_copy[i] = potential[i];
  }
#endif	
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
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

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
    double dt = block->state().dt();
    double time = block->state().time();
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
