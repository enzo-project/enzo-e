// See LICENSE_CELLO file for license and copyright information

/// @file     test_Particle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-10-12
/// @brief    Test program for the Particle class

#include "main.hpp"
#include "test.hpp"
#include <algorithm>

#include "data.hpp"

extern CProxy_Simulation proxy_simulation;
PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  proxy_simulation = CProxy_Simulation::ckNew("",0);
  
  unit_class("ParticleDescr");
  ParticleDescr * particle_descr = new ParticleDescr;
  particle_descr -> set_batch_size (1024);

  unit_class("ParticleData");
  ParticleData * particle_data = new ParticleData;

  Particle particle (particle_descr, particle_data);

  //--------------------------------------------------
  //   Type
  //--------------------------------------------------

  // new_type()
  // num_types()

  unit_func ("num_types()");
  unit_assert (particle.num_types() == 0);
  const int it_dark  = particle.new_type ("dark");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 1);
  unit_func ("num_types()");
  unit_assert (it_dark == 0);
  const int it_trace = particle.new_type ("trace");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 2);
  unit_func ("num_types()");
  unit_assert (it_trace == 1);
  const int it_star  = particle.new_type ("star");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 3);
  unit_func ("num_types()");
  unit_assert (it_star == 2);
  const int it_sink  = particle.new_type ("sink");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 4);
  unit_func ("num_types()");
  unit_assert (it_sink == 3);

  unit_func ("type_name()");
  unit_assert(particle.type_name(it_dark) == "dark");
  unit_assert(particle.type_name(it_trace) == "trace");
  unit_assert(particle.type_name(it_star) == "star");
  unit_assert(particle.type_name(it_sink) == "sink");

  unit_func ("type_index()");
  unit_assert(particle.type_index("dark")  == it_dark);
  unit_assert(particle.type_index("trace") == it_trace);
  unit_assert(particle.type_index("star")  == it_star);
  unit_assert(particle.type_index("sink")  == it_sink);

  
  //--------------------------------------------------
  // Interleaved
  //--------------------------------------------------

  unit_func("interleaved");
  unit_assert(particle.interleaved(it_dark) == false);
  unit_assert(particle.interleaved(it_trace) == false);

  unit_func("set_interleaved");
  particle.set_interleaved(it_dark,true);
  unit_assert(particle.interleaved(it_dark) == true);
  unit_assert(particle.interleaved(it_trace) == false);
  particle.set_interleaved(it_trace,true);
  unit_assert(particle.interleaved(it_dark) == true);
  unit_assert(particle.interleaved(it_trace) == true);
  particle.set_interleaved(it_dark,false);

  unit_assert(particle.interleaved(it_dark)  == false);
  unit_assert(particle.interleaved(it_trace) == true);

  //--------------------------------------------------
  //  Attribute
  //--------------------------------------------------

  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 0);

  unit_func ("new_attribute");
  const int ia_dark_x = particle.new_attribute (it_dark, "position_x", type_single);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 1);

  unit_func ("new_attribute");
  const int ia_dark_y = particle.new_attribute (it_dark, "position_y", type_single);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 2);

  unit_func ("new_attribute");
  const int ia_dark_z = particle.new_attribute (it_dark, "position_z", type_single);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 3);

  unit_func ("new_attribute");
  const int ia_dark_vx = particle.new_attribute (it_dark, "velocity_x", type_double);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 4);

  unit_func ("new_attribute");
  const int ia_dark_vy = particle.new_attribute (it_dark, "velocity_y", type_double);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 5);

  unit_func ("new_attribute");
  const int ia_dark_vz = particle.new_attribute (it_dark, "velocity_z", type_double);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 6);

  unit_func ("new_attribute");
  const int ia_dark_m = particle.new_attribute (it_dark, "mass", type_double);
  
  unit_func ("attribute_name()");
  unit_assert(particle.attribute_name(it_dark,ia_dark_x) == "position_x");
  unit_assert(particle.attribute_name(it_dark,ia_dark_y) == "position_y");
  unit_assert(particle.attribute_name(it_dark,ia_dark_z) == "position_z");
  unit_assert(particle.attribute_name(it_dark,ia_dark_vx) == "velocity_x");
  unit_assert(particle.attribute_name(it_dark,ia_dark_vy) == "velocity_y");
  unit_assert(particle.attribute_name(it_dark,ia_dark_vz) == "velocity_z");

  particle.set_position(it_dark,ia_dark_x, ia_dark_y, ia_dark_z);
  particle.set_velocity(it_dark,ia_dark_vx,ia_dark_vy,ia_dark_vz);

  unit_func ("attribute_index()");
  unit_assert(particle.attribute_index(it_dark,"position_x")  == ia_dark_x);

  const int ia_trace_x = particle.new_attribute 
    (it_trace, "position_x", type_int32);
  const int ia_trace_y = particle.new_attribute 
    (it_trace, "position_y", type_int64);
  const int ia_trace_z = particle.new_attribute 
    (it_trace, "position_z", type_int16);

  particle.set_position(it_trace,ia_trace_x, ia_trace_y, ia_trace_z);

  //--------------------------------------------------
  //  Constants
  //--------------------------------------------------

  unit_func("num_constants()");

  unit_assert (particle.num_constants(it_trace)==0);
  unit_assert (particle.num_constants(it_dark)==0);
  unit_assert (particle.num_constants(it_sink)==0);
  unit_assert (particle.num_constants(it_star)==0);

  unit_func("new_constant()");

  const int ic_dark_mass = particle.new_constant(it_dark, "mass", type_double);
  const int ic_star_mass = particle.new_constant(it_star, "mass", type_float);
  const int ic_star_type = particle.new_constant(it_star, "type", type_int8);

  unit_assert (particle.num_constants(it_trace)==0);
  unit_assert (particle.num_constants(it_dark)==1);
  unit_assert (particle.num_constants(it_sink)==0);
  unit_assert (particle.num_constants(it_star)==2);

  unit_func("constant_index()");
  unit_assert (particle.constant_index(it_dark,"mass") == ic_dark_mass);
  unit_assert (particle.constant_index(it_star,"mass") == ic_star_mass);
  unit_assert (particle.constant_index(it_star,"type") == ic_star_type);

  unit_func("constant_name()");
  unit_assert (particle.constant_name(it_dark,ic_dark_mass) == "mass");
  unit_assert (particle.constant_name(it_star,ic_star_mass) == "mass");
  unit_assert (particle.constant_name(it_star,ic_star_type) == "type");

  unit_func("constant_bytes()");
  unit_assert (particle.constant_bytes(it_dark,ic_dark_mass) == sizeof(double));
  unit_assert (particle.constant_bytes(it_star,ic_star_mass) == sizeof(float));
  unit_assert (particle.constant_bytes(it_star,ic_star_type) == 1);

  unit_func("constant_array()");
  char * dark_constant_array = particle.constant_array (it_dark);
  char * star_constant_array = particle.constant_array (it_star);
  double * dark_mass = 
    (double *) particle.constant_value  (it_dark,ic_dark_mass);
  float * star_mass = 
    (float *) particle.constant_value   (it_star,ic_star_mass);
  int8_t * star_type = 
    (int8_t *) particle.constant_value   (it_star,ic_star_type);
  *dark_mass = 10e12;
  *star_mass = 1.98855e33;
  *star_type = 'x';

  const int io_dark_mass = particle.constant_offset(it_dark,ic_dark_mass);
  const int io_star_mass = particle.constant_offset(it_star,ic_star_mass);
  const int io_star_type = particle.constant_offset(it_star,ic_star_type);

  double * dm = (double *)(&dark_constant_array[io_dark_mass]);
  float  * sm = (float *) (&star_constant_array[io_star_mass]);
  int8_t * st = (int8_t *)(&star_constant_array[io_star_type]);
  unit_assert (*dm == *dark_mass);
  unit_assert (*sm == *star_mass);
  unit_assert (*st == *star_type);
  
  //--------------------------------------------------
  //   Batch
  //--------------------------------------------------

  unit_assert (particle.batch_size() == 1024);

  unit_func("index()");
  int ib,ip;

  particle.index(0,&ib,&ip);
  unit_assert (ib==0 && ip==0);
  particle.index(1023,&ib,&ip);
  unit_assert (ib==0 && ip==1023);
  particle.index(1024,&ib,&ip);
  unit_assert (ib==1 && ip==0);
  particle.index(1024*1024,&ib,&ip);
  unit_assert (ib==1024 && ip==0);
  
  //--------------------------------------------------
  //   Type
  //--------------------------------------------------

  unit_assert(particle.attribute_type(it_dark,ia_dark_x) == type_single);
  unit_assert(particle.attribute_type(it_dark,ia_dark_y) == type_single);
  unit_assert(particle.attribute_type(it_dark,ia_dark_z) == type_single);
  unit_assert(particle.attribute_type(it_dark,ia_dark_vx) == type_double);
  unit_assert(particle.attribute_type(it_dark,ia_dark_vy) == type_double);
  unit_assert(particle.attribute_type(it_dark,ia_dark_vz) == type_double);
  unit_assert(particle.attribute_type(it_dark,ia_dark_m) == type_double);

  unit_assert(particle.attribute_type(it_trace,ia_trace_x) == type_int32);
  unit_assert(particle.attribute_type(it_trace,ia_trace_y) == type_int64);
  unit_assert(particle.attribute_type(it_trace,ia_trace_z) == type_int16);

  //--------------------------------------------------
  //   Bytes
  //--------------------------------------------------

  unit_func("attribute_bytes()");

  unit_assert(particle.attribute_bytes(it_dark,ia_dark_x) == 4);
  unit_assert(particle.attribute_bytes(it_dark,ia_dark_y) == 4);
  unit_assert(particle.attribute_bytes(it_dark,ia_dark_z) == 4);
  unit_assert(particle.attribute_bytes(it_dark,ia_dark_vx) == 8);
  unit_assert(particle.attribute_bytes(it_dark,ia_dark_vy) == 8);
  unit_assert(particle.attribute_bytes(it_dark,ia_dark_vz) == 8);
  unit_assert(particle.attribute_bytes(it_dark,ia_dark_m) == 8);

  unit_assert(particle.particle_bytes(it_dark) == 4+4+4+8+8+8+8);

  unit_assert(particle.attribute_bytes(it_trace,ia_trace_x) == 4);
  unit_assert(particle.attribute_bytes(it_trace,ia_trace_y) == 8);
  unit_assert(particle.attribute_bytes(it_trace,ia_trace_z) == 2);

  // particles have two bytes padding so 2,4,8 all divide evenly
  // to compute stride.  
  
  unit_assert(particle.particle_bytes(it_trace) == 16);

  unit_func("attribute_offset()");

  
  int mp = particle.batch_size();

  // not interleaved
  unit_assert(particle.attribute_offset(it_dark,ia_dark_x) == 0*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_dark_y) == 4*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_dark_z) == 8*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_dark_vx) == 12*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_dark_vy) == 20*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_dark_vz) == 28*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_dark_m) == 36*mp);

  // interleaved
  unit_assert(particle.attribute_bytes(it_trace,ia_trace_x) == 4);
  unit_assert(particle.attribute_bytes(it_trace,ia_trace_y) == 8);
  unit_assert(particle.attribute_bytes(it_trace,ia_trace_z) == 2);


  unit_func("stride");

  // not interleaved
  unit_assert (particle.stride(it_dark,ia_dark_x) == 1);
  unit_assert (particle.stride(it_dark,ia_dark_y) == 1);
  unit_assert (particle.stride(it_dark,ia_dark_z) == 1);
  unit_assert (particle.stride(it_dark,ia_dark_vx) == 1);
  unit_assert (particle.stride(it_dark,ia_dark_vy) == 1);
  unit_assert (particle.stride(it_dark,ia_dark_vz) == 1);

  // interleaved
  unit_assert (particle.stride(it_trace,ia_trace_x) == 16/4);
  unit_assert (particle.stride(it_trace,ia_trace_y) == 16/8);
  unit_assert (particle.stride(it_trace,ia_trace_z) == 16/2);

  //--------------------------------------------------
  //   Insert
  //--------------------------------------------------

  unit_func("num_batches");

  unit_assert(particle.num_batches(it_dark) == 0);
  unit_assert(particle.num_batches(it_trace) == 0);

  unit_func("insert_particles()");

  int i_dark_0 = particle.insert_particles (it_dark, 10000);
  int i_trace_0 = particle.insert_particles (it_trace, 20000);

  unit_assert (i_dark_0 == 0);
  unit_assert (i_trace_0 == 0);

  unit_func("num_particles()");
  unit_assert (particle.num_particles(it_dark) == 10000);
  unit_func("num_batches()");
  unit_assert (particle.num_batches(it_dark) == 10000/1024 + 1);

  unit_func("num_particles()");
  unit_assert (particle.num_particles(it_trace) == 20000);
  unit_func("num_batches()");
  unit_assert (particle.num_batches(it_trace) == 20000/1024 + 1);

  unit_func("index()");
  particle.index(10000,&ib,&ip);
  unit_assert (ib == 10000/1024);
  unit_assert (ip == 10000%1024);

  unit_func("index()");
  particle.index(20000,&ib,&ip);
  unit_assert (ib == 20000/1024);
  unit_assert (ip == 20000%1024);
  
  int i_dark_1 = particle.insert_particles (it_dark, 20000);
  int i_trace_1 = particle.insert_particles (it_trace, 10000);

  unit_assert (i_dark_1 == 10000);
  unit_assert (i_trace_1 == 20000);

  unit_func("num_particles()");
  unit_assert (particle.num_particles(it_dark) == 30000);
  unit_func("num_batches()");
  unit_assert (particle.num_batches(it_dark) == 30000/1024 + 1);

  unit_func("num_particles()");
  unit_assert (particle.num_particles(it_trace) == 30000);
  unit_func("num_batches()");
  unit_assert (particle.num_batches(it_trace) == 30000/1024 + 1);

  unit_func("index()");
  particle.index(30000,&ib,&ip);
  unit_assert (ib == 30000/1024);
  unit_assert (ip == 30000%1024);

  //--------------------------------------------------
  // dark assignment and compare
  //--------------------------------------------------

  int nb = particle.num_batches(it_dark);
  int count_particles = 0;
  int index;

  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(it_dark,ib);
    float  *  x = (float  *) particle.attribute_array(it_dark,ia_dark_x, ib);
    float  *  y = (float  *) particle.attribute_array(it_dark,ia_dark_y, ib);
    float  *  z = (float  *) particle.attribute_array(it_dark,ia_dark_z, ib);
    double * vx = (double *) particle.attribute_array(it_dark,ia_dark_vx,ib);
    double * vy = (double *) particle.attribute_array(it_dark,ia_dark_vy,ib);
    double * vz = (double *) particle.attribute_array(it_dark,ia_dark_vz,ib);
    int dx = particle.stride(it_dark,ia_dark_x);
    int dv = particle.stride(it_dark,ia_dark_vx);
    for (int ip=0,ix=0,iv=0; ip<np; ip++,ix+=dx,iv+=dv) {
      index = ip + ib*mp;
      count_particles ++;
      x[ix]  = 10*index;
      y[ix]  = 10*index+1;
      z[ix]  = 10*index+2;
      vx[iv] = 10*index+3;
      vy[iv] = 10*index+4;
      vz[iv] = 10*index+5;
    }
  }
  unit_assert(count_particles == 30000);

  // test position() and velocity()
  double xp[mp], yp[mp], zp[mp];
  double vxp[mp],vyp[mp],vzp[mp];
  int error_position=0;
  int error_velocity=0;
  for (int ib=0; ib<nb; ib++) {
    unit_func("position()");
    unit_assert(particle.position(it_dark, ib, xp, yp, zp));
    unit_func("velocity()");
    unit_assert(particle.velocity(it_dark, ib, vxp,vyp,vzp));
    const int np = particle.num_particles(it_dark,ib);
    for (int ip=0; ip<np; ip++) {
      index = ip + ib*mp;
      if (xp[ip] != 10*index) error_position++;
      if (yp[ip] != 10*index+1) error_position++;
      if (zp[ip] != 10*index+2) error_position++;
      if (vxp[ip] != 10*index+3) error_velocity++;
      if (vyp[ip] != 10*index+4) error_velocity++;
      if (vzp[ip] != 10*index+5) error_velocity++;
    }
  }
  unit_func("position()");
  unit_assert(error_position == 0);
  unit_func("velocity()");
  unit_assert(error_velocity == 0);

  // run through again and compare values before deleting
  nb = particle.num_batches(it_dark);
  int count_wrong[6];
  for (int i=0; i<6; i++) count_wrong[i] = 0;
  for (int ib=0; ib<nb; ib++) {
    index = ib*mp;
    int np = particle.num_particles(it_dark,ib);
    float  * x =  (float  *) particle.attribute_array(it_dark,ia_dark_x, ib);
    float  * y =  (float  *) particle.attribute_array(it_dark,ia_dark_y, ib);
    float  * z =  (float  *) particle.attribute_array(it_dark,ia_dark_z, ib);
    double * vx = (double *) particle.attribute_array(it_dark,ia_dark_vx,ib);
    double * vy = (double *) particle.attribute_array(it_dark,ia_dark_vy,ib);
    double * vz = (double *) particle.attribute_array(it_dark,ia_dark_vz,ib);
    int dx = particle.stride(it_dark,ia_dark_x);
    int dv = particle.stride(it_dark,ia_dark_vx);
    for (int ip=0,ix=0,iv=0; ip<np; ip++,ix+=dx,iv+=dv) {
      index = ip + ib*mp;
      if (x[ix]  != 10*index ) count_wrong[0]++;
      if (y[ix]  != 10*index+1) count_wrong[1]++;
      if (z[ix]  != 10*index+2) count_wrong[2]++;
      if (vx[iv] != 10*index+3) count_wrong[3]++;
      if (vy[iv] != 10*index+4) count_wrong[4]++;
      if (vz[iv] != 10*index+5) count_wrong[5]++;
      index++;
    }
  }
  unit_assert (count_wrong[0] == 0);
  unit_assert (count_wrong[1] == 0);
  unit_assert (count_wrong[2] == 0);
  unit_assert (count_wrong[3] == 0);
  unit_assert (count_wrong[4] == 0);
  unit_assert (count_wrong[5] == 0);

  //--------------------------------------------------
  // trace assignment and compare
  //--------------------------------------------------

  nb = particle.num_batches(it_trace);
  count_particles = 0;
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(it_trace,ib);
    int32_t * x = (int32_t *) particle.attribute_array(it_trace,ia_trace_x,ib);
    int64_t * y = (int64_t *) particle.attribute_array(it_trace,ia_trace_y,ib);
    int16_t * z = (int16_t *) particle.attribute_array(it_trace,ia_trace_z,ib);
    int dx=particle.stride(it_trace,ia_trace_x);
    int dy=particle.stride(it_trace,ia_trace_y);
    int dz=particle.stride(it_trace,ia_trace_z);
    for (int ip=0,ix=0,iy=0,iz=0; ip<np; ip++,ix+=dx,iy+=dy,iz+=dz) {
      index = ip + ib*mp;
      count_particles ++;
      x[ix]  = 3*index;
      y[iy]  = 4*index+1;
      z[iz]  = 1*index+2;
    }
  }
  unit_assert(count_particles == 30000);

  // test position() and velocity()
  error_position=0;
  for (int ib=0; ib<nb; ib++) {
    unit_func("position()");
    unit_assert(particle.position(it_trace, ib, xp, yp, zp));
    unit_func("velocity()");
    unit_assert(!particle.velocity(it_trace, ib, vxp, vyp, vzp));
    const int np = particle.num_particles(it_trace,ib);
    for (int ip=0; ip<np; ip++) {
      index = ip + ib*mp;
      if (xp[ip] != 3*index) error_position++;
      if (yp[ip] != 4*index+1) error_position++;
      if (zp[ip] != 1*index+2) error_position++;
    }
  }
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // These position tests fail because they are represented as
  // integers, which are assumed to be relative to the Block.  These
  // are converted to float [-2,2) by Particle::position(), where the
  // Block bounds are [-1,1).
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  unit_func("position()");
  unit_assert(error_position == 0);

  // run through again and compare values before deleting
  nb = particle.num_batches(it_trace);
  for (int i=0; i<3; i++) count_wrong[i] = 0;
  for (int ib=0; ib<nb; ib++) {
    index = ib*mp;
    int np = particle.num_particles(it_trace,ib);
    int32_t * x = (int32_t *) particle.attribute_array(it_trace,ia_trace_x,ib);
    int64_t * y = (int64_t *) particle.attribute_array(it_trace,ia_trace_y,ib);
    int16_t * z = (int16_t *) particle.attribute_array(it_trace,ia_trace_z,ib);
    int dx=particle.stride(it_trace,ia_trace_x);
    int dy=particle.stride(it_trace,ia_trace_y);
    int dz=particle.stride(it_trace,ia_trace_z);
    for (int ip=0,ix=0,iy=0,iz=0; ip<np; ip++,ix+=dx,iy+=dy,iz+=dz) {
      index = ip + ib*mp;
      if (x[ix]  != 3*index ) count_wrong[0]++;
      if (y[iy]  != 4*index+1) count_wrong[1]++;
      if (z[iz]  != 1*index+2) count_wrong[2]++;
      index++;
    }
  }
  unit_assert (count_wrong[0] == 0);
  unit_assert (count_wrong[1] == 0);
  unit_assert (count_wrong[2] == 0);
  unit_assert (count_wrong[3] == 0);
  unit_assert (count_wrong[4] == 0);
  unit_assert (count_wrong[5] == 0);


  //--------------------------------------------------
  //   Delete
  //--------------------------------------------------


  // initialize mask for particles to delete
  bool mask[1024];
  int count_delete;

  for (int i=0; i<1024; i++) {
    mask[i] = (i % 3 == 0);
  }

  unit_func("delete_particles()");

  count_delete = 0;
  nb = particle.num_batches(it_dark);
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(it_dark,ib);
    for (int ip=0; ip<np; ip++) if (mask[ip]) count_delete++;
    particle.delete_particles(it_dark,ib,mask);
  }
  unit_assert(particle.num_particles(it_dark) == 30000 - count_delete);

  count_delete = 0;
  nb = particle.num_batches(it_trace);
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(it_trace,ib);
    for (int ip=0; ip<np; ip++) if (mask[ip]) count_delete++;
    particle.delete_particles(it_trace,ib,mask);
  }
  
  unit_assert(particle.num_particles(it_trace) == 30000 - count_delete);

  //--------------------------------------------------
  // dark particle values after delete
  //--------------------------------------------------
  index = 0;
  nb = particle.num_batches(it_dark);
  for (int i=0; i<6; i++) count_wrong[i] = 0;
  count_particles = 0;
  for (int ib=0; ib<nb; ib++) {
    index = ib*mp;
    int np = particle.num_particles(it_dark,ib);
    float  * x  = (float  *) particle.attribute_array(it_dark,ia_dark_x, ib);
    float  * y  = (float  *) particle.attribute_array(it_dark,ia_dark_y, ib);
    float  * z  = (float  *) particle.attribute_array(it_dark,ia_dark_z, ib);
    double * vx = (double *) particle.attribute_array(it_dark,ia_dark_vx,ib);
    double * vy = (double *) particle.attribute_array(it_dark,ia_dark_vy,ib);
    double * vz = (double *) particle.attribute_array(it_dark,ia_dark_vz,ib);
    int dx = particle.stride(it_dark,ia_dark_x);
    int dv = particle.stride(it_dark,ia_dark_vx);
    for (int ip=0,ix=0,iv=0; ip<np; ip++,ix+=dx,iv+=dv) {
      count_particles ++;
      if (ip % 2 == 0) index++;
      if (x[ix]  != 10*index ) count_wrong[0]++;
      if (y[ix]  != 10*index+1) count_wrong[1]++;
      if (z[ix]  != 10*index+2) count_wrong[2]++;
      if (vx[iv] != 10*index+3) count_wrong[3]++;
      if (vy[iv] != 10*index+4) count_wrong[4]++;
      if (vz[iv] != 10*index+5) count_wrong[5]++;
      
      index++;
    }
  }
  unit_assert (count_wrong[0] == 0);
  unit_assert (count_wrong[1] == 0);
  unit_assert (count_wrong[2] == 0);
  unit_assert (count_wrong[3] == 0);
  unit_assert (count_wrong[4] == 0);
  unit_assert (count_wrong[5] == 0);

  unit_assert (30000 - count_delete == count_particles);

  //--------------------------------------------------
  // trace particle values after delete
  //--------------------------------------------------

  index = 0;
  nb = particle.num_batches(it_trace);
  for (int i=0; i<6; i++) count_wrong[i] = 0;
  count_particles = 0;
  for (int ib=0; ib<nb; ib++) {
    index = ib*mp;
    int np = particle.num_particles(it_trace,ib);
    int32_t * x = (int32_t *) particle.attribute_array(it_trace,ia_trace_x,ib);
    int64_t * y = (int64_t *) particle.attribute_array(it_trace,ia_trace_y,ib);
    int16_t * z = (int16_t *) particle.attribute_array(it_trace,ia_trace_z,ib);
    int dx=particle.stride(it_trace,ia_trace_x);
    int dy=particle.stride(it_trace,ia_trace_y);
    int dz=particle.stride(it_trace,ia_trace_z);
    for (int ip=0,ix=0,iy=0,iz=0; ip<np; ip++,ix+=dx,iy+=dy,iz+=dz) {
      count_particles ++;
      if (ip % 2 == 0) index++;
      if (x[ix]  != 3*index ) count_wrong[0]++;
      if (y[iy]  != 4*index+1) count_wrong[1]++;
      if (z[iz]  != 1*index+2) count_wrong[2]++;
      
      index++;
    }
  }
  unit_assert (count_wrong[0] == 0);
  unit_assert (count_wrong[1] == 0);
  unit_assert (count_wrong[2] == 0);

  unit_assert (30000 - count_delete == count_particles);


  //--------------------------------------------------
  // (RE)INSERT PARTICLES
  //--------------------------------------------------

  unit_func("(re)insert_particles()");

  int ip0,i2,np;

  // dark (re)insert and re-test

  nb = particle.num_batches(it_dark);
  ip0 = particle.num_particles(it_dark,nb-1);
  i2 = particle.insert_particles (it_dark, 10000);
  unit_assert (i2 == 1024*(nb-1) + ip0);

  unit_assert (particle.num_particles() == 10000 + 2*count_particles);

  unit_func("num_particles()");
  nb = particle.num_batches(it_dark);
  np=particle.num_particles(it_dark);
  int count = 0;
  for (int ib=0; ib<nb; ib++) {
    count += particle.num_particles(it_dark,ib);
  }
  unit_assert (np==count);

  // trace (re)insert 

  nb = particle.num_batches(it_trace);
  ip0 = particle.num_particles(it_trace,nb-1);
  i2 = particle.insert_particles (it_trace, 10000);
  unit_assert (i2 == 1024*(nb-1) + ip0);
  unit_assert (particle.num_particles() == 20000 + 2*count_particles);

  unit_func("num_particles()");
  nb = particle.num_batches(it_trace);
  np=particle.num_particles(it_trace);
  count = 0;
  for (int ib=0; ib<nb; ib++) {
    count += particle.num_particles(it_trace,ib);
  }
  unit_assert (np==count);

  //======================================================================

  unit_func("split_particles()");

  // ParticleData * particle_data_2 = new ParticleData;
  // Particle particle_2 (particle_descr,particle_data_2);

  count_particles = particle.num_particles(it_dark);

  count_delete = 0;
  for (int i=0; i<1024; i++) {
    mask[i] = (i % 5 == 0);
  }

  unit_func("split_particles()");
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(it_dark,ib);
    for (int ip=0; ip<np; ip++) if (mask[ip]) count_delete++;
    particle.delete_particles(it_dark,ib,mask);
  }

  unit_assert (count_particles - count_delete
	       == particle.num_particles(it_dark));


  unit_func("compress()");

  unit_assert (particle.efficiency (it_dark,0)  < 0.55);
  unit_assert (particle.efficiency (it_dark)    < 0.65);
  unit_assert (particle.efficiency (it_trace,0) < 0.70);
  unit_assert (particle.efficiency (it_trace)   < 0.80);
  unit_assert (particle.efficiency ()           < 0.65);

  particle.compress(it_dark);

  unit_assert (particle.efficiency (it_dark,0)  > 0.99);
  unit_assert (particle.efficiency (it_dark)    > 0.85);
  unit_assert (particle.efficiency (it_trace,0) < 0.70);
  unit_assert (particle.efficiency (it_trace)   < 0.80);
  unit_assert (particle.efficiency ()           > 0.85);

  particle.compress(it_trace);

  unit_assert (particle.efficiency (it_dark,0)  > 0.99);
  unit_assert (particle.efficiency (it_dark)    > 0.85);
  unit_assert (particle.efficiency (it_trace,0) > 0.99);
  unit_assert (particle.efficiency (it_trace)   > 0.99);
  unit_assert (particle.efficiency ()           > 0.90);

  //--------------------------------------------------
  //   GATHER / SCATTER
  //--------------------------------------------------

  ParticleData pd_src,*pd_array[16],pd_dst;

  Particle p_src(particle_descr,&pd_src);
  Particle * p_array[16];

  // Test gather and scatter with simulated particle exchange with
  // neighbors.  Use 4x4 array with some NULL and some duplicated
  // (for coarse-level and same-level neighbors).  
  //
  //    x---x---x---x---x       x---x---x
  //    | 4 | 5 | 6 | 7 |     4 | 5 | 6 | 7
  //    x---x---x---x---x    ---x---x---x---
  //    | 2 |       | 3 |       |       |
  //    x---x       x---x     2 |       | 3
  //    | 2 |       | 3 |       |       |
  //    x---x---x---x---x    ---X-------X---
  //    | 0 | 1 | 1 | 1 |     0 |   1  
  //    x---x---X---x---x       |

  pd_array [0] = new ParticleData;
  pd_array [1] = new ParticleData;
  pd_array [2] = pd_array[1];
  pd_array [3] = pd_array[2];

  pd_array [4] = new ParticleData;
  pd_array [5] = NULL;
  pd_array [6] = NULL;
  pd_array [7] = new ParticleData;

  pd_array [8] = pd_array[4];
  pd_array [9] = NULL;
  pd_array[10] = NULL;
  pd_array[11] = pd_array[7];

  pd_array[12] = new ParticleData;
  pd_array[13] = new ParticleData;
  pd_array[14] = new ParticleData;
  pd_array[15] = new ParticleData;

  for (int i=0; i<16; i++) {
    // Note some different Particle's with same ParticleData
    p_array[i] = (pd_array[i]) ? 
      new Particle(particle_descr,pd_array[i]) : NULL;
  }

  // initialize particles on nx x ny
  // send those outside ghost zone depth ng to neighbors
  //
  // e.g. N = 11x17
  // 66677777777777888
  // 66677777777777888
  // 66677777777777888
  // 333           555 
  // 333           555
  // 333           555
  // 333           555
  // 333           555
  // 000 111111111 222
  // 000 111111111 222
  // 000 111111111 222

  const int ndx = 166;
  const int ndy = 214;
  const int gmx = 2;
  const int gpx = 4;
  const int gmy = 3;
  const int gpy = 5;

  const int nx = ndx - (gmx + gpx);
  const int ny = ndy - (gmy + gpy);

  unit_assert ((nx/2)*2 == nx);
  unit_assert ((ny/2)*2 == ny);

  // insert uninitialized particles
  unit_func("insert()");
  p_src.insert_particles(it_dark,ndx*ndy);
  unit_assert(p_src.num_particles() == ndx*ndy);
  // set positions to (ix,iy)
  for (int ix=0; ix<ndx; ix++) {
    for (int iy=0; iy<ndy; iy++) {
      int ib,ip;
      p_src.index(ix+ndx*iy,&ib,&ip);
      float * x = (float *) p_src.attribute_array(it_dark,ia_dark_x,ib);
      float * y = (float *) p_src.attribute_array(it_dark,ia_dark_y,ib);
      x[ip] = (ix+1);
      y[ip] = (iy+1);
    }
  }

  //--------------------------------------------------
  // scatter to p_array[]
  //--------------------------------------------------

  int dx = particle.stride(it_dark,ia_dark_x);

  nb = p_src.num_batches(it_dark);

  mp = particle.batch_size();
  int index_array[mp];

  for (int ib=0; ib<nb; ib++) {
    np = p_src.num_particles(it_dark,ib);
    float * xa = (float *) p_src.attribute_array(it_dark,ia_dark_x,ib);
    float * ya = (float *) p_src.attribute_array(it_dark,ia_dark_y,ib);
    for (int ip=0; ip<mp; ip++) mask[ip] = false;
    for (int ip=0; ip<np; ip++) {
      float x = xa[ip*dx]-1;
      float y = ya[ip*dx]-1;
      // move particle only if in ghost region
      //
      //   gmx-0.5*nx  
      // < gmx+0.0*nx  0
      // < gmx+0.5*nx  1
      // < gmx+1.0*nx  2
      // < gmx+1.5*nx  3
      //
      // +---+---+---+---+
      //    gmx    gmx+nx
 
      int kx = (x-gmx+nx/2) / (0.5*nx);
      int ky = (y-gmy+ny/2) / (0.5*ny);
      mask[ip] = (kx==0 || kx==3 || ky==0 || ky==3);
      index_array[ip] = kx + 4*ky;

      // mask[ip]        = (kx || ky);
    }

    p_src.scatter(it_dark,ib,np,mask,index_array,16,pd_array);

  }

  unit_func ("scatter()");

  // as a first check, make sure number of particles moved is correct

  int ca00 = gmx*gmy;
  int ca10 = nx*gmy/2;
  int ca20 = nx*gmy/2;
  int ca30 = gpx*gmy;

  int ca01 = gmx*ny/2;
  int ca31 = gpx*ny/2;

  int ca02 = gmx*ny/2;
  int ca32 = gpx*ny/2;

  int ca03 = gmx*gpy;
  int ca13 = nx*gpy/2;
  int ca23 = nx*gpy/2;
  int ca33 = gpx*gpy;

  unit_assert(p_array[0]->num_particles(it_dark) == ca00);
  unit_assert(p_array[1]->num_particles(it_dark) == ca10+ca20+ca30);
  unit_assert(p_array[2]->num_particles(it_dark) == ca10+ca20+ca30);
  unit_assert(p_array[3]->num_particles(it_dark) == ca10+ca20+ca30);
  unit_assert(p_array[4]->num_particles(it_dark) == ca01+ca02);
  unit_assert(p_array[7]->num_particles(it_dark) == ca31+ca32);
  unit_assert(p_array[8]->num_particles(it_dark) == ca01+ca02);
  unit_assert(p_array[11]->num_particles(it_dark)== ca31+ca32);
  unit_assert(p_array[12]->num_particles(it_dark)== ca03);
  unit_assert(p_array[13]->num_particles(it_dark)== ca13);
  unit_assert(p_array[14]->num_particles(it_dark)== ca23);
  unit_assert(p_array[15]->num_particles(it_dark)== ca33);

  int count_total = 0;
  bool mask_scatter[ndx*ndy] = {0};
  int error_scatter_range=0;
  int error_scatter_duple=0;
  ParticleData * pd_array_sorted[16];
  for (int k=0; k<16; k++) pd_array_sorted[k] = pd_array[k];
  std::sort(&pd_array_sorted[0],
	    &pd_array_sorted[16]);
  for (int k=0; k<16; k++) {
    if (k>=0 && pd_array_sorted[k] == pd_array_sorted[k-1]) continue;

    if (pd_array_sorted[k]==NULL) continue;
    Particle  p(particle_descr,pd_array_sorted[k]);

    int nb = pd_array_sorted[k] ? p.num_batches(it_dark) : 0;
    for (int ib=0; ib<nb; ib++) {
      float * xa = (float *) p.attribute_array(it_dark,ia_dark_x,ib);
      float * ya = (float *) p.attribute_array(it_dark,ia_dark_y,ib);
      int np = p.num_particles(it_dark,ib);
      for (int ip=0; ip<np; ip++) {
	float x = xa[ip*dx]-1;
	float y = ya[ip*dx]-1;
	bool in_range = (0 <= x && x < ndx && 0 <= y && y < ndy);
	int ix = (int)x;
	int iy = (int)y;
	if (! in_range) {
	  ++error_scatter_range;
	} else {// in_range
	  if (mask_scatter [ix + ndx*iy]) { // duplelicate particle
	    ++ error_scatter_duple;
	  } else {
	    mask_scatter [ix+ndx*iy] = true;
	  }
	}
	++count_total;
      }
    }
  }

  unit_assert (error_scatter_range == 0);
  unit_assert (error_scatter_duple == 0);
  int error_scatter_int = 0;

  for (int iy=0; iy<ndy; iy++) {
    for (int ix=0; ix<ndx; ix++) {
      bool in_range = (0 <= ix && ix < ndx && 0 <= iy && iy < ndy);
      bool in_interior = (gmx <= ix && ix < ndx-gpx &&
			  gmy <= iy && iy < ndy-gpy);
      // all particles are not interior zones 
      if (! (mask_scatter[ix+ndx*iy] == 
	     (in_range && (! in_interior)))) {
	error_scatter_int++;
      }
      //	printf ("%d",mask_scatter[ix+ndx*iy]?1:0);
    }
    //     printf ("\n");
  }

  unit_assert (error_scatter_int == 0);

  //--------------------------------------------------
  // gather to p_dst[]
  //--------------------------------------------------

  unit_func ("gather()");

  Particle p_dst(particle_descr,&pd_dst);

  p_dst.gather(it_dark,16,pd_array);

  
  // printf ("%d %d\n",p_dst.num_particles(it_dark) , count_total);

  // FIRST FAIL
  unit_assert (p_dst.num_particles(it_dark) == count_total);

  // ensure all ghost particles are in p_dst, and only ghost particles
  bool mask_array[ndx*ndy] = {0};
  nb = p_dst.num_batches(it_dark);
  int error_gather_range = 0;
  int error_gather_duple = 0;
  for (int ib=0; ib<nb; ib++) {
    int np = p_dst.num_particles(it_dark,ib);
    float * xa = (float *) p_dst.attribute_array(it_dark,ia_dark_x,ib);
    float * ya = (float *) p_dst.attribute_array(it_dark,ia_dark_y,ib);
    for (int ip=0; ip<np; ip++) {
      float x = xa[ip*dx]-1;
      float y = ya[ip*dx]-1;
      bool in_range = (0 <= x && x < ndx && 0 <= y && y < ndy);
      int ix = (int)x;
      int iy = (int)y;
      if (! in_range) {
	printf ("ERROR gather range %d %d\n",ix,iy);
	++error_gather_range;
      } else {// in_range
	if (mask_array [ix + ndx*iy]) { // duplicate particle
	  printf ("ERROR gather duple %d %d\n",ix,iy);
	  ++ error_gather_duple;
	} else {
	  mask_array [ix+ndx*iy] = true;
	}
      }
    }
  }

  // printf ("error_gather_range %d\n",error_gather_range);
  // printf ("error_gather_duple %d\n",error_gather_duple);
  unit_assert (error_gather_range == 0);
  unit_assert (error_gather_duple == 0);
  int error_gather_int = 0;

  for (int iy=0; iy<ndy; iy++) {
    for (int ix=0; ix<ndx; ix++) {
      bool in_range = (0 <= ix && ix < ndx && 0 <= iy && iy < ndy);
      bool in_interior = (gmx <= ix && ix < ndx-gpx &&
			  gmy <= iy && iy < ndy-gpy);
      // all particles are not interior zones 
      if (! (mask_array[ix+ndx*iy] == 
	     (in_range && (! in_interior)))) {
	error_gather_int++;
      }
	// printf ("%d",mask_array[ix+ndx*iy]?1:0);
    }
     // printf ("\n");
  }

  unit_assert (error_gather_int == 0);

  //--------------------------------------------------
  // data_size(), save_data(), load_data() 
  //--------------------------------------------------

  unit_func("data_size()");
  int n = p_dst.data_size();

  unit_func("save_data()");
  char * buffer = new char[n];
  char * buffer_next = p_dst.save_data(buffer);
  unit_assert (buffer_next - buffer == n);
  if (buffer_next - buffer != n)
    printf ("buffer size mismatch: %d %d\n",buffer_next - buffer,n);

  unit_func("load_data()");
  ParticleData new_p_data;
  Particle new_p (particle_descr,&new_p_data);
  buffer_next = new_p.load_data(buffer);
  if (buffer_next - buffer != n)
    printf ("buffer size mismatch: %d %d\n",buffer_next - buffer,n);
  unit_assert (buffer_next - buffer == n);
  unit_assert (p_dst == new_p);


  // printf ("error_gather_int %d\n",error_gather_int);

  //--------------------------------------------------
  //   Grouping
  //--------------------------------------------------

  // See test_Grouping.cpp

  unit_finalize();

  exit_();


}


PARALLEL_MAIN_END

