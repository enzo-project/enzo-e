// See LICENSE_CELLO file for license and copyright information

/// @file     test_Particle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-10-12
/// @brief    Test program for the Particle class

#include "main.hpp"
#include "test.hpp"

#include "data.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ParticleDescr");
  ParticleDescr * particle_descr = new ParticleDescr(1024);

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
  const int i_star  = particle.new_type ("star");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 3);
  unit_func ("num_types()");
  unit_assert (i_star == 2);
  const int i_sink  = particle.new_type ("sink");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 4);
  unit_func ("num_types()");
  unit_assert (i_sink == 3);

  unit_func ("type_name()");
  unit_assert(particle.type_name(it_dark) == "dark");
  unit_assert(particle.type_name(it_trace) == "trace");
  unit_assert(particle.type_name(i_star) == "star");
  unit_assert(particle.type_name(i_sink) == "sink");

  unit_func ("type_index()");
  unit_assert(particle.type_index("dark")  == it_dark);
  unit_assert(particle.type_index("trace") == it_trace);
  unit_assert(particle.type_index("star")  == i_star);
  unit_assert(particle.type_index("sink")  == i_sink);

  
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

  unit_assert(particle.interleaved(it_dark) == false);

  unit_assert(particle.interleaved(it_trace) == true);

  //--------------------------------------------------
  //  Attribute
  //--------------------------------------------------

  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 0);
  unit_func ("new_attribute");
  const int ia_x = particle.new_attribute (it_dark, "position_x", 4);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 1);
  unit_func ("new_attribute");
  const int ia_y = particle.new_attribute (it_dark, "position_y", 4);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 2);
  unit_func ("new_attribute");
  const int ia_z = particle.new_attribute (it_dark, "position_z", 4);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 3);
  unit_func ("new_attribute");
  const int ia_vx = particle.new_attribute (it_dark, "velocity_x", 8);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 4);
  unit_func ("new_attribute");
  const int ia_vy = particle.new_attribute (it_dark, "velocity_y", 8);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 5);
  unit_func ("new_attribute");
  const int ia_vz = particle.new_attribute (it_dark, "velocity_z", 8);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(it_dark) == 6);
  unit_func ("new_attribute");
  const int ia_m = particle.new_attribute (it_dark, "mass", 8);
  
  unit_func ("attribute_name()");
  unit_assert(particle.attribute_name(it_dark,ia_x) == "position_x");
  unit_assert(particle.attribute_name(it_dark,ia_y) == "position_y");
  unit_assert(particle.attribute_name(it_dark,ia_z) == "position_z");
  unit_assert(particle.attribute_name(it_dark,ia_vx) == "velocity_x");
  unit_assert(particle.attribute_name(it_dark,ia_vy) == "velocity_y");
  unit_assert(particle.attribute_name(it_dark,ia_vz) == "velocity_z");

  unit_func ("attribute_index()");
  unit_assert(particle.attribute_index(it_dark,"position_x")  == ia_x);

  const int it_x = particle.new_attribute (it_trace, "position_x", 4);
  const int it_y = particle.new_attribute (it_trace, "position_y", 8);
  const int it_z = particle.new_attribute (it_trace, "position_z", 2);

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
  //   Bytes
  //--------------------------------------------------

  unit_func("attribute_bytes()");

  unit_assert(particle.attribute_bytes(it_dark,ia_x) == 4);
  unit_assert(particle.attribute_bytes(it_dark,ia_y) == 4);
  unit_assert(particle.attribute_bytes(it_dark,ia_z) == 4);
  unit_assert(particle.attribute_bytes(it_dark,ia_vx) == 8);
  unit_assert(particle.attribute_bytes(it_dark,ia_vy) == 8);
  unit_assert(particle.attribute_bytes(it_dark,ia_vz) == 8);
  unit_assert(particle.attribute_bytes(it_dark,ia_m) == 8);

  unit_assert(particle.attribute_bytes(it_trace,it_x) == 4);
  unit_assert(particle.attribute_bytes(it_trace,it_y) == 8);
  unit_assert(particle.attribute_bytes(it_trace,it_z) == 2);


  unit_func("attribute_offset()");

  // not interleaved
  
  int mp = particle.batch_size();

  unit_assert(particle.attribute_offset(it_dark,ia_x) == 0*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_y) == 4*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_z) == 8*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_vx) == 12*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_vy) == 20*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_vz) == 28*mp);
  unit_assert(particle.attribute_offset(it_dark,ia_m) == 36*mp);

  unit_assert(particle.attribute_bytes(it_trace,it_x) == 4);
  unit_assert(particle.attribute_bytes(it_trace,it_y) == 8);
  unit_assert(particle.attribute_bytes(it_trace,it_z) == 2);

  unit_func("stride");

  // not interleaved
  unit_assert (particle.stride(it_dark,ia_x) == 1);
  unit_assert (particle.stride(it_dark,ia_y) == 1);
  unit_assert (particle.stride(it_dark,ia_z) == 1);
  unit_assert (particle.stride(it_dark,ia_vx) == 1);
  unit_assert (particle.stride(it_dark,ia_vy) == 1);
  unit_assert (particle.stride(it_dark,ia_vz) == 1);

  // interleaved
  unit_assert (particle.stride(it_trace,it_x) == 16/4);
  unit_assert (particle.stride(it_trace,it_y) == 16/8);
  unit_assert (particle.stride(it_trace,it_z) == 16/2);

  //--------------------------------------------------
  //   Insert
  //--------------------------------------------------

  unit_func("num_batches");

  unit_assert(particle.num_batches(it_dark) == 0);
  unit_assert(particle.num_batches(it_trace) == 0);

  unit_func("insert_particles()");

  int i0 = particle.insert_particles (it_dark, 10000);

  unit_assert (i0 == 0);
  unit_func("num_particles()");
  unit_assert (particle.num_particles(it_dark) == 10000);
  unit_func("num_batches()");
  unit_assert (particle.num_batches(it_dark) == 10000/1024 + 1);

  unit_func("index()");
  particle.index(10000,&ib,&ip);
  unit_assert (ib == 10000/1024);
  unit_assert (ip == 10000%1024);
  
  int i1 = particle.insert_particles (it_dark, 10000);

  unit_assert (i1 == 10000);
  unit_func("num_particles()");
  unit_assert (particle.num_particles(it_dark) == 20000);
  unit_func("num_batches()");
  unit_assert (particle.num_batches(it_dark) == 20000/1024 + 1);

  unit_func("index()");
  particle.index(20000,&ib,&ip);
  unit_assert (ib == 20000/1024);
  unit_assert (ip == 20000%1024);

  //--------------------------------------------------
  //   Delete
  //--------------------------------------------------

  int nb = particle.num_batches(it_dark);

  int count_particles = 0;
  int index;
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(it_dark,ib);
    float * x = (float *) particle.attribute_array(it_dark,ib,ia_x);
    float * y = (float *) particle.attribute_array(it_dark,ib,ia_y);
    float * z = (float *) particle.attribute_array(it_dark,ib,ia_z);
    double * vx = (double *) particle.attribute_array(it_dark,ib,ia_vx);
    double * vy = (double *) particle.attribute_array(it_dark,ib,ia_vy);
    double * vz = (double *) particle.attribute_array(it_dark,ib,ia_vz);
    int dx = particle.stride(it_dark,ia_x);
    int dv = particle.stride(it_dark,ia_vx);
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
  unit_assert(count_particles == 20000);

  // run through again and compare values before deleting
  nb = particle.num_batches(it_dark);
  int count_wrong[6];
  for (int i=0; i<6; i++) count_wrong[i] = 0;
  for (int ib=0; ib<nb; ib++) {
    index = ib*mp;
    int np = particle.num_particles(it_dark,ib);
    float * x = (float *) particle.attribute_array(it_dark,ib,ia_x);
    float * y = (float *) particle.attribute_array(it_dark,ib,ia_y);
    float * z = (float *) particle.attribute_array(it_dark,ib,ia_z);
    double * vx = (double *) particle.attribute_array(it_dark,ib,ia_vx);
    double * vy = (double *) particle.attribute_array(it_dark,ib,ia_vy);
    double * vz = (double *) particle.attribute_array(it_dark,ib,ia_vz);
    int dx = particle.stride(it_dark,ia_x);
    int dv = particle.stride(it_dark,ia_vx);
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
  // initialize mask for particles to delete
  bool mask[1024];
  int count_delete = 0;

  for (int i=0; i<1024; i++) {
    mask[i] = (i % 3 == 0);
  }

  unit_func("delete_particles()");
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(it_dark,ib);
    for (int ip=0; ip<np; ip++) if (mask[ip]) count_delete++;
    particle.delete_particles(it_dark,ib,mask);
  }

  unit_assert(particle.num_particles(it_dark) == 20000 - count_delete);

  index = 0;
  nb = particle.num_batches(it_dark);
  for (int i=0; i<6; i++) count_wrong[i] = 0;
  count_particles = 0;
  for (int ib=0; ib<nb; ib++) {
    index = ib*mp;
    int np = particle.num_particles(it_dark,ib);
    float  * x =   (float *) particle.attribute_array(it_dark,ib,ia_x);
    float  * y =   (float *) particle.attribute_array(it_dark,ib,ia_y);
    float  * z =   (float *) particle.attribute_array(it_dark,ib,ia_z);
    double * vx = (double *) particle.attribute_array(it_dark,ib,ia_vx);
    double * vy = (double *) particle.attribute_array(it_dark,ib,ia_vy);
    double * vz = (double *) particle.attribute_array(it_dark,ib,ia_vz);
    int dx = particle.stride(it_dark,ia_x);
    int dv = particle.stride(it_dark,ia_vx);
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

  unit_assert (20000 - count_delete == count_particles);

  unit_func("(re)insert_particles()");

  nb = particle.num_batches(it_dark);
  int ip0 = particle.num_particles(it_dark,nb-1);
  int i2 = particle.insert_particles (it_dark, 10000);
  unit_assert (i2 == 1024*(nb-1) + ip0);

  unit_assert (particle.num_particles() == 10000 + count_particles);

  unit_func("num_particles()");
  nb = particle.num_batches(it_dark);
  int np=particle.num_particles(it_dark);
  int count = 0;
  for (int ib=0; ib<nb; ib++) {
    count += particle.num_particles(it_dark,ib);
  }
  unit_assert (np==count);

  //======================================================================

  unit_func("split_particles()");


  ParticleData * particle_data_2 = new ParticleData;
  Particle particle_2 (particle_descr,particle_data_2);

  count_particles = particle.num_particles(it_dark);

  count_delete = 0;
  for (int i=0; i<1024; i++) {
    mask[i] = (i % 5 == 0);
  }

  unit_func("split_particles()");
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(it_dark,ib);
    for (int ip=0; ip<np; ip++) if (mask[ip]) count_delete++;
    particle.split_particles(it_dark,ib,mask,particle_2);
  }

  unit_assert (count_particles - count_delete
	       == particle.num_particles(it_dark));

  unit_assert (particle_2.num_particles(it_dark) == count_delete);
  


  unit_func("compress()");
  unit_assert (false);
  unit_func("num_particles");
  unit_assert (false);

  //--------------------------------------------------
  //   Grouping
  //--------------------------------------------------


  unit_finalize();

  exit_();


}


PARALLEL_MAIN_END

