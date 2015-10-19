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
  ParticleDescr particle_descr;

  unit_class("ParticleData");
  ParticleData particle_data;

  Particle particle (&particle_descr, &particle_data);

  //--------------------------------------------------
  //   Type
  //--------------------------------------------------

  // new_type()
  // num_types()

  unit_func ("num_types()");
  unit_assert (particle.num_types() == 0);
  const int i_dark  = particle.new_type ("dark");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 1);
  unit_func ("num_types()");
  unit_assert (i_dark == 0);
  const int i_trace = particle.new_type ("trace");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 2);
  unit_func ("num_types()");
  unit_assert (i_trace == 1);
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
  unit_assert(particle.type_name(i_dark) == "dark");
  unit_assert(particle.type_name(i_trace) == "trace");
  unit_assert(particle.type_name(i_star) == "star");
  unit_assert(particle.type_name(i_sink) == "sink");

  unit_func ("type_index()");
  unit_assert(particle.type_index("dark")  == i_dark);
  unit_assert(particle.type_index("trace") == i_trace);
  unit_assert(particle.type_index("star")  == i_star);
  unit_assert(particle.type_index("sink")  == i_sink);

  
  //--------------------------------------------------
  // Interleaved
  //--------------------------------------------------

  unit_func("interleaved");
  unit_assert(particle.interleaved(i_dark) == false);
  unit_assert(particle.interleaved(i_trace) == false);

  unit_func("set_interleaved");
  particle.set_interleaved(i_dark,true);
  unit_assert(particle.interleaved(i_dark) == true);
  unit_assert(particle.interleaved(i_trace) == false);
  particle.set_interleaved(i_trace,true);
  unit_assert(particle.interleaved(i_dark) == true);
  unit_assert(particle.interleaved(i_trace) == true);
  particle.set_interleaved(i_dark,false);
  unit_assert(particle.interleaved(i_dark) == false);
  unit_assert(particle.interleaved(i_trace) == true);

  //--------------------------------------------------
  //  Attribute
  //--------------------------------------------------

  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 0);
  unit_func ("new_attribute");
  const int id_x = particle.new_attribute (i_dark, "position_x", 4);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 1);
  unit_func ("new_attribute");
  const int id_y = particle.new_attribute (i_dark, "position_y", 4);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 2);
  unit_func ("new_attribute");
  const int id_z = particle.new_attribute (i_dark, "position_z", 4);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 3);
  unit_func ("new_attribute");
  const int id_vx = particle.new_attribute (i_dark, "velocity_x", 8);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 4);
  unit_func ("new_attribute");
  const int id_vy = particle.new_attribute (i_dark, "velocity_y", 8);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 5);
  unit_func ("new_attribute");
  const int id_vz = particle.new_attribute (i_dark, "velocity_z", 8);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 6);
  unit_func ("new_attribute");
  const int id_m = particle.new_attribute (i_dark, "mass", 8);
  
  unit_func ("attribute_name()");
  unit_assert(particle.attribute_name(i_dark,id_x) == "position_x");
  unit_assert(particle.attribute_name(i_dark,id_y) == "position_y");
  unit_assert(particle.attribute_name(i_dark,id_z) == "position_z");
  unit_assert(particle.attribute_name(i_dark,id_vx) == "velocity_x");
  unit_assert(particle.attribute_name(i_dark,id_vy) == "velocity_y");
  unit_assert(particle.attribute_name(i_dark,id_vz) == "velocity_z");

  unit_func ("attribute_index()");
  unit_assert(particle.attribute_index(i_dark,"position_x")  == id_x);

  const int it_x = particle.new_attribute (i_trace, "position_x", 4);
  const int it_y = particle.new_attribute (i_trace, "position_y", 8);
  const int it_z = particle.new_attribute (i_trace, "position_z", 2);

  //--------------------------------------------------
  //   Bytes
  //--------------------------------------------------

  unit_func("attribute_bytes(it,ia)");
  unit_assert(particle.attribute_bytes(i_dark,id_x) == 4);
  unit_assert(particle.attribute_bytes(i_dark,id_y) == 4);
  unit_assert(particle.attribute_bytes(i_dark,id_z) == 4);
  unit_assert(particle.attribute_bytes(i_dark,id_vx) == 8);
  unit_assert(particle.attribute_bytes(i_dark,id_vy) == 8);
  unit_assert(particle.attribute_bytes(i_dark,id_vz) == 8);
  unit_assert(particle.attribute_bytes(i_dark,id_m) == 8);
  unit_func("attribute_bytes(it)");
  // needs to be divisible by 8 for stride computation:
  // 4+4+4+8+8+8+8 rounded up to closest multiple of 8
  unit_assert(particle.attribute_bytes(i_dark) == 48);

  unit_func("attribute_offset(it,ia)");
  unit_assert(particle.attribute_offset(i_dark,id_x) == 0);
  unit_assert(particle.attribute_offset(i_dark,id_y) == 4);
  unit_assert(particle.attribute_offset(i_dark,id_z) == 8);
  unit_assert(particle.attribute_offset(i_dark,id_vx) == 16);
  unit_assert(particle.attribute_offset(i_dark,id_vy) == 24);
  unit_assert(particle.attribute_offset(i_dark,id_vz) == 32);
  unit_assert(particle.attribute_offset(i_dark,id_m) == 40);

  unit_func("attribute_bytes(it,ia)");
  unit_assert(particle.attribute_bytes(i_trace,it_x) == 4);
  unit_assert(particle.attribute_bytes(i_trace,it_y) == 8);
  unit_assert(particle.attribute_bytes(i_trace,it_z) == 2);
  unit_func("attribute_bytes(it)");
  unit_assert(particle.attribute_bytes(i_trace) == 16);

  unit_func("stride");
  unit_assert (particle.stride(i_dark,id_x) == 1);
  unit_assert (particle.stride(i_dark,id_y) == 1);
  unit_assert (particle.stride(i_dark,id_z) == 1);
  unit_assert (particle.stride(i_dark,id_vx) == 1);
  unit_assert (particle.stride(i_dark,id_vy) == 1);
  unit_assert (particle.stride(i_dark,id_vz) == 1);
  unit_assert (particle.stride(i_trace,it_x) == 16/4);
  unit_assert (particle.stride(i_trace,it_y) == 16/8);
  unit_assert (particle.stride(i_trace,it_z) == 16/2);

  //--------------------------------------------------
  //   Batch
  //--------------------------------------------------

  unit_func("batch_size()");
  unit_assert (particle.batch_size() == 1);
  unit_func("set_batch_size()");
  particle.set_batch_size(1024);
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
  //   Insert
  //--------------------------------------------------

  unit_func("num_batches");

  unit_assert(particle.num_batches(i_dark) == 0);
  unit_assert(particle.num_batches(i_trace) == 0);

  unit_func("insert_particles()");

  int i0 = particle.insert_particles (i_dark, 10000);

  unit_assert (i0 == 0);
  unit_func("num_particles()");
  unit_assert (particle.num_particles(i_dark) == 10000);
  unit_func("num_batches()");
  unit_assert (particle.num_batches(i_dark) == 10000/1024 + 1);

  unit_func("index()");
  particle.index(10000,&ib,&ip);
  unit_assert (ib == 10000/1024);
  unit_assert (ip == 10000%1024);
  
  int i1 = particle.insert_particles (i_dark, 10000);

  unit_assert (i1 == 10000);
  unit_func("num_particles()");
  unit_assert (particle.num_particles(i_dark) == 20000);
  unit_func("num_batches()");
  unit_assert (particle.num_batches(i_dark) == 20000/1024 + 1);

  unit_func("index()");
  particle.index(20000,&ib,&ip);
  unit_assert (ib == 20000/1024);
  unit_assert (ip == 20000%1024);

  //--------------------------------------------------
  //   Delete
  //--------------------------------------------------

  int nb = particle.num_batches(i_dark);

  int count = 0;
  int index=0;
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(i_dark,ib);
    float * x = (float *) particle.attribute_array(i_dark,ib,id_x);
    float * y = (float *) particle.attribute_array(i_dark,ib,id_y);
    float * z = (float *) particle.attribute_array(i_dark,ib,id_z);
    double * vx = (double *) particle.attribute_array(i_dark,ib,id_vx);
    double * vy = (double *) particle.attribute_array(i_dark,ib,id_vy);
    double * vz = (double *) particle.attribute_array(i_dark,ib,id_vz);
    int dx = particle.stride(i_dark,id_x);
    int dv = particle.stride(i_dark,id_vx);
    for (int ip=0,ix=0,iv=0; ip<np; ip++,ix+=dx,iv+=dv) {
      count ++;
      x[ix]  = 10*index;
      y[ix]  = 10*index+1;
      z[ix]  = 10*index+2;
      vx[iv] = 10*index+3;
      vy[iv] = 10*index+4;
      vz[iv] = 10*index+5;
      index++;
    }
  }
  unit_assert(count == 20000);

  // run through again and compare values before deleting
  nb = particle.num_batches(i_dark);
  index = 0;
  printf ("nb %d\n",nb);
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(i_dark,ib);
    printf ("ib np %d %d\n",ib,np);
    float * x = (float *) particle.attribute_array(i_dark,ib,id_x);
    float * y = (float *) particle.attribute_array(i_dark,ib,id_y);
    float * z = (float *) particle.attribute_array(i_dark,ib,id_z);
    double * vx = (double *) particle.attribute_array(i_dark,ib,id_vx);
    double * vy = (double *) particle.attribute_array(i_dark,ib,id_vy);
    double * vz = (double *) particle.attribute_array(i_dark,ib,id_vz);
    int dx = particle.stride(i_dark,id_x);
    int dv = particle.stride(i_dark,id_vx);
    for (int ip=0,ix=0,iv=0; ip<np; ip++,ix+=dx,iv+=dv) {
      if (x[ix] != 10*index)    printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,x[ix],1.0*10*index);
      if (y[ix] != 10*index+1)  printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,y[ix],1.0*10*index+1);
      if (z[ix] != 10*index+2)  printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,z[ix],1.0*10*index+2);
      if (vx[iv] != 10*index+3) printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,vx[ix],1.0*10*index+3);
      if (vy[iv] != 10*index+4) printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,vy[ix],1.0*10*index+4);
      if (vz[iv] != 10*index+5) printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,vz[ix],1.0*10*index+5);
      index++;
    }
  }

  // initialize mask for particles to delete
  bool mask[1024];
  count = 0;
  for (int i=0; i<1024; i++) {
    mask[i] = (i % 3 == 0);
    if (mask[i]) count++;
  }

  printf ("deleting %d particles\n",count);
  unit_func("delete_particles()");
  particle.delete_particles(i_dark,0,mask);

  printf ("num particles = %d\n",particle.num_particles(i_dark));

  unit_assert(particle.num_particles(i_dark) == 20000 - count);

  int count2 = 0;
  index = 0;
  nb = particle.num_batches(i_dark);
  printf ("nb %d\n",nb);
  for (int ib=0; ib<nb; ib++) {
    int np = particle.num_particles(i_dark,ib);
    printf ("ib np %d %d\n",ib,np);
    float * x = (float *) particle.attribute_array(i_dark,ib,id_x);
    float * y = (float *) particle.attribute_array(i_dark,ib,id_y);
    float * z = (float *) particle.attribute_array(i_dark,ib,id_z);
    double * vx = (double *) particle.attribute_array(i_dark,ib,id_vx);
    double * vy = (double *) particle.attribute_array(i_dark,ib,id_vy);
    double * vz = (double *) particle.attribute_array(i_dark,ib,id_vz);
    int dx = particle.stride(i_dark,id_x);
    int dv = particle.stride(i_dark,id_vx);
    for (int ip=0,ix=0,iv=0; ip<np; ip++,ix+=dx,iv+=dv) {
      count2 ++;
      if (index % 3) index++;
      if (x[ix] != 10*index)    printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,x[ix],1.0*10*index);
      if (y[ix] != 10*index+1)  printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,y[ix],1.0*10*index+1);
      if (z[ix] != 10*index+2)  printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,z[ix],1.0*10*index+2);
      if (vx[iv] != 10*index+3) printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,vx[ix],1.0*10*index+3);
      if (vy[iv] != 10*index+4) printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,vy[ix],1.0*10*index+4);
      if (vz[iv] != 10*index+5) printf ("%s:%d ib %d index %d value %g expected %g\n",__FILE__,__LINE__,ib,index,vz[ix],1.0*10*index+5);
      index++;
    }
  }

  unit_assert (count == count2);
  unit_func("split_particles()");
  unit_assert (false);
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

