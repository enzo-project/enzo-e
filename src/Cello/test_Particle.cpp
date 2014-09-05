// See LICENSE_CELLO file for license and copyright information

/// @file     test_Particle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 15:36:02 PDT 2014
/// @brief    Test program for the Particle class

#include "main.hpp"
#include "test.hpp"

#include "particle.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ParticleDescr");
  ParticleDescr pd (3);

  // unit_class("ParticleBlock");
  // ParticleBlock pb = new ParticleBlock;
  // unit_assert (pb != NULL);

  //--------------------------------------------------

  unit_func ("count()");

  unit_func ("set_attribute_size()");
  pd.set_attribute_size(particle_attribute_id, 8);
  pd.set_attribute_size(particle_attribute_position,6);
  pd.set_attribute_size(particle_attribute_velocity,6);
  pd.set_attribute_size(particle_attribute_mass,8);

  unit_func ("attribute_size()");
  unit_assert(pd.attribute_size(particle_attribute_id) == 8);
  unit_assert(pd.attribute_size(particle_attribute_position) == 6);
  unit_assert(pd.attribute_size(particle_attribute_velocity) == 6);
  unit_assert(pd.attribute_size(particle_attribute_mass) == 8);

  unit_assert (pd.num_types() == 0);
  const int it = pd.new_type("Tracer");
  unit_assert (pd.num_types() == 1);
  const int id = pd.new_type("DarkMatter");
  unit_assert (pd.num_types() == 2);
  const int is = pd.new_type("Sink");
  unit_assert (pd.num_types() == 3);

  unit_func ("type()");
  unit_assert(pd.type(it)=="Tracer");
  unit_assert(pd.type(id)=="DarkMatter");
  unit_assert(pd.type(is)=="Sink");

  unit_func ("add_attribute()");
  pd.add_attribute(it,particle_attribute_id);
  pd.add_attribute(it,particle_attribute_position);
  unit_func ("num_attributes()");
  unit_assert(pd.num_attributes(it)==2);

  pd.add_attribute(id,particle_attribute_id);
  pd.add_attribute(id,particle_attribute_position);
  pd.add_attribute(id,particle_attribute_velocity);
  pd.add_attribute(id,particle_attribute_mass);
  unit_func ("num_attributes()");
  unit_assert(pd.num_attributes(id)==4);

  pd.add_attribute(is,particle_attribute_position);
  pd.add_attribute(is,particle_attribute_velocity);
  pd.add_attribute(is,particle_attribute_mass);
  unit_func ("num_attributes()");
  unit_assert(pd.num_attributes(is)==3);

  unit_func ("particle_size()");
  unit_assert(pd.particle_size(it)== 1 + 8 + 6);
  unit_assert(pd.particle_size(id)== 1 + 8 + 6 + 6 + 8);
  unit_assert(pd.particle_size(is)== 1 +     6 + 6 + 8);

  unit_func ("attribute_offset()");
  unit_assert(pd.attribute_offset(it,0)==1);
  unit_assert(pd.attribute_offset(it,1)==1+8);

  unit_assert(pd.attribute_offset(id,0)==1);
  unit_assert(pd.attribute_offset(id,1)==1+8);
  unit_assert(pd.attribute_offset(id,2)==1+8+6);
  unit_assert(pd.attribute_offset(id,3)==1+8+6+6);

  unit_assert(pd.attribute_offset(is,0)==1);
  unit_assert(pd.attribute_offset(is,1)==1+6);
  unit_assert(pd.attribute_offset(is,2)==1+6+6);

  //--------------------------------------------------

  const double lower[3] = {-0.5,0.0,-1.0};
  const double upper[3] = { 0.5,2.0, 0.0};
  ParticleBlock pb(&pd);

  // Create actual particles of each type

  pb.create(&pd,it,1000);
  pb.create(&pd,id,500);
  pb.create(&pd,is,20);
 

  std::vector <char> data;
  pb.data(it,&data);
  
  printf ("data size = %ld\n",data.size());

  unit_func("ParticleBlock::data()");
  unit_assert(data.size() == pb.particle_count(it)*pd.particle_size(it));

  // local positions

  unsigned ** u;
  u = new unsigned * [3];
  u[0] = new unsigned[1000];
  u[1] = new unsigned[1000];
  u[2] = new unsigned[1000];

  // global positions

  long double ** a;
  a = new long double * [3];
  a[0] = new long double[1000];
  a[1] = new long double[1000];
  a[2] = new long double[1000];

  long unsigned k = RAND_MAX;
  int count = 0;
  while (k) {
    k = k >> 1;
    count ++;
  }
  
  printf ("RAND_MAX = %d bits = %d\n",RAND_MAX,count);
  for (int i=0; i<1000; i++) {
    u[0][i] = rand();
    u[1][i] = rand();
    u[2][i] = rand();

    a[0][i] = lower[0] + ((upper[0]-lower[0])*u[0][i])/RAND_MAX;
    a[1][i] = lower[1] + ((upper[1]-lower[1])*u[1][i])/RAND_MAX;
    a[2][i] = lower[2] + ((upper[2]-lower[2])*u[2][i])/RAND_MAX;

  }

  unit_func ("set_global_positions()");
  pb.set_global_positions (&pd,it,(const double **)a);

  unit_func ("local_positions()");
  pb.local_positions (&pd,it,(unsigned **)u);

  //--------------------------------------------------

  unit_func ("set_local_positions()");
  pb.set_local_positions (&pd,it,(const unsigned **)u);

  unit_func ("global_positions()");
  pb.global_positions (&pd,it,(double **)a);


  delete [] u[0];
  delete [] u[1];
  delete [] u[2];
  delete [] u;

  delete [] a[0];
  delete [] a[1];
  delete [] a[2];
  delete [] a;

  unit_finalize();

  exit_();


}


PARALLEL_MAIN_END

