// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the FieldBlock class

#include "cello.hpp"

#include "error.hpp"
#include "test.hpp"
#include "field.hpp"

int main()
{

  //----------------------------------------------------------------------
  unit_init();
  //----------------------------------------------------------------------

  FieldDescr field_descr;

  int index_density      = field_descr.insert_field("density");
  int index_velocity_x   = field_descr.insert_field("velocity_x");
  int index_velocity_y   = field_descr.insert_field("velocity_y");
  int index_velocity_z   = field_descr.insert_field("velocity_z");
  int index_total_energy = field_descr.insert_field("total_energy");

  field_descr.set_precision(index_density,     precision_single);
  field_descr.set_precision(index_velocity_x,  precision_double);
  field_descr.set_precision(index_velocity_y,  precision_double);
  field_descr.set_precision(index_velocity_z,  precision_double);
  field_descr.set_precision(index_total_energy,precision_single);

  field_descr.set_ghosts(index_density,      1,1,1);
  field_descr.set_ghosts(index_velocity_x,   2,2,2);
  field_descr.set_ghosts(index_velocity_y,   3,2,1);
  field_descr.set_ghosts(index_velocity_z,   1,2,3);
  field_descr.set_ghosts(index_total_energy, 0,1,2);

  field_descr.set_centering(index_velocity_x, false, true,  true);
  field_descr.set_centering(index_velocity_y, true,  false, true);
  field_descr.set_centering(index_velocity_z, true,  true,  false);

  //  printf ("sizeof(half) = %d\n",sizeof(float16));
  printf ("sizeof(single) = %lu\n",sizeof(float));
  printf ("sizeof(double) = %lu\n",sizeof(double));
  printf ("sizeof(extended) = %lu\n",sizeof(long double));

  //----------------------------------------------------------------------
  unit_class ("FieldBlock");
  //----------------------------------------------------------------------

  FieldBlock field_block;

  //----------------------------------------------------------------------

  unit_func("field_descr");

  field_block.set_field_descr(&field_descr);
  unit_assert (field_block.field_descr() == & field_descr);

  //----------------------------------------------------------------------

  unit_func("dimensions");

  int nx,ny,nz;
  nx=4; ny=5; nz=6;
  field_block.set_dimensions(nx,ny,nz);
  int dimensions[3];
  field_block.dimensions(&dimensions[0],&dimensions[1],&dimensions[2]);
  unit_assert(dimensions[0]==nx && dimensions[1]==ny && dimensions[2]==nz);

  nx=5; ny=3; nz=4;
  field_block.set_dimensions(nx,ny,nz);
  field_block.dimensions(&dimensions[0],&dimensions[1],&dimensions[2]);
  unit_assert(dimensions[0]==nx && dimensions[1]==ny && dimensions[2]==nz);

  //----------------------------------------------------------------------
  // allocate / deallocate
  //----------------------------------------------------------------------

  unit_func("array_allocated");
  unit_assert( ! field_block.array_allocated());

  // Allocate

  unit_func("allocate_array");
  field_block.allocate_array();
  unit_assert(field_block.array() != 0);

  unit_func("array_allocated");
  unit_assert( field_block.array_allocated());

  // Deallocate

  unit_func("deallocate_array");
  field_block.deallocate_array();
  unit_assert(field_block.array() == 0);

  unit_func("array_allocate");
  unit_assert( ! field_block.array_allocated());

  // Allocate

  unit_func("allocate_array");
  field_block.allocate_array();
  unit_assert(field_block.array() != 0);
  
  unit_func("array_allocated");
  unit_assert( field_block.array_allocated());

  
  //----------------------------------------------------------------------

  float *  values_density;
  double * values_velocity_x;
  double * values_velocity_y;
  double * values_velocity_z;
  float * values_total_energy;

  float *  unknowns_density;
  double * unknowns_velocity_x;
  double * unknowns_velocity_y;
  double * unknowns_velocity_z;
  float * unknowns_total_energy;

  int bytes_density;
  int bytes_velocity_x;
  int bytes_velocity_y;
  int bytes_velocity_z;


  unit_func("field_values");  // without ghosts
  
  values_density    = 
    (float *) field_block.field_values(index_density);
  values_velocity_x = 
    (double *) field_block.field_values(index_velocity_x);
  values_velocity_y = 
    (double *) field_block.field_values(index_velocity_y);
  values_velocity_z = 
    (double *) field_block.field_values(index_velocity_z);
  values_total_energy =
    (float *) field_block.field_values(index_total_energy);
  
  unit_assert(values_density != 0);
  unit_assert(values_velocity_x != 0);
  unit_assert(values_velocity_y != 0);
  unit_assert(values_velocity_z != 0);
  unit_assert(values_total_energy != 0);

  bytes_density =    (char *)values_velocity_x - (char *)values_density;
  bytes_velocity_x = (char *)values_velocity_y - (char *)values_velocity_x;
  bytes_velocity_y = (char *)values_velocity_z - (char *)values_velocity_y;
  bytes_velocity_z = (char *)values_total_energy-(char *)values_velocity_z;

  // ghost zone depths
  int gdnx = 1, gdny = 1, gdnz = 1;
  int gvxx = 2, gvxy = 2, gvxz = 2;
  int gvyx = 3, gvyy = 2, gvyz = 1;
  int gvzx = 1, gvzy = 2, gvzz = 3;
  int gtex = 0, gtey = 1, gtez = 2;

  // field sizes without ghosts
  int num_unknowns_density      = (nx)*(ny)*(nz);
  int num_unknowns_velocity_x   = (nx+1)*(ny)*(nz);
  int num_unknowns_velocity_y   = (nx)*(ny+1)*(nz);
  int num_unknowns_velocity_z   = (nx)*(ny)*(nz+1);

  // field sizes with ghosts

  int num_values_density      = (nx+2*gdnx)*(ny+2*gdny)*(nz+2*gdnz);
  int num_values_velocity_x   = (nx+1+2*gvxx)*(ny+2*gvxy)*(nz+2*gvxz);
  int num_values_velocity_y   = (nx+2*gvyx)*(ny+1+2*gvyy)*(nz+2*gvyz);
  int num_values_velocity_z   = (nx+2*gvzx)*(ny+2*gvzy)*(nz+1+2*gvzz);
  
  field_descr.set_ghosts(index_density,      gdnx,gdny,gdnz);
  field_descr.set_ghosts(index_velocity_x,   gvxx,gvxy,gvxz);
  field_descr.set_ghosts(index_velocity_y,   gvyx,gvyy,gvyz);
  field_descr.set_ghosts(index_velocity_z,   gvzx,gvzy,gvzz);
  field_descr.set_ghosts(index_total_energy, gtex,gtey,gtez);

  unit_assert (bytes_density    == 
	       (int)sizeof (float) * num_unknowns_density);
  unit_assert (bytes_velocity_x == 
	       (int)sizeof (double)* num_unknowns_velocity_x);
  unit_assert (bytes_velocity_y == 
	       (int)sizeof (double)* num_unknowns_velocity_y);
  unit_assert (bytes_velocity_z == 
	       (int)sizeof (double)* num_unknowns_velocity_z);

  //----------------------------------------------------------------------

  unit_func("field_unknowns");  // without ghosts

  unknowns_density    = 
    (float *) field_block.field_unknowns(index_density);
  unknowns_velocity_x = 
    (double *) field_block.field_unknowns(index_velocity_x);
  unknowns_velocity_y = 
    (double *) field_block.field_unknowns(index_velocity_y);
  unknowns_velocity_z = 
    (double *) field_block.field_unknowns(index_velocity_z);
  unknowns_total_energy =
    (float *) field_block.field_unknowns(index_total_energy);

  unit_assert(unknowns_density != 0);
  unit_assert(unknowns_velocity_x != 0);
  unit_assert(unknowns_velocity_y != 0);
  unit_assert(unknowns_velocity_z != 0);
  unit_assert(unknowns_total_energy != 0);

  bytes_density =    (char *)unknowns_velocity_x - (char *)unknowns_density;
  bytes_velocity_x = (char *)unknowns_velocity_y - (char *)unknowns_velocity_x;
  bytes_velocity_y = (char *)unknowns_velocity_z - (char *)unknowns_velocity_y;
  bytes_velocity_z = (char *)unknowns_total_energy-(char *)unknowns_velocity_z;

  unit_assert (bytes_density    == 
	       (int)sizeof (float) * num_unknowns_density);
  unit_assert (bytes_velocity_x == 
	       (int)sizeof (double)* num_unknowns_velocity_x);
  unit_assert (bytes_velocity_y == 
	       (int)sizeof (double)* num_unknowns_velocity_y);
  unit_assert (bytes_velocity_z == 
	       (int)sizeof (double)* num_unknowns_velocity_z);


  //----------------------------------------------------------------------

  unit_func("allocate_ghosts");  // with ghosts

  field_block.allocate_ghosts();
  
  values_density    = 
    (float *) field_block.field_values(index_density);

  values_velocity_x = 
    (double *) field_block.field_values(index_velocity_x);
  values_velocity_y = 
    (double *) field_block.field_values(index_velocity_y);
  values_velocity_z = 
    (double *) field_block.field_values(index_velocity_z);
  values_total_energy =
    (float *) field_block.field_values(index_total_energy);
  
  unit_assert(values_density != 0);
  unit_assert(values_velocity_x != 0);
  unit_assert(values_velocity_y != 0);
  unit_assert(values_velocity_z != 0);
  unit_assert(values_total_energy != 0);

  bytes_density =    (char *)values_velocity_x - (char *)values_density;
  bytes_velocity_x = (char *)values_velocity_y - (char *)values_velocity_x;
  bytes_velocity_y = (char *)values_velocity_z - (char *)values_velocity_y;
  bytes_velocity_z = (char *)values_total_energy-(char *)values_velocity_z;

  unit_assert (bytes_density    == 
	       (int)sizeof (float) * num_values_density);
  unit_assert (bytes_velocity_x == 
	       (int)sizeof (double)* num_values_velocity_x);
  unit_assert (bytes_velocity_y == 
	       (int)sizeof (double)* num_values_velocity_y);
  unit_assert (bytes_velocity_z == 
	       (int)sizeof (double)* num_values_velocity_z);

  unit_func("field_unknowns");  // with ghosts

  unknowns_density    = 
    (float *) field_block.field_unknowns(index_density);

  unknowns_velocity_x = 
    (double *) field_block.field_unknowns(index_velocity_x);
  unknowns_velocity_y = 
    (double *) field_block.field_unknowns(index_velocity_y);
  unknowns_velocity_z = 
    (double *) field_block.field_unknowns(index_velocity_z);
  unknowns_total_energy =
    (float *) field_block.field_unknowns(index_total_energy);

  unit_assert(unknowns_density != 0);
  unit_assert(unknowns_velocity_x != 0);
  unit_assert(unknowns_velocity_y != 0);
  unit_assert(unknowns_velocity_z != 0);
  unit_assert(unknowns_total_energy != 0);

  bytes_density =    (char *)unknowns_velocity_x - (char *)unknowns_density;
  bytes_velocity_x = (char *)unknowns_velocity_y - (char *)unknowns_velocity_x;
  bytes_velocity_y = (char *)unknowns_velocity_z - (char *)unknowns_velocity_y;
  bytes_velocity_z = (char *)unknowns_total_energy-(char *)unknowns_velocity_z;

  // a,b fields  u unknowns  v values  g ghosts
  // bu - au = (bv + bg) - (av + ag) 
  //         = (bv + (bu-bv)) - (av + (au-av))
  //         = (bv - av) + (bu-bv) - (au-av)


  
  int ghost_offset_density    = 
    (int)sizeof (float) * (gdnx + (nx+2*gdnx) *(gdny + (ny+2*gdny)*gdnz));
  int ghost_offset_velocity_x = 
    (int)sizeof (double) * (gvxx + (nx+1+2*gvxx)*(gvxy + (ny+2*gvxy)*gvxz));
  int ghost_offset_velocity_y = 
    (int)sizeof (double) * (gvyx + (nx+2*gvyx)*(gvyy + (ny+1+2*gvyy)*gvyz));
  int ghost_offset_velocity_z = 
    (int)sizeof (double) * (gvzx + (nx+2*gvzx)*(gvzy + (ny+2*gvzy)*gvzz));
  int ghost_offset_total_energy = 
    (int)sizeof (float) * (gtex + (nx+2*gtex)*(gtey + (ny+2*gtey)*gtez));

  unit_assert (bytes_density    == 
	       (int)sizeof (float) * num_values_density + 
	       ghost_offset_velocity_x - ghost_offset_density);
  unit_assert (bytes_velocity_x == 
	       (int)sizeof (double)* num_values_velocity_x +
	       ghost_offset_velocity_y - ghost_offset_velocity_x);
  unit_assert (bytes_velocity_y == 
	       (int)sizeof (double)* num_values_velocity_y +
	       ghost_offset_velocity_z - ghost_offset_velocity_y);
  unit_assert (bytes_velocity_z == 
	       (int)sizeof (double)* num_values_velocity_z +
	       ghost_offset_total_energy - ghost_offset_velocity_z);

  //----------------------------------------------------------------------
  unit_func("box_extent");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("cell_width");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("clear");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("ghosts_allocated");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("allocate_ghosts");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("deallocate_ghosts");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("split");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("merge");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("read");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("write");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------
}
