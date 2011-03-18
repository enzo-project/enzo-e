// $Id: test_EnzoMethodPpm.hpp 1942 2011-01-20 00:53:45Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef TEST_ENZO_METHOD_PPM_HPP
#define TEST_ENZO_METHOD_PPM_HPP

/// @file     test_EnzoMethodPpm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Feb  4 13:13:50 PST 2011
/// @todo     Salvage code for e.g. method_Output.cpp or enzo_OutputEnzo.cpp
/// @brief    [\ref Enzo] Convenience functions for [defunct] test_EnzoMethodPpm.cpp


//----------------------------------------------------------------------

void output_images(int cycle,
		   FieldBlock * field_block,
		   int cycle_image = 1)
{

  if (! (cycle_image && cycle % cycle_image == 0)) return;

  const FieldDescr * field_descr = field_block->field_descr();
  Monitor * monitor = Monitor::instance();
  int nx,ny,nz;
  int gx,gy,gz;
  int mx,my,mz;
  field_block->enforce_boundary(boundary_reflecting);
  field_block->size(&nx,&ny,&nz);
  int count = field_descr->field_count();
  for (int index = 0; index < count; index++) {
    field_descr->ghosts(index,&gx,&gy,&gz);
    mx=nx+2*gx;
    my=ny+2*gy;
    mz=nz+2*gz;
    char filename[80];
    std::string field_name = field_descr->field_name(index);
    Scalar * field_values = (Scalar *)field_block->field_values(index);
    sprintf (filename,"ppm-%s-%05d.png",field_name.c_str(),cycle);
    monitor->image (filename, field_values, mx,my,mz, mx,my,mz, 0,0,0, 2, reduce_sum, 0.0, 1.0);
  }
}

//----------------------------------------------------------------------

void output_dump(FileHdf5 & hdf5,
		 int cycle,
		 FieldBlock * field_block,
		 int cycle_dump = 1)
{

  // Exit if we don't dump data this cycle
  if (! (cycle_dump && cycle % cycle_dump == 0)) return;

  // Refresh boundary conditions 
  // (should have check to not do it more than once)

  const FieldDescr * field_descr = field_block->field_descr();
  field_block->enforce_boundary(boundary_reflecting);

  // Open file

  char filename[80];
  sprintf (filename,"ppm-%05d.h5",cycle);
  hdf5.open_file (filename,"w");

  // Get block size
  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  // Loop over fields in block

  int count = field_descr->field_count();

  for (int index = 0; index < count; index++) {

    // Get field's ghost zone depth
    int gx,gy,gz;
    field_descr->ghosts(index,&gx,&gy,&gz);

    // Get field's total size including ghosts

    int mx,my,mz;
    mx=nx+2*gx;
    my=ny+2*gy;
    mz=nz+2*gz;

    // Get field's name for filename

    std::string field_name = field_descr->field_name(index);

    // Write block, including ghosts
    // (should include option to omit ghosts)
    // (what if ghosts aren't allocated in general?)

    // prepare to write the field to the file

    hdf5.open_dataset (field_name,
		       field_descr->precision(index),
		       mx,my,mz);

    // write the field to the file

    hdf5.write((char *)field_block->field_values(index),
		field_descr->precision(index));

    // close the field in the file
    hdf5.close_dataset ();
  }

  // close the file

  hdf5.close_file();

}

//----------------------------------------------------------------------

void output_progress (Monitor * monitor,
		      int cycle,
		      double time,
		      double dt,
		      int cycle_progress = 1
		      )
{
  if (! (cycle_progress && cycle % cycle_progress == 0)) return;
  char buffer[100];
  sprintf (buffer," cycle = %05d  sim-time = %10.8f  dt = %10.8f",
	   cycle,time,dt);
  monitor->print (buffer);
}


#endif /* TEST_ENZO_METHOD_PPM_HPP */
