// $Id: user_MethodEnzoControl.cpp 1388 2010-04-20 23:57:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     user_MethodEnzoControl.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of MethodEnzoControl user-dependent class member functions

#include "method.hpp"
#include "user.hpp"
#include "error.hpp"
#include "cello_hydro.h"

//----------------------------------------------------------------------

void MethodEnzoControl::initialize (DataDescr * data_descr) throw()
{
}

//----------------------------------------------------------------------

void MethodEnzoControl::finalize (DataDescr * data_descr) throw()
{
}

//----------------------------------------------------------------------

void MethodEnzoControl::initialize_block ( DataBlock * data_block ) throw ()
{

  FieldBlock * field_block = data_block->field_block();
  
  double xm,xp,ym,yp,zm,zp;

  field_block->box_extent(&xm,&xp,&ym,&yp,&zm,&zp);

  GridLeftEdge[0]    = xm;
  GridLeftEdge[1]    = ym;
  GridLeftEdge[2]    = zm;

  // Get block dimensions

  int nx,ny,nz;
  field_block -> dimensions (&nx,&ny,&nz);

  CycleNumber = 0;
  dt          = 1;

  // Grid dimensions

  GridDimension[0]  = nx + 2*ghost_depth[0];
  GridDimension[1]  = ny + 2*ghost_depth[1];
  GridDimension[2]  = nz + 2*ghost_depth[2];
  GridStartIndex[0] = ghost_depth[0];
  GridStartIndex[1] = ghost_depth[1];
  GridStartIndex[2] = ghost_depth[2];
  GridEndIndex[0]   = nx + ghost_depth[0] - 1;
  GridEndIndex[1]   = ny + ghost_depth[1] - 1;
  GridEndIndex[2]   = nz + ghost_depth[2] - 1;

  // Initialize CellWidth[].  Should be converted to constants hx,hy,hz

  double h3[3];
  field_block->cell_width(&h3[0],&h3[1],&h3[2]);

  for (int dim=0; dim<GridRank; dim++) {
    CellWidth[dim] = new ENZO_FLOAT[GridDimension[dim]];
    for (int i=0; i<GridDimension[dim]; i++) {
      CellWidth[dim][i] = h3[dim];
    }
  }

  // Initialize BaryonField[] pointers

  for (int field = 0; field < NumberOfBaryonFields; field++) {
    BaryonField[field] = (float *)field_block->field_values(field);
    printf ("field %d %g\n",field,BaryonField[field][0]);
  }
}


//----------------------------------------------------------------------

void MethodEnzoControl::finalize_block ( DataBlock * data_block ) throw ()
{
  for (int dim=0; dim<GridRank; dim++) {
    delete [] CellWidth[dim];
  }
}

//----------------------------------------------------------------------

void MethodEnzoControl::refresh_block(DataBlock * data_block) throw()
{
  INCOMPLETE_MESSAGE ("MethodEnzoControl::refresh_block","");
}

//======================================================================

