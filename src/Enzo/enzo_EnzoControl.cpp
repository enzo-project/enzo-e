// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoControl.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoControl user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoControl::initialize (DataDescr * data_descr) throw()
{



}

//----------------------------------------------------------------------

void EnzoControl::finalize (DataDescr * data_descr) throw()
{
}

//----------------------------------------------------------------------

void EnzoControl::initialize_block ( DataBlock * data_block ) throw ()
{

  FieldBlock * field_block = data_block->field_block();
  
  double xm,xp,ym,yp,zm,zp;

  field_block->extent(&xm,&xp,&ym,&yp,&zm,&zp);

  enzo_->GridLeftEdge[0]    = xm;
  enzo_->GridLeftEdge[1]    = ym;
  enzo_->GridLeftEdge[2]    = zm;

  // Grid dimensions

  int nx,ny,nz;
  field_block -> size (&nx,&ny,&nz);

  enzo_->GridDimension[0]  = nx + 2*enzo_->ghost_depth[0];
  enzo_->GridDimension[1]  = ny + 2*enzo_->ghost_depth[1];
  enzo_->GridDimension[2]  = nz + 2*enzo_->ghost_depth[2];
  enzo_->GridStartIndex[0] = enzo_->ghost_depth[0];
  enzo_->GridStartIndex[1] = enzo_->ghost_depth[1];
  enzo_->GridStartIndex[2] = enzo_->ghost_depth[2];
  enzo_->GridEndIndex[0]   = nx + enzo_->ghost_depth[0] - 1;
  enzo_->GridEndIndex[1]   = ny + enzo_->ghost_depth[1] - 1;
  enzo_->GridEndIndex[2]   = nz + enzo_->ghost_depth[2] - 1;

  // Initialize CellWidth[].  Should be converted to constants hx,hy,hz

  double h3[3];
  field_block->cell_width(&h3[0],&h3[1],&h3[2]);

  for (int dim=0; dim<enzo_->GridRank; dim++) {
    enzo_->CellWidth[dim] = new ENZO_FLOAT[enzo_->GridDimension[dim]];
    for (int i=0; i<enzo_->GridDimension[dim]; i++) {
      enzo_->CellWidth[dim][i] = h3[dim];
    }
  }

  // Initialize BaryonField[] pointers

  for (int field = 0; field < enzo_->NumberOfBaryonFields; field++) {
    enzo_->BaryonField[field] = (float *)field_block->field_values(field);
  }

 
  //   // Boundary
  //   /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
  //      set it to one. */
 
  //   /* 1) Compute Courant condition for baryons. */
 
  //    // boundary
 
  //    BoundaryRank = 2;
  //    BoundaryDimension[0] = GridDimension[0];
  //    BoundaryDimension[1] = GridDimension[1];

  //    for (int field=0; field<NumberOfBaryonFields; field++) {
  //      BoundaryFieldType[field] = enzo_->FieldType[field];
  //      for (int dim = 0; dim < 3; dim++) {
  //        for (int face = 0; face < 2; face++) {
  //  	int n1 = GridDimension[(dim+1)%3];
  //  	int n2 = GridDimension[(dim+2)%3];
  //  	int size = n1*n2;
  //  	BoundaryType [field][dim][face] = new bc_type [size];
  //  	BoundaryValue[field][dim][face] = NULL;
  //  	for (int i2 = 0; i2<n2; i2++) {
  //  	  for (int i1 = 0; i1<n1; i1++) {
  //  	    int i = i1 + n1*i2;
  //  	    BoundaryType[field][dim][face][i] = bc_reflecting;
  //  	  }
  //  	}
  //        }
  //      }
  //    }

}


//----------------------------------------------------------------------

void EnzoControl::finalize_block ( DataBlock * data_block ) throw ()
{
  ++ enzo_->CycleNumber;
  for (int dim=0; dim < enzo_->GridRank; dim++) {
    delete [] enzo_->CellWidth[dim];
  }
}

//----------------------------------------------------------------------

void EnzoControl::refresh_ghost(DataBlock * data_block, 
				bool xm, bool xp, 
				bool ym, bool yp, 
				bool zm, bool zp) throw()
{

}

//======================================================================

