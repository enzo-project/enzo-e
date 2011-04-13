// $Id: enzo_EnzoBoundary.cpp 2085 2011-03-10 01:16:42Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBoundary.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @todo     Merge similar functions enforce_[reflecting|outflow...]_()
/// @brief    Implements the EnzoBoundary class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoBoundary::EnzoBoundary (boundary_type_enum boundary_type) throw()
  : Boundary(),
    boundary_type_ (boundary_type)
{ 
}

//----------------------------------------------------------------------

void EnzoBoundary::enforce 
(
 const FieldDescr * field_descr,
 Block * block,
 face_enum     face,
 axis_enum     axis 
 ) const throw()
{
  if (face == face_all) {
    // WARNING: recursive
    enforce(field_descr,block,face_lower,axis);
    // WARNING: recursive
    enforce(field_descr,block,face_upper,axis);
  } else if (axis == axis_all) {
    // WARNING: recursive
    enforce(field_descr,block,face,axis_x);
    // WARNING: recursive
    enforce(field_descr,block,face,axis_y);
    // WARNING: recursive
    enforce(field_descr,block,face,axis_z);
  } else {
    FieldBlock * field_block = block->field_block();
    if ( ! field_block->ghosts_allocated() ) {
      ERROR("EnzoBoundary::enforce",
	    "Function called with ghosts not allocated");
    }
    switch (boundary_type_) {
    case boundary_type_reflecting:
      enforce_reflecting_(field_descr, field_block,face,axis);
      break;
    case boundary_type_outflow:
      enforce_outflow_   (field_descr, field_block,face,axis);
      break;
    case boundary_type_inflow:
      enforce_inflow_    (field_descr, field_block,face,axis);
      break;
    case boundary_type_periodic:
      // Periodic handled by ghost refresh
      break;
    default:
      ERROR("EnzoBoundary::enforce",
	    "Undefined boundary type");
      break;
    }
  }
}

//----------------------------------------------------------------------

void EnzoBoundary::enforce_reflecting_
(
 const FieldDescr * field_descr,
 FieldBlock       * field_block,
 face_enum face,
 axis_enum axis
 ) const throw()
{

  int nx,ny,nz;
  int gx,gy,gz;
  field_block->size(&nx,&ny,&nz);
  for (int field = 0; field < field_descr->field_count(); field++) {
    field_descr->ghosts(field,&gx,&gy,&gz);
    void * array = field_block->field_values(field);
    bool vx =      (field_descr->field_name(field) == "velocity_x");
    bool vy =      (field_descr->field_name(field) == "velocity_y");
    bool vz =      (field_descr->field_name(field) == "velocity_z");
    precision_enum precision = field_descr->precision(field);
    switch (precision) {
    case precision_single:
      enforce_reflecting_precision_(face,axis, (float *)array,
				    nx,ny,nz, gx,gy,gz,  vx,vy,vz);
      break;
    case precision_double:
      enforce_reflecting_precision_(face,axis, (double *)array,
				    nx,ny,nz, gx,gy,gz,  vx,vy,vz);
      break;
    case precision_extended80:
    case precision_extended96:
    case precision_quadruple:
      enforce_reflecting_precision_(face,axis, (long double *)array,
				    nx,ny,nz, gx,gy,gz, vx,vy,vz);
      break;
    default:
      break;
    }
  }
}

//----------------------------------------------------------------------

template<class T>
void EnzoBoundary::enforce_reflecting_precision_
(
 face_enum face, 
 axis_enum axis,
 T * array,
 int nx,int ny,int nz,
 int gx,int gy,int gz,
 bool vx,bool vy,bool vz
 ) const throw()
{
  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  int ix,iy,iz,ig;
  T sign;

  if (face == face_lower && 
      axis == axis_x) {
    if (nx > 1) {
      ix = gx;
      sign = vx ? -1.0 : 1.0;
      for (ig=0; ig<gx; ig++) {
	for (iy=0; iy<my; iy++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix+ig,  iy,iz,mx,my);
	    int i_external = INDEX(ix-ig-1,iy,iz,mx,my);
	    array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && 
	     axis == axis_x) {
    if (nx > 1) {
      ix = nx+gx-1;
      sign = vx ? -1.0 : 1.0;
      for (ig=0; ig<gx; ig++) {
	for (iy=0; iy<my; iy++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix-ig,  iy,iz,mx,my);
	    int i_external = INDEX(ix+ig+1,iy,iz,mx,my);
	    array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_lower && 
	     axis == axis_y) {
    if (ny > 1) {
      iy = gy;
      sign = vy ? -1.0 : 1.0;
      for (ig=0; ig<gy; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,iy+ig,  iz,mx,my);
	    int i_external = INDEX(ix,iy-ig-1,iz,mx,my);
	    array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && 
	     axis == axis_y) {
    if (ny > 1) {
      iy = ny+gy-1;
      sign = vy ? -1.0 : 1.0;
      for (ig=0; ig<gy; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,iy-ig,  iz,mx,my);
	    int i_external = INDEX(ix,iy+ig+1,iz,mx,my);
	    array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_lower && 
	     axis == axis_z) {
    if (nz > 1) {
      iz = gz;
      sign = vz ? -1.0 : 1.0;
      for (ig=0; ig<gz; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iy=0; iy<my; iy++) {
	    int i_internal = INDEX(ix,iy,iz+ig,  mx,my);
	    int i_external = INDEX(ix,iy,iz-ig-1,mx,my);
	    array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && 
	     axis == axis_z) {
    if (nz > 1) {
      iz = nz+gz-1;
      sign = vz ? -1.0 : 1.0;
      for (ig=0; ig<gz; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iy=0; iy<my; iy++) {
	    int i_internal = INDEX(ix,iy,iz-ig,  mx,my);
	    int i_external = INDEX(ix,iy,iz+ig+1,mx,my);
	    array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  } else {
    ERROR("EnzoBoundary::enforce_reflecting_precision_",
	  "Cannot be called with face_all");
  }
}

//----------------------------------------------------------------------

void EnzoBoundary::enforce_outflow_
(
 const FieldDescr * field_descr,
 FieldBlock       * field_block,
 face_enum face,
 axis_enum axis
 ) const throw()
{

  int nx,ny,nz;
  int gx,gy,gz;
  field_block->size(&nx,&ny,&nz);
  for (int field = 0; field < field_descr->field_count(); field++) {
    field_descr->ghosts(field,&gx,&gy,&gz);
    void * array = field_block->field_values(field);
    precision_enum precision = field_descr->precision(field);
    switch (precision) {
    case precision_single:
      enforce_outflow_precision_(face,axis, (float *)array,
				 nx,ny,nz, gx,gy,gz);
      break;
    case precision_double:
      enforce_outflow_precision_(face,axis, (double *)array,
				 nx,ny,nz, gx,gy,gz);
      break;
    case precision_extended80:
    case precision_extended96:
    case precision_quadruple:
      enforce_outflow_precision_(face,axis, (long double *)array,
				 nx,ny,nz, gx,gy,gz);
      break;
    default:
      break;
    }
  }
}

//----------------------------------------------------------------------

template<class T>
void EnzoBoundary::enforce_outflow_precision_
(
 face_enum face, 
 axis_enum axis,
 T * array,
 int nx,int ny,int nz,
 int gx,int gy,int gz) const throw()
{
  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  int ix,iy,iz,ig;

  if (face == face_lower && 
      axis == axis_x) {
    if (nx > 1) {
      ix = gx;
      for (ig=0; ig<gx; ig++) {
	for (iy=0; iy<my; iy++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,     iy,iz,mx,my);
	    int i_external = INDEX(ix-ig-1,iy,iz,mx,my);
	    array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && 
	     axis == axis_x) {
    if (nx > 1) {
      ix = nx+gx-1;
      for (ig=0; ig<gx; ig++) {
	for (iy=0; iy<my; iy++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,     iy,iz,mx,my);
	    int i_external = INDEX(ix+ig+1,iy,iz,mx,my);
	    array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_lower && 
	     axis == axis_y) {
    if (ny > 1) {
      iy = gy;
      for (ig=0; ig<gy; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,iy,     iz,mx,my);
	    int i_external = INDEX(ix,iy-ig-1,iz,mx,my);
	    array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && 
	     axis == axis_y) {
    if (ny > 1) {
      iy = ny+gy-1;
      for (ig=0; ig<gy; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,iy,     iz,mx,my);
	    int i_external = INDEX(ix,iy+ig+1,iz,mx,my);
	    array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_lower && 
	     axis == axis_z) {
    if (nz > 1) {
      iz = gz;
      for (ig=0; ig<gz; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iy=0; iy<my; iy++) {
	    int i_internal = INDEX(ix,iy,iz,     mx,my);
	    int i_external = INDEX(ix,iy,iz-ig-1,mx,my);
	    array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && 
	     axis == axis_z) {
    if (nz > 1) {
      iz = nz+gz-1;
      for (ig=0; ig<gz; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iy=0; iy<my; iy++) {
	    int i_internal = INDEX(ix,iy,iz,     mx,my);
	    int i_external = INDEX(ix,iy,iz+ig+1,mx,my);
	    array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else {
    ERROR("EnzoBoundary::enforce_outflow_precision_",
	  "Cannot be called with face_all");
  }
}

//----------------------------------------------------------------------

void EnzoBoundary::enforce_inflow_
(
 const FieldDescr * field_descr,
 FieldBlock * field_block,
 face_enum face,
 axis_enum axis
 ) const throw()
{
  INCOMPLETE("EnzoBoundary::enforce_inflow");
}

//----------------------------------------------------------------------

