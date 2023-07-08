// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBoundary.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoBoundary class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoBoundary::EnzoBoundary 
(axis_enum axis, face_enum face, std::shared_ptr<Mask> mask,
 boundary_type type) throw()
  : Boundary(axis,face,mask),
    boundary_type_ (type)
{  }

//----------------------------------------------------------------------

void EnzoBoundary::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Boundary::pup(p);

  p | boundary_type_;

}

//----------------------------------------------------------------------

void EnzoBoundary::enforce 
(
 Block   * block,
 face_enum face,
 axis_enum axis 
 ) const throw()
{
  if ( ! applies_(axis,face)) {
    return;
  }

  if (face == face_all) {
    enforce(block,face_lower,axis);
    enforce(block,face_upper,axis);
  } else if (axis == axis_all) {
    enforce(block,face,axis_x);
    enforce(block,face,axis_y);
    enforce(block,face,axis_z);
  } else {

    Data * data = block->data();

    Field field = data->field();

    if ( ! field.ghosts_allocated() ) {
      ERROR("EnzoBoundary::enforce",
	    "Function called with ghosts not allocated");
    }

    switch (boundary_type_) {
    case boundary_type_reflecting:
      enforce_reflecting_(field,block,face,axis);
      break;
    case boundary_type_outflow:
      enforce_outflow_   (field,block,face,axis);
      break;
    default:
      ERROR("EnzoBoundary::enforce",
	    "Undefined boundary type");
      break;
    }
  }
}

//----------------------------------------------------------------------
bool EnzoBoundary::has_vector_name_(std::string field_name,
				    std::string component) const throw()
{
  return (field_name == ("bfield_"  + component) ||
	  field_name == ("bfieldi_" + component) ||
	  field_name == ("velocity_" + component));
}


//----------------------------------------------------------------------

void EnzoBoundary::enforce_reflecting_
(
 Field     field,
 Block   * block,
 face_enum face,
 axis_enum axis
 ) const throw()
{

  Data * data = block->data();

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  double * xc = new double [nx];
  double * yc = new double [ny];
  double * zc = new double [nz];

  data->field_cells(xc,yc,zc);

  double * xf = new double [nx+1];
  double * yf = new double [ny+1];
  double * zf = new double [nz+1];

  data->field_cell_faces(xf,yf,zf);

  double xm,ym,zm;
  double xp,yp,zp;
  data -> lower(&xm,&ym,&zm);
  data -> upper(&xp,&yp,&zp);

  double t = block->time();

  const int rank = cello::rank();

  // @@@
  // @@@ BUG: loops through all fields; should only use fields in field_list
  // @@@
  for (int index = 0; index < field.field_count(); index++) {
    if (field.is_temporary(index)){
      continue;
    }
    int gx,gy,gz;
    field.ghost_depth(index,&gx,&gy,&gz);
    int cx,cy,cz;
    field.centering(index, &cx,&cy,&cz);
    double * x = (cx == 1) ? xf : xc;
    double * y = (cy == 1) ? yf : yc;
    double * z = (cz == 1) ? zf : zc;
 
    enzo_float * array = (enzo_float * ) field.values(index);
    bool vx = (rank >= 1) && has_vector_name_(field.field_name(index), "x");
    bool vy = (rank >= 2) && has_vector_name_(field.field_name(index), "y");
    bool vz = (rank >= 3) && has_vector_name_(field.field_name(index), "z");
    enforce_reflecting_precision_(face,axis, array,
				  nx,ny,nz, gx,gy,gz, cx,cy,cz, vx,vy,vz,
				  x,y,z,    xm,ym,zm, xp,yp,zp, t);
  }

  delete [] xc;
  delete [] yc;
  delete [] zc;
  delete [] xf;
  delete [] yf;
  delete [] zf;
}

//----------------------------------------------------------------------

void EnzoBoundary::enforce_reflecting_precision_
(
 face_enum face, 
 axis_enum axis,
 enzo_float * array,
 int nx,int ny,int nz,
 int gx,int gy,int gz,
 int cx,int cy,int cz,
 bool vx,bool vy,bool vz,
 double *x, double *y, double *z,
 double xm, double ym, double zm,
 double xp, double yp, double zp,
 double t
 ) const throw()
{
  int mx = nx + 2*gx + cx;
  int my = ny + 2*gy + cy;
  int mz = nz + 2*gz + cz;

  int ix,iy,iz,ig;
  enzo_float sign;

  if (nx > 1) {
    if (face == face_lower && axis == axis_x) {
      ix = gx+cx;
      sign = vx ? -1.0 : 1.0;
      for (ig=0; ig<gx; ig++) {
	for (iy=0; iy<my; iy++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix+ig,     iy,iz,mx,my);
	    int i_external = INDEX(ix-ig-1-cx,iy,iz,mx,my);
	    if (! mask_ || mask_->evaluate(t,xm,y[iy],z[iz]))
	      array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
    if (face == face_upper &&  axis == axis_x) {
      ix = nx+gx-1;
      sign = vx ? -1.0 : 1.0;
      for (ig=0; ig<gx; ig++) {
	for (iy=0; iy<my; iy++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix-ig,     iy,iz,mx,my);
	    int i_external = INDEX(ix+ig+1+cx,iy,iz,mx,my);
	    if (! mask_ || mask_->evaluate(t,xp,y[iy],z[iz]))
	      array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  } 

  if (ny > 1) {
    if (face == face_lower && axis == axis_y) {
      iy = gy+cy;
      sign = vy ? -1.0 : 1.0;
      for (ig=0; ig<gy; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,iy+ig,     iz,mx,my);
	    int i_external = INDEX(ix,iy-ig-1-cy,iz,mx,my);
	    if (! mask_ || mask_->evaluate(t,x[ix],ym,z[iz]))
	      array[i_external] = sign*array[i_internal];
	  }
	}
      }
    } 
    if (face == face_upper && axis == axis_y) {
      iy = ny+gy-1;
      sign = vy ? -1.0 : 1.0;
      for (ig=0; ig<gy; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,iy-ig,     iz,mx,my);
	    int i_external = INDEX(ix,iy+ig+1+cy,iz,mx,my);
	    if (! mask_ || mask_->evaluate(t,x[ix],yp,z[iz]))
	      array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  }

  if (nz > 1) {
    if (face == face_lower && axis == axis_z) {
      iz = gz+cz;
      sign = vz ? -1.0 : 1.0;
      for (ig=0; ig<gz; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iy=0; iy<my; iy++) {
	    int i_internal = INDEX(ix,iy,iz+ig,     mx,my);
	    int i_external = INDEX(ix,iy,iz-ig-1-cz,mx,my);
	    if (! mask_ || mask_->evaluate(t,x[ix],y[iy],zm))
	      array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
    if (face == face_upper && axis == axis_z) {
      iz = nz+gz-1;
      sign = vz ? -1.0 : 1.0;
      for (ig=0; ig<gz; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iy=0; iy<my; iy++) {
	    int i_internal = INDEX(ix,iy,iz-ig,     mx,my);
	    int i_external = INDEX(ix,iy,iz+ig+1+cz,mx,my);
	    if (! mask_ || mask_->evaluate(t,x[ix],y[iy],zp))
	      array[i_external] = sign*array[i_internal];
	  }
	}
      }
    }
  }

  if (face == face_all) {
    ERROR("EnzoBoundary::enforce_reflecting_precision_",
	  "Cannot be called with face_all");
  }

}

//----------------------------------------------------------------------

void EnzoBoundary::enforce_outflow_
(
 Field     field,
 Block   * block,
 face_enum face,
 axis_enum axis
 ) const throw()
{

  Data * data = block->data();

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  double * xc = new double [nx];
  double * yc = new double [ny];
  double * zc = new double [nz];

  data->field_cells(xc,yc,zc);

  double * xf = new double [nx+1];
  double * yf = new double [ny+1];
  double * zf = new double [nz+1];

  data->field_cell_faces(xf,yf,zf);

  double xm,ym,zm;
  double xp,yp,zp;
  data -> lower(&xm,&ym,&zm);
  data -> upper(&xp,&yp,&zp);

  double t = block->time();

  // @@@
  // @@@ BUG: loops through all fields; should only use fields in field_list
  // @@@

  for (int index = 0; index < field.field_count(); index++) {
    if (field.is_temporary(index)){
      continue;
    }
    
    int gx,gy,gz;
    field.ghost_depth(index,&gx,&gy,&gz);

    int cx,cy,cz;
    field.centering(index, &cx,&cy,&cz);
    double * x = (cx == 1) ? xf : xc;
    double * y = (cy == 1) ? yf : yc;
    double * z = (cz == 1) ? zf : zc;
    
    enzo_float * array = (enzo_float * ) field.values(index);

    enforce_outflow_precision_(face,axis, array,
			       nx,ny,nz, gx,gy,gz, cx,cy,cz,
			       x,y,z,    xm,ym,zm, xp,yp,zp, t);
    
  }
  delete [] xc;
  delete [] yc;
  delete [] zc;

  delete [] xf;
  delete [] yf;
  delete [] zf;
}

//----------------------------------------------------------------------

void EnzoBoundary::enforce_outflow_precision_
(
 face_enum face, 
 axis_enum axis,
 enzo_float * array,
 int nx,int ny,int nz,
 int gx,int gy,int gz,
 int cx,int cy,int cz,
 double * x, double * y, double * z,
 double xm, double ym, double zm,
 double xp, double yp, double zp,
 double t
) const throw()
{
  int mx = nx + 2*gx + cx;
  int my = ny + 2*gy + cy;
  int mz = nz + 2*gz + cz;

  int ix,iy,iz,ig;

  if (face == face_lower && axis == axis_x) {
    if (nx > 1) {
      ix = gx;
      for (ig=0; ig<gx; ig++) {
	for (iy=0; iy<my; iy++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,     iy,iz,mx,my);
	    int i_external = INDEX(ix-ig-1,iy,iz,mx,my);
	    if (! mask_ || mask_->evaluate(t,xm,y[iy],z[iz]))
	      array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && axis == axis_x) {
    if (nx > 1) {
      ix = nx+gx-1+cx;
      for (ig=0; ig<gx; ig++) {
	for (iy=0; iy<my; iy++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,     iy,iz,mx,my);
	    int i_external = INDEX(ix+ig+1,iy,iz,mx,my);
	    if (! mask_ || mask_->evaluate(t,xp,y[iy],z[iz]))
	      array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_lower && axis == axis_y) {
    if (ny > 1) {
      iy = gy;
      for (ig=0; ig<gy; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,iy,     iz,mx,my);
	    int i_external = INDEX(ix,iy-ig-1,iz,mx,my);
	    if (! mask_ || mask_->evaluate(t,x[ix],ym,z[iz]))
	      array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && axis == axis_y) {
    if (ny > 1) {
      iy = ny+gy-1+cy;
      for (ig=0; ig<gy; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iz=0; iz<mz; iz++) {
	    int i_internal = INDEX(ix,iy,     iz,mx,my);
	    int i_external = INDEX(ix,iy+ig+1,iz,mx,my);
	    if (! mask_ || mask_->evaluate(t,x[ix],yp,z[iz]))
	      array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_lower && axis == axis_z) {
    if (nz > 1) {
      iz = gz;
      for (ig=0; ig<gz; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iy=0; iy<my; iy++) {
	    int i_internal = INDEX(ix,iy,iz,     mx,my);
	    int i_external = INDEX(ix,iy,iz-ig-1,mx,my);
	    if (! mask_ || mask_->evaluate(t,x[ix],y[iy],zm))
	      array[i_external] = array[i_internal];
	  }
	}
      }
    }
  } else if (face == face_upper && axis == axis_z) {
    if (nz > 1) {
      iz = nz+gz-1+cz;
      for (ig=0; ig<gz; ig++) {
	for (ix=0; ix<mx; ix++) {
	  for (iy=0; iy<my; iy++) {
	    int i_internal = INDEX(ix,iy,iz,     mx,my);
	    int i_external = INDEX(ix,iy,iz+ig+1,mx,my);
	    if (! mask_ || mask_->evaluate(t,x[ix],y[iy],zp))
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

