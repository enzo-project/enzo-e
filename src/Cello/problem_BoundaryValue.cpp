// See LICENSE_CELLO file for license and copyright information

/// @file     problem_BoundaryValue.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-04-02
/// @brief    Implementation of the default BoundaryValue boundary value class

#include "problem.hpp"

//----------------------------------------------------------------------

std::string BoundaryValue::debug_string() const throw(){
  // this function is very inefficient! But it's useful to return a string
  // while debugging
  auto vec_to_str = [](const std::vector<std::string>& vec) -> std::string
    {
      std::string out = "[";
      const std::size_t len = vec.size();
      for (std::size_t i = 0; i < len; i++) {
        out += ('"' + vec[i] + '"');
        if ((i+1) != len) { out += ", "; }
      }
      return out + ']';
    };

  std::string out = "BoundaryValue {\n";
  out            += ("  axis_: " + std::to_string(axis_) + ",\n");
  out            += ("  face_: " + std::to_string(face_) + ",\n");
  out            += ("  pairs_: [\n");
  for (std::size_t i = 0; i < pairs_.size(); i++) {
    const Value& value_obj = pairs_[i].first;
    const std::vector<std::string>& field_list = pairs_[i].second;
    out          += ("    {" + value_obj.debug_string() + ",\n");
    out          += ("     " + vec_to_str(field_list) + "}");
    if ((i + 1) != pairs_.size()) { out += ','; }
    out += '\n';
  }
  out            += ("  ]\n");
  return out;
}

//----------------------------------------------------------------------

void BoundaryValue::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change
  Boundary::pup(p); 
  TRACEPUP; 

  // historically there have been some issues with puping certain types of
  // Value objects. These issues are now always addressed within Value::pup
  p | pairs_;
}

//----------------------------------------------------------------------

static std::string ensure_trailing_colon_(const std::string &my_str)
{ return (my_str.back() == ':') ? my_str : (my_str + ':'); }

//----------------------------------------------------------------------

std::vector<BoundaryValue::ValueFListPair>
BoundaryValue::construct_ValueFList_pairs_
(Parameters& p, const std::string& parameter_group) throw()
{
  // param_param looks like: "Boundary:<boundary-name>:"
  const std::string param_prefix = ensure_trailing_colon_(parameter_group);

  const std::string value_param_str = param_prefix + "value";
  const std::string flist_param_str = param_prefix + "field_list";

  // pair_vec will hold a list of (value, field list) pairs
  std::vector<BoundaryValue::ValueFListPair> pair_vec;

  // in the following if-else statement, we append entries to pair_vec
  if (p.param(value_param_str) == nullptr) {
    // value_param_str isn't a parameter-name, so we retrieve a list of all
    // leaf parameters (if any) named: "<value_param_str>:<leaf_par>"
    // - all of these leaf parameters are assigned value-expressions and
    //   are named after the field that the expression applies to.
    std::vector<std::string> p_list = p.leaf_parameter_names(value_param_str);

    if ((p_list.size() > 0) && (p.param(flist_param_str) != nullptr)) {
      ERROR2("BoundaryValue::from_parameters",
             "Because \"%s\" is a parameter-group, the \"%s\" parameter is "
             "NOT allowed to be specified.",
             value_param_str.c_str(), flist_param_str.c_str());
    }

    pair_vec.reserve(p_list.size());
    const std::string grp_prefix = value_param_str + ":";
    for (const std::string& field_name : p_list) {
      std::vector<std::string> field_list = {field_name};
      pair_vec.push_back( { Value(&p, grp_prefix + field_name), field_list } );
    }

  } else {
    /// value_param_str is a parameter-name. Construct value object from it and
    /// load field_list parameter

    std::vector<std::string> field_list =
      p.value_full_strlist(flist_param_str, true,
                           true); // last arg says: return empty vec on error
    if (field_list.size() == 0) {
      ERROR2("BoundaryValue::from_parameters",
             "Since \"%s\" is a parameter with an asigned value-expression, "
             "\"%s\" must be assigned a field name or list of 1 or more field "
             "names.",
             value_param_str.c_str(), flist_param_str.c_str());
    }

    pair_vec.push_back( { Value(&p, value_param_str), field_list } );
  }

  if (pair_vec.size() == 0){
    ERROR1("BoundaryValue::from_parameters",
           "Error with \"%s\". It is expected to either be a parameter "
           "holding a value-expression or the name of a parameter-group, "
           "in which each parameter specifies a value-expression",
           value_param_str.c_str());
  }

  return pair_vec;
}

//----------------------------------------------------------------------

void BoundaryValue::enforce 
(Block * block, face_enum face, axis_enum axis) const throw()
{
  if ( ! applies_(axis,face)) return;

  if (face == face_all) {
    enforce(block,face_lower,axis);
    enforce(block,face_upper,axis);
    return;
  } else if (axis == axis_all) {
    enforce(block,face,axis_x);
    enforce(block,face,axis_y);
    enforce(block,face,axis_z);
    return;
  }

  Data * data = block->data();
  Field field = data->field();

  if ( ! field.ghosts_allocated() ) {
    ERROR("EnzoBoundary::enforce",
          "Function called with ghosts not allocated");
  }

  double xm,ym,zm;
  double xp,yp,zp;
  data -> lower(&xm,&ym,&zm);
  data -> upper(&xp,&yp,&zp);

  double t = block->state()->time();

  for (const BoundaryValue::ValueFListPair& cur_pair : pairs_) {
    const Value& value = cur_pair.first;
    const std::vector<std::string>& field_list = cur_pair.second;
    for (const std::string& field_name : field_list) {

      int nx,ny,nz;
      field.size(&nx,&ny,&nz);

      int index_field = field.field_id(field_name);
      int gx,gy,gz;
      field.ghost_depth(index_field,&gx,&gy,&gz);

      int cx,cy,cz;
      field.centering(index_field, &cx,&cy,&cz);

      int ndx=nx+2*gx+cx;
      int ndy=ny+2*gy+cy;
      int ndz=nz+2*gz+cz;

      double * x = new double [ndx];
      double * y = new double [ndy];
      double * z = new double [ndz];

      data->field_cell_faces(x,y,z,gx,gy,gz,cx,cy,cz);

      void * array = field.values(index_field);

      precision_type precision = field.precision(index_field);

      int ix0=0 ,iy0=0,iz0=0;

      nx = ndx;
      ny = ndy;
      nz = ndz;

      if (axis == axis_x) nx=gx;
      if (axis == axis_y) ny=gy;
      if (axis == axis_z) nz=gz;

      if (face == face_upper) {
	if (axis == axis_x) ix0 = ndx - gx;
	if (axis == axis_y) iy0 = ndy - gy;
	if (axis == axis_z) iz0 = ndz - gz;
      }

      int i0=ix0 + ndx*(iy0 + ndy*iz0);

      bool * mask = 0;

      if (mask_ != nullptr) mask = new bool [nx*ny*nz];

      switch (precision) {
      case precision_single:
	{
	  float * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (float *)array;
	    array = new float [ndx*ndy*ndz];
	  }
	  
	  value.evaluate((float *)array+i0, t,
                         ndx,nx,x+ix0,
                         ndy,ny,y+iy0,
                         ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) ((float *)temp)[i]=((float *)array)[i];
	    delete [] ((float*)array);
	    array = temp;
	  }
	}
       	break;
      case precision_double:
	{
	  double * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (double *)array;
	    array = new double [ndx*ndy*ndz];
	  }
	  value.evaluate((double *)array+i0, t,
                         ndx,nx,x+ix0,
                         ndy,ny,y+iy0,
                         ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) ((double *)temp)[i]=((double *)array)[i];
	    delete [] ((double *)array);
	    array = temp;
	  }
	}
       	break;
      case precision_extended80:
      case precision_extended96:
      case precision_quadruple:
	{
	  long double * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (long double *)array;
	    array = new long double [ndx*ndy*ndz];
	  }
	  value.evaluate((long double *)array+i0, t,
                         ndx,nx,x+ix0,
                         ndy,ny,y+iy0,
                         ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) 
	      ((long double *)temp)[i]=((long double *)array)[i];
	    delete [] ((long double *)array);
	    array = temp;
	  }
	}
       	break;
      }

      delete [] x;
      delete [] y;
      delete [] z;
      delete [] mask;
    } // for field_name in cur_pair.fields
  } // for cur_pair in pairs_
}

//----------------------------------------------------------------------

template <class T>
void BoundaryValue::copy_(T * field, double * value,
			  int ndx, int ndy, int ndz,
			  int nx,  int ny,  int nz,
			  int ix0, int iy0, int iz0) const throw()
{
  for (int ix=ix0; ix<ix0+nx; ix++) {
    for (int iy=iy0; iy<iy0+ny; iy++) {
      for (int iz=iz0; iz<iz0+nz; iz++) {
	int iv = (ix-ix0) + nx*((iy-iy0) + ny*(iz-iz0));
	int ib = ix + ndx*(iy + ndy*(iz));
	field[ib] = (T) value[iv];
      }
    }
  }
}

//----------------------------------------------------------------------
