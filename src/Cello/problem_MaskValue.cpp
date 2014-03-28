// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MaskValue.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    Implementation of the MaskValue class

#include "problem.hpp"

//----------------------------------------------------------------------
#ifdef BLAH_BLAH

MaskValue::MaskValue(Parameters * parameters,
		     std::string parameter_name,
		     ) throw ()
  : parameters_(parameters)
{

  if (parameters_->type(parameter_name) == parameter_list) {
    int num_values = parameters_->list_length(parameter_name);
    if (num_values > 1) {
      num_masks_ =  num_values / 2;
      mask_ = new bool * [num_masks_];
      nx_   = new int [num_masks_];
      ny_   = new int [num_masks_];
      // loop through value masks
      for (int index_mask=0; 
	   index_mask < num_masks_;
	   index_mask++) {
	int index_value = index_mask*2+1;
	if (parameters_->list_type
	    (index_value,parameter_name) == parameter_string) {
	  std::string file = parameters_->list_value_string
	    (index_value,parameter_name,"default");
	  create_mask_ (&mask_[index_mask],
			&nx_[index_mask],
			&ny_[index_mask],
			file);
			   
	  
	}
      }
    }
  }
}

//----------------------------------------------------------------------

MaskValue::~MaskValue() throw ()
{
  INCOMPLETE("MaskValue::~MaskValue");
}

//----------------------------------------------------------------------

MaskValue::MaskValue(const MaskValue & mask) throw ()
/// @param     MaskValue  Object being copied
{
  INCOMPLETE("MaskValue::MaskValue(MaskValue)");
}

//----------------------------------------------------------------------

MaskValue & MaskValue::operator= (const MaskValue & mask) throw ()
/// @param     MaskValue  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE("MaskValue::operator=");
  return *this;
}

void MaskValue::evaluate
  (const Hierarchy * hierarchy,
   const CommBlock * comm_block,
   FieldBlock * field_block, 
   int index_field,  int index_value,
   std::string field_name, 
   const FieldDescr * field_descr,
   int n, bool * mask, bool * deflt,
   double * x, double * y, double * z, double t) throw ()
{
  // parameter: Initial : <field> : value

  int index_mask = index_value / 2;

  parameter_type value_type = 
    parameters_->list_type(index_value,"value");

  bool v;
  std::string file;
  int nxb,nyb,nzb;

  switch (value_type) {
  case parameter_logical:

    // logical constant

    v = parameters_->list_value_logical(index_value,"value",false);
    for (int i=0; i<n; i++) mask[i] = v;
    break;

  case parameter_logical_expr:

    // logical expression

    parameters_->list_evaluate_logical
      (index_value,"value",n,mask,deflt,x,y,z,t);
    break;

  case parameter_string:

    {
      field_block->size(&nxb,&nyb,&nzb);
      ASSERT1("InitialValue::evaluate_logical",
	      "mask file %s requires problem to be 2D",
	      field_name.c_str(),
	      nyb > 1 && nzb == 1);

      bool * mask_png = mask_[index_mask];
      int nx_png = nx_[index_mask];
      int ny_png = ny_[index_mask];
      evaluate_mask_png_
	(mask,nxb,nyb,
	 mask_png,nx_png, ny_png,
	 hierarchy,comm_block,field_descr);
    }
    break;

  default:
    ERROR3("InitialValue::evaluate_mask",
	   "Even-index element %d of %s is of illegal type %d",
	   index_value,field_name.c_str(),value_type);
    break;
  }
}
//======================================================================

#endif
