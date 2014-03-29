// See LICENSE_CELLO file for license and copyright information

/// @file     method_InitialValue.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Method] Default initialization method

#ifndef METHOD_INITIAL_DEFAULT_HPP
#define METHOD_INITIAL_DEFAULT_HPP

class InitialValue : public Initial {

  /// @class    InitialValue
  /// @ingroup  Method
  /// @brief    [\ref Method] Default initialization method

public: // interface

  /// CHARM++ constructor
  InitialValue() throw() { }
  
  /// Constructor
  InitialValue(Parameters * parameters, 
		 const FieldDescr * field_descr,
		 int cycle, double time) throw();

  /// Destructor
  virtual ~InitialValue() throw();

  PUPable_decl(InitialValue);

  InitialValue(CkMigrateMessage *m) : Initial (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Read initialization values from Initial group in parameter file

  virtual void enforce_block (CommBlock * block,
			      const FieldDescr * field_descr,
			      const Hierarchy * hierarchy
			      ) throw();

private: // functions
  
  void allocate_xyzt_(CommBlock * block,
		      int index_field,
		      const FieldBlock * field_block,
		      const FieldDescr * field_descr,
		      int * mx, int * my, int * mz,
		      double ** value, double ** vdeflt,
		      bool ** mask, bool ** rdeflt,
		      double ** x, double ** y, double ** z,
		      double * t) throw();

  void copy_values_ (
		     const FieldDescr * field_descr,
		     FieldBlock * field_block,
		     double * value, bool * mask,
		     int index_field,
		     int nx, int ny, int nz) throw();

  void evaluate_float_ (FieldBlock * field_block, int index_field, 
			std::string field_name,
			int n, double * value, double * vdeflt,
			double * x, double * y, double * z, double t) throw();

  void evaluate_mask_ (const Hierarchy * hierarchy,
		       const CommBlock * block,
		       FieldBlock * field_block,
		       int index_field, int index_value,
		       std::string field_name,
		       const FieldDescr * field_descr,			
		       int n, bool * value, bool * vdeflt,
		       double * x, double * y, double * z, double t) throw();

  /// Read in a PNG file and create an integer array using r + b + g values
  void create_mask_png_ (bool ** mask, int * nx, int * ny ,
			 std::string pngfile);

  /// Create a mask for the block given mask_[][]
  void evaluate_mask_png_ ( bool            * mask_block, int nxb, int nyb,
			    bool            * mask_png,   int nx,  int ny,
			    const Hierarchy * hierarchy,
			    const CommBlock * comm_block,
			    const FieldDescr * field_descr);

  template<class T>
  void copy_precision_
  (T * field, bool * mask, int offset, double * value, int nx, int ny, int nz);

  template<class T>
  void copy_precision_
  (T * field, int offset, double * value, int nx, int ny, int nz);

private: // attributes

  Parameters * parameters_;

  /// Field descriptor
  const FieldDescr * field_descr_;

  /// number of fields
  int num_fields_;

  /// number of masked values per field
  int * num_masks_;

  /// Masks for fields and values: mask_[index_field][index_value]
  bool *** mask_;

  Mask *** mask_list_;

  /// Size of the masks
  int **nx_;
  int **ny_;

};

#endif /* METHOD_INITIAL_DEFAULT_HPP */

