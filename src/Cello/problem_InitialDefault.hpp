// See LICENSE_CELLO file for license and copyright information

/// @file     method_InitialDefault.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Method] Default initialization method

#ifndef METHOD_INITIAL_DEFAULT_HPP
#define METHOD_INITIAL_DEFAULT_HPP

class InitialDefault : public Initial {

  /// @class    InitialDefault
  /// @ingroup  Method
  /// @brief    [\ref Method] Default initialization method

public: // interface

  /// Constructor
  InitialDefault(Parameters * parameters, int cycle, double time) throw();

  /// Destructor
  virtual ~InitialDefault() throw()
  {}

  /// Read initialization values from Initial group in parameter file

  virtual void enforce (Block * block,
			const FieldDescr * field_descr,
			const Hierarchy * hierarchy
			) throw();

private: // functions
  
  void allocate_xyzt_(Block * block,
		      int index_field,
		      int * mx, int * my, int * mz,
		      double ** value, double ** vdeflt,
		      bool ** mask, bool ** rdeflt,
		      double ** x, double ** y, double ** z,
		      double ** t) throw();

  void copy_values_ (
		     const FieldDescr * field_descr,
		     FieldBlock * field_block,
		     double * value, bool * mask,
		     int index_field,
		     int nx, int ny, int nz) throw();

  void evaluate_float_ (FieldBlock * field_block, int index_field, 
			std::string field_name,
			int n, double * value, double * vdeflt,
			double * x, double * y, double * z, double * t) throw();

  void evaluate_logical_ (const Hierarchy * hierarchy,
			  const Block * block,
			  FieldBlock * field_block, int index_field, 
			 std::string field_name,
			 int n, bool * value, bool * vdeflt,
			 double * x, double * y, double * z, double * t) throw();

  /// Read in a PNG file and create an integer array using r + b + g values
  void create_png_mask_ (bool * mask,
			 const Hierarchy * hierarchy,
			 const Block * block,
			 const char * pngfile,
			 int   nxb, int nyb,
			 int * nx, int * ny );

private: // attributes

  Parameters * parameters_;

};

#endif /* METHOD_INITIAL_DEFAULT_HPP */

