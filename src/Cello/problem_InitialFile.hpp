// See LICENSE_CELLO file for license and copyright information

/// @file     method_InitialFile.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Problem] Declaration of the InitialFile class
///
/// 

#ifndef METHOD_INITIAL_FILE_HPP
#define METHOD_INITIAL_FILE_HPP

class InitialFile : public Initial {

  /// @class    InitialFile
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Declaration of the InitialFile class
  ///
  /// This class is used to define initial conditions by reading in
  /// data from files


public: // interface

  /// Constructor
  InitialFile(Parameters * parameters, int cycle, double time) throw();

  /// Destructor
  virtual ~InitialFile() throw();

  /// Read initialization values from Initial group in parameter file

  virtual void enforce (const Hierarchy * hierarchy,
			const FieldDescr * field_descr,
			Block * block) throw();

private: // functions

  void get_filename_(std::string * file_name,
		     std::vector<std::string> * file_args) throw();

private: // attributes

  Parameters * parameters_;

  Input * input_;
};

#endif /* METHOD_INITIAL_FILE_HPP */

