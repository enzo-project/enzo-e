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

  /// CHARM++ constructor
  InitialFile() throw() { }

  /// Constructor
  InitialFile(Parameters * parameters, 
	      const GroupProcess * group_process,
	      int cycle, double time) throw();

  /// Destructor
  virtual ~InitialFile() throw();

#ifdef CONFIG_USE_CHARM

  PUPable_decl(InitialFile);

  InitialFile(CkMigrateMessage *m) : Initial (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#endif

  /// Read initialization values from Initial group in parameter file

  /// Enforce initial conditions for the given Block

  virtual void enforce (Block * block,
			const FieldDescr * field_descr,
			const Hierarchy * hierarchy) throw();

  /// Override default: InitialFile expects blocks not to be allocated
  virtual bool expects_blocks_allocated() const throw()
  { return false; }

private: // functions

  void get_filename_(std::string * file_name,
		     std::vector<std::string> * file_args) throw();

private: // attributes

  /// Parameters object
  Parameters * parameters_;

  /// Parallel GroupProcess object for creating Patches
  GroupProcess * group_process_;

  /// Associated Input object
  Input * input_;

#ifdef CONFIG_USE_CHARM
  /// Counter for reading blocks from current patch
  Loop block_loop_;
#endif
};

#endif /* METHOD_INITIAL_FILE_HPP */

