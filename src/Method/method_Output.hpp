// $Id: method_Output.hpp 1896 2010-12-03 23:54:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Output.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Method] Declaration for the Output component

#ifndef METHOD_OUTPUT_HPP
#define METHOD_OUTPUT_HPP

class Output {

  /// @class    Output
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate functions for simulation output

public: // interface

  /// Create a new Output
  Output() throw()
  {};

public: // virtual functions

  /// Write projections of fields to png images
  virtual void write_image (Mesh * mesh) throw(), bool top_level=true;

  /// Write projections of fields to png images
  virtual void write_image (Patch * patch) throw(), bool top_level=true;

  /// Write projections of fields to png images
  virtual void write_image (Block * block) throw(), bool top_level=true;

  //--------------------------------------------------

  /// Write data dumps to disk
  virtual void write_data (Mesh * mesh, bool top_level=true) throw();

  /// Write data dumps to disk
  virtual void write_data (Patch * patch, bool top_level=true) throw();

  /// Write data dumps to disk
  virtual void write_data (Block * block, bool top_level=true) throw();

protected: // functions

  /// Return whether to write this cycle
  bool is_active_cycle_(int cycle) throw();

  /// Return whether to write this time
  bool is_active_time_(double time) throw();

protected: // attributes

  std::string filename_;

  /// cycle start,stop,step,start,stop,step,etc.  Empty if not active
  std::vector<int> cycle_range_;

  /// time start,stop,step,start,stop,step,etc.  Empty if not active
  std::vector<double> time_range_;
  
};

#endif /* METHOD_OUTPUT_HPP */
