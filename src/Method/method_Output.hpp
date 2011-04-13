// $Id: method_Output.hpp 1896 2010-12-03 23:54:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Output.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @todo     Create and use Schedule object for scheduling output
/// @brief    [\ref Method] Declaration for the Output component

#ifndef METHOD_OUTPUT_HPP
#define METHOD_OUTPUT_HPP

enum output_schedule_enum {
  output_schedule_unknown,
  output_schedule_cycle_interval,
  output_schedule_cycle_list,
  output_schedule_time_interval,
  output_schedule_time_list
};

class Output {

  /// @class    Output
  /// @ingroup  Method
  /// @brief    [\ref Method] define interface for various types of simulation output

public: // functions

  /// Create an uninitialized Output object with the given file_name format
  Output() throw();

  /// Set file_name
  void set_file_name (std::string file_name) throw()
  { file_name_ = file_name; };

  /// Set cycle interval (start, step, stop)
  void set_cycle_interval
  (int cycle_start, int cycle_step, int cycle_stop) throw();

  /// Set cycle list
  void set_cycle_list (std::vector<int> cycle_list) throw();

  /// Set field list
  void set_field_list (std::vector<int> field_list) throw();

  /// Set time interval (start, step, stop)
  void set_time_interval
  (double time_start, double time_step, double time_stop) throw();

  /// Set time list
  void set_time_list (std::vector<double> time_list) throw();

  /// Set whether the Output object is active or not
  void set_active(bool active) throw()
  { active_ = active; };
    
  /// Reduce timestep if next write time is greater than time + dt
  double update_timestep(double time, double dt) const throw();

  /// Return true if output should be performed this cycle.
  bool write_this_cycle(int cycle, double time) throw();

  /// Return whether the output object is active
  bool is_active() const throw()
  { return active_; };

  /// Write mesh-related data to disk if scheduled
  void scheduled_write
  (  const FieldDescr * field_descr,
     Mesh * mesh, int cycle, double time, bool root_call=true) throw();

  /// Write a patch-related data to disk if scheduled
  void scheduled_write
  (  const FieldDescr * field_descr,
     Patch * patch, Mesh * mesh,
     int cycle, double time, bool root_call=true) throw();

  /// Write a block-related to disk if scheduled
  void scheduled_write
  (  const FieldDescr * field_descr,
     Block * block, Patch * patch, Mesh * mesh,
     int cycle, double time, bool root_call=true) throw();

  std::string expand_file_name (int cycle, double time) const throw();

public: // virtual functions

  /// Write mesh data to disk
  virtual void write 
  ( const FieldDescr * field_descr,
    int index, Mesh * mesh, 
    int cycle, double time,
    bool root_call=true, int ix0=0, int iy0=0, int iz0=0) throw() = 0;

  /// Write patch data to disk; may be called by write (Mesh)
  virtual void write 
  ( const FieldDescr * field_descr,
    int index, Patch * patch, Mesh * mesh,
    int cycle, double time, 
    bool root_call=true, int ix0=0, int iy0=0, int iz0=0) throw() = 0;

  /// Write block data to disk; may be called by write (Patch)
  virtual void write 
  ( const FieldDescr * field_descr,
    int index, Block * block, Patch * patch, Mesh * mesh, 
    int cycle, double time, 
    bool root_call=true, int ix0=0, int iy0=0, int iz0=0) throw() = 0;

protected: // attributes

  /// Name of the file to write, including printf-type format
  std::string file_name_;

  /// Whether Output is currently active
  bool active_;

  /// Whether Output is scheduled for next call to scheduled write
  bool scheduled_;
  
  /// Schedule type of the Output object
  output_schedule_enum output_schedule_;

  /// cycle start, step, and stop of output
  std::vector<int> cycle_interval_;

  /// List of cycles to perform output
  std::vector<int> cycle_list_;

  /// time start, step, and stop of output
  std::vector<double> time_interval_;

  /// List of times to perform output
  std::vector<double> time_list_;

  /// List of fields to output
  std::vector<int> field_list_;

  /// Index of time or cycle interval or list for next output
  size_t index_;

};

#endif /* METHOD_OUTPUT_HPP */
