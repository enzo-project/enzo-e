// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     monitor_Monitor.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-05
/// @brief    [\ref Monitor] Declaration of the Monitor class
//----------------------------------------------------------------------

#ifndef MONITOR_MONITOR_HPP
#define MONITOR_MONITOR_HPP

//----------------------------------------------------------------------
/// @def    MONITOR_LENGTH 255
/// @brief  Maximum length of monitor text output

#define MONITOR_LENGTH 255
   
//----------------------------------------------------------------------
class Timer; 

class Monitor {

  /// @class    Monitor
  /// @ingroup  Monitor
  /// @todo     Make calling image() easier
  /// @todo     Add configuration parameters, e.g. monitoring level, date/time output, format, etc.
  /// @todo     Add calling component to print
  /// @brief    [\ref Monitor] User monitoring of simulation execution status
  ///
  /// The Monitor component is used to communicate information about
  /// the running simulation to the user. Information can be output in
  /// several forms, including text files, HTML files, plots, or other
  /// (generally small) image files. Information is assumed to be from
  /// a correctly-running simulation: anomalous errors or warnings are
  /// output by the Error component.

  //----------------------------------------------------------------------

private:

  friend class Simulation;

  /// Private constructor of the Monitor object [singleton design pattern]
  Monitor();

  /// Private destructor  of the Monitor object [singleton design pattern]
  ~Monitor();

//----------------------------------------------------------------------

public: // interface

  /// Return an instance of a Monitor object
  static Monitor * instance()
  { 
    return & instance_;
    // if ( instance_ == NULL )
    //   instance_ = new Monitor;

    // return instance_;
  };

  /// Set whether the monitor is active for text output.  Useful for
  /// parallel, e.g. "monitor->set_active(parallel->is_root())"
  void set_active(bool active) { active_ = active; };

  /// Return whether monitoring is active
  bool is_active() const throw () { return active_; };

  /// Print the Cello header 
  void header () const;

  /// Print a message to stdout
  void print (const char * component, const char * buffer, ...) const;

  /// Write a message to file
  void write (FILE * fp, 
	      const char * component, const char * buffer, ...) const;

private: // functions


  //----------------------------------------------------------------------

private: // attributes

  /// Timer for keeping track of time for output
  Timer * timer_; 

  /// Whether monitoring is activated.  Used for e.g. ip != 0.
  bool active_;

  //----------------------------------------------------------------------

private: // static attributes

  /// Single instance of the Monitor object [singleton design pattern]
  // static Monitor * instance_;
  static Monitor instance_;

};

#endif /* MONITOR_MONITOR_HPP */

