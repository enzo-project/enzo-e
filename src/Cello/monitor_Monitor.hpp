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

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change

    p |  *timer_;
    p |  active_;
    WARNING("Monitor::pup","ip_ may change");
    p |  ip_;
    p |  group_default_;
    p |  group_active_;

  }
#endif

  /// Return an instance of a Monitor object
  static Monitor * instance()
  { 
    return & instance_;
  };

  /// Set the processor's rank 
  /// (note cannot use GroupProcess since Monitor object created at startup)
  void set_process_rank (int ip)
  { ip_ = ip; }

  /// Set whether the monitor is active for text output.  Useful for
  /// parallel, e.g. "monitor->set_active(parallel->is_root())"
  void set_active(bool active) { active_ = active; };

  /// Return whether monitoring is active
  bool is_active() const throw () { return active_; };

  /// Set whether to monitor for the given component
  void set_active(const char * component, bool active) 
  { group_active_[component] = active; };

  /// Return whether monitoring is active for this component
  bool is_active(const char *) const throw ();

  /// Print the Cello header 
  void header () const;

  /// Write a message to file
  void write (FILE * fp, 
	      const char * component, const char * buffer, ...) const;
  /// Write a message to file
  void write_verbatim (FILE * fp, 
	      const char * component, const char * buffer) const;

  /// Print a message with possible format specifications to stdout
  void print (const char * component, const char * buffer, ...) const;

  /// Print a message without format specifications to stdout
  void print_verbatim (const char * component, const char * buffer) const;


private: // functions

  //----------------------------------------------------------------------

private: // attributes

  /// Timer for keeping track of time for output
  Timer * timer_; 

  /// Whether monitoring is activated.  Used for e.g. ip != 0.
  bool active_;

  /// Owning processor's rank
  int ip_;

  /// Whether default is to output all groups or output no groups
  bool group_default_;

  /// Override default of group_active_ for specific groups
  std::map<std::string,bool> group_active_;


  //----------------------------------------------------------------------

private: // static attributes

  /// Single instance of the Monitor object [singleton design pattern]
  // static Monitor * instance_;
  static Monitor instance_;

};

#endif /* MONITOR_MONITOR_HPP */

