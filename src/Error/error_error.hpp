// $Id: error_error.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ERROR_ERROR_HPP
#define ERROR_ERROR_HPP

/// @file     error_error.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Add Parallel support
/// @todo     Use singleton design pattern
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the Error class

class Error {

  /// @class    Error
  /// @ingroup  Error
  /// @brief    Singleton class for reporting errors and warnings

public: // functions

  /// Get single instance of the Error object
  static Error * instance() throw ()
  { return & instance_; };

  /// Warning message

  void warning (const char * file,
		int          line,
		const char * function,
		const char * message)
  {
    message_("WARNING",file,line,function,message);
  };

  /// Incomplete message

  void incomplete (const char * file,
		   int          line,
		   const char * function,
		   const char * message)
  {
    message_("INCOMPLETE",file,line,function,message);
  };

  /// Error message

  void error (const char * file,
	      int          line,
	      const char * function,
	      const char * message)
  {
    message_("ERROR",file,line,function,message);
    exit(1);
  };

  void trace (const char * file,
	      int          line)
  {
    if (trace_) {
      printf ("TRACE %s:%d\n",__FILE__,__LINE__); 
      fflush(stdout);
    }
  };

  void tracing (bool trace) 
  { trace_ = trace; };

private: // functions

  /// Initialize the Error object (singleton design pattern)
  Error() 
    : trace_(true)
  {};

  /// Delete the Error object (singleton design pattern)
  ~Error()
  {};

private: // functions

  void message_  
  (
   const char * type,
   const char * file,
   int          line,
   const char * function,
   const char * message)
  {
    printf ("\n");
    printf ("     %s File:     %s:%d\n",type,file,line);
    printf ("     %s Function: %s()\n", type,function);
    printf ("     %s Message:  %s\n",   type,message);
    printf ("\n");
    fflush(stdout);
  };

private: // attributes

  /// Single instance of the Error object (singleton design pattern)
  static Error instance_;

  /// Whether tracing is on or off
  bool trace_;

};

/// @def      ERROR_MESSAGE_LENGTH
/// @brief    Maximum length of error and warning messages
#define ERROR_MESSAGE_LENGTH 255

/// @def      WARNING_MESSAGE
/// @brief    Handle a (non-lethal) warning message
#define WARNING_MESSAGE(FUNCTION,MESSAGE) \
  Error::instance()->warning(__FILE__,__LINE__,FUNCTION,MESSAGE)

/// @def      ERROR_MESSAGE
/// @brief    Handle a (lethal) error message
#define ERROR_MESSAGE(FUNCTION,MESSAGE) \
  Error::instance()->error(__FILE__,__LINE__,FUNCTION,MESSAGE)

/// @def      INCOMPLETE_MESSAGE
/// @brief    Placeholder for code that is incomplete
#define INCOMPLETE_MESSAGE(FUNCTION,MESSAGE) \
  Error::instance()->incomplete(__FILE__,__LINE__,FUNCTION,MESSAGE)

/// @def      TRACE
/// @brief    Trace file name and location to stdout
#define TRACE					\
  Error::instance()->trace(__FILE__,__LINE__)


#endif /* ERROR_ERROR_HPP */

