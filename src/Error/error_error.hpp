// $Id: error_error.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ERROR_ERROR_HPP
#define ERROR_ERROR_HPP

/// @file     error_error.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Add Parallel support
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the Error class

/// @def      ERROR_MESSAGE_LENGTH
/// @brief    Maximum length of error and warning messages
#define ERROR_MESSAGE_LENGTH 255

/// @def      WARNING_MESSAGE
/// @brief    Handle a (non-lethal) warning message
#define WARNING_MESSAGE(FUNCTION,MESSAGE) \
  Error::instance()->warning_(__FILE__,__LINE__,FUNCTION,MESSAGE)

/// @def      ERROR_MESSAGE
/// @brief    Handle a (lethal) error message
#define ERROR_MESSAGE(FUNCTION,MESSAGE) \
  Error::instance()->error_(__FILE__,__LINE__,FUNCTION,MESSAGE)

/// @def      INCOMPLETE_MESSAGE
/// @brief    Placeholder for code that is incomplete
#define INCOMPLETE_MESSAGE(FUNCTION,MESSAGE) \
  Error::instance()->incomplete_(__FILE__,__LINE__,FUNCTION,MESSAGE)

/// @def      TRACE
/// @brief    Trace file name and location to stdout
#define TRACE					\
  Error::instance()->trace_(__FILE__,__LINE__)

/// @def      ASSERT
/// @brief    Equivalent to assert()
#define ASSERT(FUNCTION,MESSAGE,ASSERTION) \
  Error::instance()->assert_(__FILE__,__LINE__,FUNCTION,MESSAGE,ASSERTION)

//----------------------------------------------------------------------

class Error {

  /// @class    Error
  /// @ingroup  Error
  /// @brief    Singleton class for reporting errors and warnings

public: // functions

  /// Get single instance of the Error object
  static Error * instance() throw ()
  { return & instance_; };

  /// Set whether to trace on this processor
  void set_tracing (bool tracing) 
  { tracing_ = tracing; };

public: // functions

  //----------------------------------------------------------------------
  /// Warning message
  void warning_ (const char * file,
		 int          line,
		 const char * function,
		 const char * message)
  {
    message_(stdout,"WARNING",file,line,function,message);
  };

  //----------------------------------------------------------------------
  /// Incomplete message
  void incomplete_ (const char * file,
		    int          line,
		    const char * function,
		    const char * message)
  {
    message_(stdout,"INCOMPLETE",file,line,function,message);
  };

  //----------------------------------------------------------------------
  /// Error message
  void error_ (const char * file,
	       int          line,
	       const char * function,
	       const char * message)
  {
    message_(stderr,"ERROR",file,line,function,message);
    exit(1);
  };

  //----------------------------------------------------------------------
  void trace_ (const char * file,
	       int          line)
  {
    if (tracing_) {
      printf ("TRACE %s:%d\n",file,line); 
      fflush(stdout);
    }
  };

  //----------------------------------------------------------------------
  void assert_ (const char * file,
		int          line,
		const char * function,
		const char * message,
		bool         assertion)
  {
    if (!assertion) {
      message_(stderr,"ASSERT",file,line,function,message);
      exit(1);
    }
  };

private: // functions

  //----------------------------------------------------------------------
  /// Initialize the Error object (singleton design pattern)
  Error() 
    : tracing_(true)
  {};

private: // functions

  void message_  
  (
   FILE *       fp,
   const char * type,
   const char * file,
   int          line,
   const char * function,
   const char * message)
  {
    fprintf (fp,"\n");
    fprintf (fp,"     %s File:     %s:%d\n",type,file,line);
    fprintf (fp,"     %s Function: %s()\n", type,function);
    fprintf (fp,"     %s Message:  %s\n",   type,message);
    fprintf (fp,"\n");
  };

private: // attributes

  /// Single instance of the Error object (singleton design pattern)
  static Error instance_;

  /// Whether tracing is on or off
  bool tracing_;

};


#endif /* ERROR_ERROR_HPP */

