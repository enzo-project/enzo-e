// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_CountersPapi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-20
/// @brief    Implementation of the CountersPapi class

#include "lcaperf.hpp"

#ifdef CONFIG_USE_PAPI

namespace lca {

const char * papi_counter_name[] = {
  "PAPI_FP_OPS",
};

//----------------------------------------------------------------------

CountersPapi::CountersPapi() throw ()
: CountersUser(),
  is_papi_active_(false),
  event_set_(PAPI_NULL),
  papi_counters_(NULL),
  vtime_begin_(0)
{
  create("papi-vtime-begin",  counter_type_absolute);
  create("papi-vtime-end",    counter_type_absolute);
  create("papi-vtime-active", counter_type_relative);
  create("papi-fp-ops",       counter_type_relative);

  // Initialize the PAPI library

  int retval;
  _CALL_PAPI(PAPI_library_init(PAPI_VER_CURRENT),"library_init",
	     PAPI_VER_CURRENT,retval);

  // Add PAPI counters

  papi_counters_ = new long long [num_papi_counters];

  for (int i=0; i<num_papi_counters; i++) {
    papi_insert_counter_((char *)papi_counter_name[i]);
  }
  papi_start_counters_();
}

//----------------------------------------------------------------------

CountersPapi::~CountersPapi() throw ()
{
  papi_stop_counters_();
  papi_delete_eventset_();
  delete [] papi_counters_;
  papi_counters_ = 0;
}

//----------------------------------------------------------------------

CountersPapi::CountersPapi(const CountersPapi & counters) throw ()
  : CountersUser(counters)
/// @param     counters  Object being copied
{
}

//----------------------------------------------------------------------

CountersPapi & CountersPapi::operator= (const CountersPapi & counters) throw ()
/// @param     counters  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}

//======================================================================

long long * CountersPapi::start_()
{
  long long * counters = CountersUser::start_();
  TRACE("void CountersPapi::start_");

#ifdef CONFIG_USE_PAPI
  vtime_begin_ = PAPI_get_virt_usec();
#else
  vtime_begin_ = 0;
#endif

  assign("papi-vtime-begin",  vtime_begin_);
  assign("papi-vtime-end",    0);
  assign("papi-vtime-active", 0);
  assign("papi-fp-ops",       0);
  return counters;

}

//----------------------------------------------------------------------

void CountersPapi::stop_(long long * counters)
{
  TRACE("void CountersPapi::stop_");
  papi_read_counters_();
#ifdef CONFIG_USE_PAPI
  long long vtime_end = PAPI_get_virt_usec();
#else
  long long vtime_end = 0;
#endif

  assign("papi-vtime-end",    vtime_end);
  increment("papi-vtime-active", vtime_end);
  increment("papi-fp-ops",  papi_counters_[index_papi_fp_ops]);

  CountersUser::stop_(counters);

}

//----------------------------------------------------------------------

void CountersPapi::update_(std::string key, long long * counters)
{
  TRACE("void CountersPapi::update_");

  bool first_call = (global_[key] == 0);

  CountersUser::update_(key,counters);

  if (first_call) {
    global_[key][0] = vtime_begin_; // Reset vtime_begin
  }
}

//======================================================================

void CountersPapi::papi_create_eventset_()
{
  if (event_set_ == PAPI_NULL) {
    int retval = 0;
    // _CALL_PAPI(PAPI_multiplex_init(),"multiplex_init",PAPI_OK,retval);
    _CALL_PAPI(PAPI_create_eventset(&event_set_), "create_eventset", 
	       PAPI_OK,retval);
    // _CALL_PAPI(PAPI_set_multiplex(event_set_),"set_multiplex",PAPI_OK,retval);
    // if (retval) {
    //   fprintf (stderr, _ERROR  "PAPI_create_eventset() returned 'retval=%d'\n",
    // 	       __FILE__,__LINE__,retval);
    //   fflush(stderr);
    // }
  } else {
    fprintf (stderr,_DEBUG "Trying to re-create existing PAPI eventset!\n",
	     __FILE__,__LINE__);
    fflush(stderr);
  }
}

//----------------------------------------------------------------------

void CountersPapi::papi_delete_eventset_()

{
  if (event_set_ != PAPI_NULL) {
    int retval = 0;
    _CALL_PAPI(PAPI_cleanup_eventset( PAPI3_ARG_ADJUST event_set_),
	       "cleanup_eventset",PAPI_OK,retval);
    if (retval) {
      fprintf (stderr, _ERROR 
	       "PAPI_cleanup_eventset() returned 'retval=%d'\n",
	       __FILE__,__LINE__,retval);
      fflush(stderr);
    }
    _CALL_PAPI(PAPI_destroy_eventset(&event_set_),
	       "destroy_eventset", PAPI_OK,retval);
    if (retval) {
      fprintf (stderr, _ERROR 
	       "PAPI_destroy_eventset() returned 'retval=%d'\n",
	       __FILE__,__LINE__,retval);
      fflush(stderr);
    }
    event_set_ = PAPI_NULL;
  }
}

//----------------------------------------------------------------------

int CountersPapi::papi_insert_counter_(char * counter)

{
  int retval = 0;

  if (! is_papi_active_) {

    // Determine if counter corresponds to a PAPI event

    int code;

    _CALL_PAPI(PAPI_event_name_to_code (counter,&code),
	       "event_name_to_code",PAPI_OK,retval);
    _CALL_PAPI(PAPI_query_event (code),"query_event",PAPI_OK,retval);

    if (retval == PAPI_OK) {

      // Create PAPI eventset if needed

      if (event_set_ == PAPI_NULL) {
	papi_create_eventset_();
      }

      // Insert the event

      _CALL_PAPI(PAPI_add_event (PAPI3_ARG_ADJUST event_set_,code),
		 "add_event",PAPI_OK,retval);
    } else {

      fprintf (stderr, _ERROR 
	       "papi_insert_counter_(%s) called with undefined PAPI event\n",
	       __FILE__,__LINE__,counter);
      fflush(stderr);

    }
  } else {
      fprintf (stderr, _ERROR 
	       "papi_insert_counter_(%s) called with active PAPI eventset!\n",
	       __FILE__,__LINE__,counter);
      fflush(stderr);
  }
  return retval;
}

//----------------------------------------------------------------------

int CountersPapi::papi_delete_counter_(char * counter)

{
  int retval = 0;

  if (! is_papi_active_) {

    // Determine if counter corresponds to a PAPI event

    int code;

    _CALL_PAPI(PAPI_event_name_to_code (counter,&code),
	       "event_name_to_code",PAPI_OK,retval);
    _CALL_PAPI(PAPI_query_event (code),"query_event",PAPI_OK,retval);

    // delete the event

    _CALL_PAPI(PAPI3_remove_event (PAPI3_ARG_ADJUST event_set_,code),
	       "rem_event",PAPI_OK,retval);

  } else {
    fprintf (stderr, _ERROR 
	     "papi_delete_counter_(%s) called with active PAPI eventset!\n",
	     __FILE__,__LINE__,counter);
    fflush(stderr);
  }
  return retval;

}

//----------------------------------------------------------------------

void CountersPapi::papi_start_counters_ ()

{
  if (num_counters_ > 0) {
    if (is_papi_active_) {
      fprintf (stderr, _ERROR 
	       "papi_start_counters_() called with active event set!\n",
	       __FILE__,__LINE__);
      fflush(stderr);
      return;
    }

    int retval = 0;
    _CALL_PAPI(PAPI_start(event_set_),"start",PAPI_OK,retval);
    if (retval) {
      fprintf (stderr, _ERROR "PAPI_start() returned 'retval=%d'\n",
	       __FILE__,__LINE__,retval);
      fflush(stderr);
    }
    is_papi_active_ = true;
  }
  // printf ("multiplexed = %s\n",
  // 	  PAPI_get_multiplex(event_set_) ? "true" : "false");
}

//----------------------------------------------------------------------

void CountersPapi::papi_stop_counters_ ()

{
  if (num_counters_ > 0) {

    if (! is_papi_active_) {
      fprintf (stderr, _ERROR 
	       "papi_stop_counters_() called with inactive event set!\n",
	       __FILE__,__LINE__);
      fflush(stderr);
      return;
    }

    int retval = 0;

    is_papi_active_ = false;

    _CALL_PAPI(PAPI_stop(event_set_, papi_counters_),"stop",PAPI_OK,retval);

    if (retval) {
      fprintf (stderr, _ERROR "PAPI_stop() returned 'retval=%d'\n",
	       __FILE__,__LINE__,retval);
      fflush(stderr);
    }
  }
}

//----------------------------------------------------------------------

void CountersPapi::papi_read_counters_ ()

{
  if (num_counters_ > 0) {

    if (! is_papi_active_) {
      fprintf (stderr, _ERROR 
	       "papi_read_counters_() called with inactive event set!\n",
	       __FILE__,__LINE__);
      fflush(stderr);
      return;
    }

    int retval = 0;

    _CALL_PAPI(PAPI_read(event_set_, papi_counters_),"read",PAPI_OK,retval);

    if (retval) {
      fprintf (stderr, _ERROR "PAPI_read() returned 'retval=%d'\n",
	       __FILE__,__LINE__,retval);
      fflush(stderr);
    }
  }
}

}
#endif
