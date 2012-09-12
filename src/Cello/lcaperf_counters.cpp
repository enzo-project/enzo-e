// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_Counters.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-05-19
/// @brief    Implementation of the Counters abstract base class

#include "lcaperf.hpp"

//----------------------------------------------------------------------

Counters::Counters(int num_counters) throw ()
  : num_counters_(num_counters),
    is_tracing_active_(false),
    is_logging_active_(false),
    global_(),
    index_(),
    name_(),
    type_(),
    frame_(),
    counters_()
{
  TRACE("Counters::Counters");
  name_.resize(num_counters_);
  type_.resize(num_counters_);
}

//----------------------------------------------------------------------

Counters::~Counters() throw ()
{
  // delete global counters
  std::map<std::string,long long *>::iterator iter_global;
  for (iter_global=global_.begin(); iter_global!=global_.end(); ++iter_global) {
    delete [] iter_global->second;
  }
  global_.clear();
  // delete key counters
  while (! counters_.empty()) {
    delete [] counters_.top();
    counters_.pop();
  }
}

//----------------------------------------------------------------------

Counters::Counters(const Counters & counters) throw ()
/// @param     counters  Object being copied
{
  TRACE("Counters::Counters");
}

//----------------------------------------------------------------------

Counters & Counters::operator= (const Counters & counters) throw ()
/// @param     counters  Source object of the assignment
/// @return    The target assigned object
{
  TRACE("Counters & Counters::operator=");
  return *this;
}

//======================================================================

void Counters::begin ()
{
  TRACE("void Counters::begin");
  if (is_tracing_active_) {
    fprintf (stderr, _ERROR "begin() called while in active mode!\n",
	     __FILE__,__LINE__);
    fflush(stderr);
    return;
  }

  // Change tracing state to active

  is_tracing_active_ = true;

}

//----------------------------------------------------------------------

void Counters::end ()
{
  TRACE("void Counters::end");
  if (! is_tracing_active_) {
    fprintf (stderr,_ERROR "end() called after event counting stopped\n",
	     __FILE__,__LINE__);
    fflush(stderr);
  }

  // Change tracing state to inactive

  is_tracing_active_ = false;

  // Output performance summary
  //  print();

  // Clear global counter values
  std::map<std::string,long long *>::iterator iter_global;
  for (iter_global=global_.begin(); iter_global!=global_.end(); ++iter_global) {
    delete [] iter_global->second;
  }
  global_.clear();

}

//----------------------------------------------------------------------

void Counters::start (std::string region, 
		      const Attributes * attributes)
{
  TRACE("void Counters::start");
  // Make sure we are in active state

  if (!is_tracing_active_) {
    fprintf (stderr,_ERROR "start(\"%s\") called before begin()\n",
             __FILE__,__LINE__,region.c_str());
    fflush(stderr);
    return;
  }

  // Generate key given region and current attributes

  std::string key = generate_key_(region,attributes);

  // Save key on frame stack

  frame_.push(key);

  // Create and start a new set of counters

  long long * counters = start_();

  // Push local counters on stack

  counters_.push(counters);

}

//----------------------------------------------------------------------

void Counters::stop (std::string region,
		     const Attributes * attributes)
{
  TRACE("void Counters::stop");

  // Make sure we are in active state

  if (!is_tracing_active_) {
    fprintf (stderr,_ERROR "stop(\"%s\") called before begin()\n",
             __FILE__,__LINE__,region.c_str());
    fflush(stderr);
    return;
  }

  // Pop local counters from stack

  long long * counters = counters_.top();

  counters_.pop();

  // Stop local counters

  stop_ (counters);

  // Generate key given region and current attributes
  
  std::string key = generate_key_(region, attributes);

  // Check that key is expected

  std::string key_start = frame_.top();
  frame_.pop();

  if (key_start != key) {
    fprintf (stderr,_ERROR 
	     "key changed from \"%s\" to \"%s\" in stop(%s)\n",
             __FILE__,__LINE__,key_start.c_str(), key.c_str(),region.c_str());
    fflush(stderr);
    exit(1);
    //    key = merge_keys_(key_start,key);
  }

  // Update global counters for the key, and delete current counters

  update_ (key, counters);
 
}

//----------------------------------------------------------------------

long long Counters::value (std::string key, std::string counter)
{
  TRACE("long long Counters::value");
  if (global_.find(key) == global_.end()) {
    fprintf (stderr,_ERROR 
	     "value(\"%s\", \"%s\") called with unknown key \"%s\"\n",
             __FILE__,__LINE__,key.c_str(), counter.c_str(),key.c_str());
    fflush(stderr);
    return 0;
  }
  if (index_.find(counter) == index_.end()) {
    fprintf (stderr,_ERROR 
	     "value(\"%s\",\"%s\") called with unknown counter \"%s\"\n",
             __FILE__,__LINE__,key.c_str(), counter.c_str(), counter.c_str());
    fflush(stderr);
    return 0;
  }
  
  return global_[key][index_[counter]];

}

//----------------------------------------------------------------------

void Counters::clear ()
{
  TRACE("void Counters::clear");
  bool error = false;
  if (counters_.size() > 0) {
    int size = counters_.size();
    fprintf (stderr, _ERROR 
	     "clear() called with counters_.size() = %d > 0 !\n",
	     __FILE__,__LINE__,size);
    error = true;
  }
  if (frame_.size() > 0) {
    int size = frame_.size();
    fprintf (stderr, _ERROR 
	     "clear() called with frame_.size() = %d > 0 !\n",
	     __FILE__,__LINE__,size);
    error = true;
  }
  if (! error) {
  // Clear global counter values
    std::map<std::string,long long *>::iterator iter_global;
    for (iter_global=global_.begin(); iter_global!=global_.end(); ++iter_global) {
      delete [] iter_global->second;
    }
    global_.clear();
  }
}

//----------------------------------------------------------------------

void Counters::print (FILE * fp)
{
  TRACE("void Counters::print");
  std::map<std::string,long long *>::iterator iter;
  for (iter = global_.begin();
       iter != global_.end();
       ++iter) {
    long long * counters = iter->second;
    if (counters) {
      for (int i=0; i<num_counters_; i++) {
	if (fp) fprintf (fp,"%s %s %lld\n",
			 iter->first.c_str(),name_[i].c_str(),counters[i]);
      }
    }
  }
}

//======================================================================

std::string Counters::generate_key_ 
(
 std::string        region,
 const Attributes * attributes
 ) const
{
  return region + attributes->get_key();
  // std::string key;

  // key = std::string(region);

  // // Set key = region + [":" + attribute-value]

  // for (int i=0; i<attributes->size(); i++) {
  //   key = key + ":" + attributes->value(i);
  // }

  // return key;
}

//----------------------------------------------------------------------

std::string Counters::merge_keys_ (std::string key_1, std::string key_2) const
{
  TRACE("merge_keys IN");
  std::string new_key;

  // Check that region portions match

  if (key_1.substr(0,key_1.find(":")) != key_2.substr(0,key_2.find(":"))) {
    fprintf (stderr, _ERROR 
	     "merge_augmented_regions_() regions %s and %s do not match!\n",
	     __FILE__,__LINE__,
	     (key_1.substr(0,key_1.find(":"))).c_str(),
	     (key_2.substr(0,key_2.find(":"))).c_str());
    fflush(stderr);
  } else {
    new_key = key_1.substr(0,key_1.find(":"));
  }

  // Change differing Attribute values to LCAP_STRING_NULL

  std::string tail1 
    = key_1.substr(key_1.find(":")+1,key_1.size()-key_1.find(":"));
  std::string tail2 
    = key_2.substr(key_2.find(":")+1,key_2.size()-key_2.find(":"));

  while (tail1 != "") {
    if (tail1.substr(0,tail1.find(":")) != tail2.substr(0,tail2.find(":"))) {
      new_key = new_key + ":" + LCAP_STRING_NULL;
    } else {
      new_key = new_key + ":" + tail1.substr(0,tail1.find(":"));
    }

    tail1 = tail1.substr(tail1.find(":")+1,tail1.size()-tail1.find(":"));
    tail2 = tail2.substr(tail2.find(":")+1,tail2.size()-tail2.find(":"));
  }

  TRACE("merge_keys OUT");
  return new_key;
  
}

//======================================================================

// ItCounterKeys::ItCounterKeys ( Counters * counters ) throw ()
//   : counters_(counters),
//     iter_(counters->global_.begin()),
//     value_(0)
// {}

// //----------------------------------------------------------------------

// ItCounterKeys::~ItCounterKeys ( ) throw ()
// {}

// //----------------------------------------------------------------------

// const char * ItCounterKeys::operator++ () throw()
// {
//   const char * key;
//   if (iter_ == counters_->global_.end()) {
//     key = 0;
//     value_ = 0;
//     iter_ = counters_->global_.begin();
//   } else {
//     key    = iter_->first.c_str();
//     value_ = iter_->second;
//     ++iter_;
//   }
//   return key;
// }

// //----------------------------------------------------------------------

// bool ItCounterKeys::done () const throw()
// {
//   return iter_ == counters_->global_.begin();
// }

// //----------------------------------------------------------------------

// long long * ItCounterKeys::value () const throw()
// {
//   return value_;
// }




