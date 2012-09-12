// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_Lcaperf.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2004-03-12
/// @brief    Parallel performance monitoring class

#include "lcaperf.hpp"

//----------------------------------------------------------------------

LcaPerf::LcaPerf ()
  : attributes_(),
    counters_()
{
  TRACE("lcaPerf");
  counters_["basic"]  = new CountersBasic;
  TRACE1("counters = %p",counters_["basic"]);
  counters_["user"]   = new CountersUser;
#ifdef CONFIG_TRACE_MEM
  counters_["mem"] = new CountersMem;
#endif
#ifdef CONFIG_USE_PAPI
  counters_["papi"]   = new CountersPapi;
#endif
#ifdef CONFIG_USE_MPI
  counters_["mpi"]    = new CountersMpi;
#endif

}

//----------------------------------------------------------------------

LcaPerf::~LcaPerf () // + Finalize the LcaPerf object
{
  TRACE("~LcaPerf");
  std::map<std::string,Counters *>::iterator iter;
  
  for (iter = counters_.begin(); iter != counters_.end(); ++iter) {
    delete iter->second;
  }
}

//----------------------------------------------------------------------

void LcaPerf::new_region  (const char * region)
{
  TRACE("new_region");
  // NOTE: does not check whether region is already in regions_ vector
  regions_.push_back(region);
}

//----------------------------------------------------------------------

void LcaPerf::new_attribute  (const char * name, int type = 0)
{
  TRACE("new_attribute");
  attributes_.create(name);
}

//----------------------------------------------------------------------

void LcaPerf::new_counter (const char * name, 
			   counters_type type)
{
  TRACE("new_counter");
  CountersUser * counters_user = 
    dynamic_cast<CountersUser *> (counters_["user"]);

  if (counters_user) counters_user->create(name,type);
}

//----------------------------------------------------------------------

long long LcaPerf::value (const char * set, 
			  const char * key, 
			  const char * counter)
{
  TRACE("value");
  return counters_[set]->value (key,counter);
}

//----------------------------------------------------------------------

void LcaPerf::delete_attribute (const char * attribute_name)
{
  TRACE("delete_attribute");
  attributes_.remove(attribute_name);
}

//----------------------------------------------------------------------

void LcaPerf::delete_counter (const char * name)
{
  TRACE("delete_counter");
  CountersUser * counters_user = 
    dynamic_cast<CountersUser *> (counters_["user"]);

  if (counters_user) counters_user->remove(name);
}

//----------------------------------------------------------------------

void LcaPerf::initialize (const char * filename)
{
  TRACE("initialize");
  // Print header (only on root process)

  int ip = 0;

#ifdef CONFIG_USE_MPI
  int init_called = 0;
  MPI_Initialized (&init_called);
  if (init_called) {
    MPI_Comm_rank(MPI_COMM_WORLD,&ip);
  }
#endif

  // Display header
  if (ip == 0) {
    printf ("---------------------------\n");
    printf ("lcaperf: Version %d.%d\n",LCAP_VERSION_MAJOR,LCAP_VERSION_MINOR);
    printf ("lcaperf: Copyright 2011, James Bordner and the Regents of the\n");
    printf ("lcaperf: University of California.  All rights reserved.\n");
    printf ("---------------------------\n");
  }

  std::map<std::string,Counters *>::iterator iter;
  
  for (iter = counters_.begin(); iter != counters_.end(); ++iter) {
    iter->second->initialize();
  }
}

//----------------------------------------------------------------------

void LcaPerf::finalize ()
{
  TRACE("finalize");
}

//----------------------------------------------------------------------

void LcaPerf::begin ()
{
  TRACE("begin");
// #ifdef CONFIG_USE_MPI
//   // Define mpi-size and mpi-rank attributes
//   int ip,np;
//   MPI_Comm_size (MPI_COMM_WORLD,&np);
//   MPI_Comm_rank (MPI_COMM_WORLD,&ip);
//   new_attribute("mpi-size");
//   new_attribute("mpi-rank");
//   attribute("mpi-size",&np, LCAP_INT);
//   attribute("mpi-rank",&ip, LCAP_INT);
// #endif

  std::map<std::string,Counters *>::iterator iter;
  
  for (iter = counters_.begin(); iter != counters_.end(); ++iter) {
    iter->second->begin();
  }
}

//----------------------------------------------------------------------

void LcaPerf::end ()
{
  TRACE("end");
  attributes_.remove("*");

  std::map<std::string,Counters *>::iterator iter;
  
  for (iter = counters_.begin(); iter != counters_.end(); ++iter) {
    iter->second->end();
  }
}

//----------------------------------------------------------------------

void LcaPerf::start (const char * region_base)
{
  TRACE1("start(%s)",region_base);
  std::map<std::string,Counters *>::iterator iter;
  
  for (iter = counters_.begin(); iter != counters_.end(); ++iter) {
    iter->second->start(region_base,&attributes_);
  }
}

//----------------------------------------------------------------------

void LcaPerf::stop (const char * region_base)
{
  TRACE1("stop(%s)",region_base);
  std::map<std::string,Counters *>::iterator iter;
  
  for (iter = counters_.begin(); iter != counters_.end(); ++iter) {
    iter->second->stop(region_base,&attributes_);
  }
}

//----------------------------------------------------------------------

void LcaPerf::increment (const char *name, long long value)
{
  TRACE("increment");
  CountersUser * counters_user = 
    dynamic_cast<CountersUser *> (counters_["user"]);

  if (counters_user) counters_user->increment(name,value);
}

//----------------------------------------------------------------------

void LcaPerf::assign (const char *name, long long value)
{
  TRACE("assign");
  CountersUser * counters_user = 
    dynamic_cast<CountersUser *> (counters_["user"]);

  if (counters_user) counters_user->assign(name,value);
}

//----------------------------------------------------------------------

void LcaPerf::attribute (const char * name, 
			 const void * value, 
			 int          type)
{
  TRACE("attribute");
  char attribute_string[80];

  if (value == 0) {

    strcpy (attribute_string,"*");

  } else {

    switch (type) {
    case LCAP_STRING: 
      sprintf (attribute_string,"%s",(char *) value);
      break;
    case LCAP_INT:
      sprintf (attribute_string,"%d",*((int *)value));
      break;
    case LCAP_DOUBLE:
      sprintf (attribute_string,"%lf",*((double *)value));
      break;
    default:
      fprintf (stderr,_ERROR "Unknown type %d in attribute(\"%s\")\n",
	       __FILE__,__LINE__,type,name);
      fflush(stderr);
      break;
    }

  }
  attributes_.assign(name, attribute_string);
}

#ifdef CONFIG_USE_MPI
const bool l_mpi = true;
#else
const bool l_mpi = false;
#endif

#ifdef CONFIG_USE_PAPI
const bool l_papi = true;
#else
const bool l_papi = false;
#endif

#ifdef USE_MEM
const bool l_mem = true;
#else
const bool l_mem = false;
#endif

#ifdef USE_HDF5
const bool l_disk = true;
#else
const bool l_disk = false;
#endif

//----------------------------------------------------------------------

void LcaPerf::header ()
{
    printf ("lcaperf: ");
    printf             ("   time(s)   " "   ");
    if (l_mpi)  printf ("   MPI(s)    " "   ");
    if (l_papi) printf (" flops(GF)   " "   ");
    if (l_mem)  printf (" memory(GB)  " "   ");
    if (l_disk) printf ("  disk(GB)   " "   ");
    //    printf ("cycle ");
    //    printf ("region\n");
}

//----------------------------------------------------------------------

void LcaPerf::print ()
{
  DEBUG("LcaPerf::print");
  TRACE("print");

  // Get comm size, required for computing average times etc. over all
  // processes

  int np = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD,&np);
#endif

  int         cycle_index  = attributes_.index("cycle");
  std::string cycle_value  = attributes_.value(cycle_index);
  int         level_index  = attributes_.index("level");
  std::string level_value  = attributes_.value(level_index);

  const int i_avg = 0, i_max = 1;
  double time_avg, time_max, time_eff;

  long long counter_array_reduce[2];

  // Line to print
  std::string line;
  char field[80];

  // Loop over regions to print

  for (size_t i_region = 0; i_region < regions_.size(); ++i_region) {

    bool empty = true;
    
    line = "lcaperf: ";

    std::string region = regions_[i_region];

    // Create the region key to check counter keys against

    //    NOTE: EvolveLevel must be handled differently since it is recursive

    std::string level_string = "*";

    //    WARNING: dependency on number of attributes and attribute ordering


    std::stringstream convert;
    convert << region << ":" << cycle_value << ":" << level_value;
    std::string region_key = convert.str();
    DEBUG1 ("region = %s",region_key.c_str());

  
    //--------------------------------------------------
    // TOTAL TIME
    //--------------------------------------------------

    // Loop over keys in the Counters object to sum counters over levels

    counter_array_reduce[i_avg] = 0;
    counter_array_reduce[i_max] = 0;

    ItCounterKeys itKeysBasic (counters_["basic"]);
    int i_time = counters_["basic"]->index("time");

    DEBUG0;

    while (const char * key = ++itKeysBasic) {

      DEBUG0;

      // Select matching keys

      DEBUG2 ("key = %s  region_key = %s",key,region_key.c_str());
      bool keys_match = attributes_.keys_match(key,region_key);

      if (keys_match) {
 	// Get array of counters
 	long long * counter_array = itKeysBasic.value();
 	// Sum over levels
 	counter_array_reduce[i_avg] += counter_array[i_time];
 	counter_array_reduce[i_max] += counter_array[i_time];
      }
    }

#ifdef USE_MPI
    MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_avg],1,MPI_LONG_LONG,
 		   MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_max],1,MPI_LONG_LONG,
 		   MPI_MAX,MPI_COMM_WORLD);
#endif

    time_avg = 0.0;
    time_eff = 1.0;

    if (counter_array_reduce[i_max] != 0) {

      empty = false;

      time_avg = 1e-6*counter_array_reduce[i_avg]/np;
      time_max = 1e-6*counter_array_reduce[i_max];
      time_eff = time_avg / time_max;
    }

    sprintf (field, "%06.2f %5.4f   ",time_avg,time_eff);

    line = line + field;

    //--------------------------------------------------
    // MPI TIME
    //--------------------------------------------------

    if (counters_.find("mpi") != counters_.end()) {

      counter_array_reduce[i_avg] = 0;
      counter_array_reduce[i_max] = 0;

      ItCounterKeys itKeysMpi (counters_["mpi"]);
      i_time = counters_["mpi"]->index("mpi-time");

	DEBUG0;
      while (const char * key = ++itKeysMpi) {
	DEBUG0;

	// Select matching keys

	bool keys_match = attributes_.keys_match(key,region_key);

	if (keys_match) {
	  // Get array of counters
	  long long * counter_array = itKeysMpi.value();
	  // Sum over levels
	  counter_array_reduce[i_avg] += counter_array[i_time];
	  counter_array_reduce[i_max] += counter_array[i_time];
	}
      }

#ifdef USE_MPI
      MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_avg],1,MPI_LONG_LONG,
		     MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE,&counter_array_reduce[i_max],1,MPI_LONG_LONG,
		     MPI_MAX,MPI_COMM_WORLD);
#endif

      time_avg = 0.0;
      time_eff = 1.0;

      if (counter_array_reduce[i_max] != 0) {

	empty = false;

	time_avg = 1e-6*counter_array_reduce[i_avg]/np;
	time_max = 1e-6*counter_array_reduce[i_max];
	time_eff = time_max ? (time_avg / time_max) : 1.0;

      }

      sprintf (field, "%06.2f %5.4f   ",time_avg,time_eff);

      line = line + field;
    }

    if (! empty) {

      sprintf (field , "%s   %s\n",cycle_value.c_str(),region.c_str());

      line = line + field;

      printf ("%s",line.c_str());
    }
  }

  printf ("\n");

  if (counters_.find("basic") != counters_.end()) counters_["basic"]->clear();
  if (counters_.find("mpi") != counters_.end())   counters_["mpi"]->clear();

}
