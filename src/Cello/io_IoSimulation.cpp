// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoSimulation class

#include "io.hpp"

//----------------------------------------------------------------------

IoSimulation::IoSimulation(const Simulation * s) throw ()
  : Io(),
    rank_ (s->rank_), 
    cycle_(s->cycle_),
    time_ (s->time_),
    dt_   (s->dt_)
{
  meta_name_.push_back("rank");
  meta_name_.push_back("cycle");
  meta_name_.push_back("time");
  meta_name_.push_back("dt");
}


//----------------------------------------------------------------------

void IoSimulation::meta_value
(int index,
 void ** buffer, std::string * name, int * type,
 int * mx, int * my, int * mz) throw()
{
  Io::meta_value(index,buffer,name,type,mx,my,mz);

  int count = 0;

  if (index == count++) {

    (*buffer) = (void *) &rank_;
    (*type) = type_int;
    (*mx)  = 1;

  } else if (index == count++) {

    (*buffer) = (void *) &cycle_;
    (*type) = type_int;
    (*mx)  = 1;

  } else if (index == count++) {

    (*buffer) = (void *) &time_;
    (*type) = type_double;
    (*mx)  = 1;

  } else if (index == count++) {

    (*buffer) = (void *) &dt_;
    (*type) = type_double;
    (*mx)  = 1;

  }
}

//----------------------------------------------------------------------

void IoSimulation::save_to (void * v)
{
  Simulation * s = static_cast<Simulation *>(v);

  s->rank_  = rank_;
  s->cycle_ = cycle_;
  s->time_  = time_;
  s->dt_    = dt_;;
}


