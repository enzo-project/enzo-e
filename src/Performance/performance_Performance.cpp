// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      performance_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2009-10-15
/// @brief     Implementation of the Performance class

#include "cello.hpp"

#include "performance.hpp"

Performance::Performance 
(
 unsigned num_attributes,
 unsigned num_counters,
 unsigned num_groups,
 unsigned num_regions)
  : counters_(),
    num_attributes_            (num_attributes),
    attribute_names_           (NULL),
    attribute_types_           (NULL),
    num_counters_              (num_counters),
    counter_names_             (NULL),
    num_groups_                (num_groups),
    current_group_             (0),
    group_names_               (NULL),
    region_names_              (NULL),
    current_region_            (0)
{
  attribute_names_.resize(num_attributes + 1);
  attribute_types_.resize(num_attributes + 1);

  counter_names_.resize(num_counters + 1);

  group_names_.resize(num_groups + 1);

  region_names_.resize(num_regions + 1);

  // Create initial Counters object

  counters_.push_back(new Counters(num_attributes,num_counters));

}

//----------------------------------------------------------------------

Performance::~Performance()
{

  for (unsigned i=0; i<counters_.size(); i++) {
    delete counters_.at(i);
  }

}

//----------------------------------------------------------------------

void Performance::start () throw ()
{
  timer.start();
  papi.start();
}

//----------------------------------------------------------------------

void Performance::stop () throw ()
{
  timer.stop();
  papi.stop();
}

//----------------------------------------------------------------------

void Performance::print () const throw ()
{
  timer.print();
  papi.print();
}

//----------------------------------------------------------------------

Performance::Performance(const Performance & classname) throw()
{
  INCOMPLETE("Performance::Performance","");
}

//----------------------------------------------------------------------

Performance & Performance::operator= (const Performance & classname) throw()
{
  INCOMPLETE("Performance::operator =","");
  return *this;
}

//----------------------------------------------------------------------

void Performance::new_attribute(unsigned            id_attribute, 
				std::string         attribute_name,
				attribute_type_enum type)
/// @param    id_attribute
/// @param    attribute_name
/// @param    type
{
  new_item_ (attribute_names_, 
	     id_attribute, 
	     attribute_name);

  attribute_types_[id_attribute] = type;
}

//----------------------------------------------------------------------

int Performance::attribute(unsigned id_attribute)
{
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_attribute(unsigned id_attribute)
{
}

//----------------------------------------------------------------------

void Performance::new_group(unsigned    id_group, 
			    std::string group_name)
{
  new_item_ (group_names_,
	     id_group, 
	     group_name);
}

//----------------------------------------------------------------------

int Performance::group(unsigned id_group)
{
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_group(unsigned id_group)
{
}

//----------------------------------------------------------------------

void Performance::begin_group(unsigned group_id)
{
  
  if ( current_group_ ){
    char message [ ERROR_LENGTH ];
    sprintf (message, "Performance group started when one already active");
    WARNING("Performance::begin_group",message);
  }

  current_group_ = group_id;

}

//----------------------------------------------------------------------

void Performance::end_group(unsigned id_group)
{
  if (id_group != current_group_) {
    char message [ ERROR_LENGTH ];
    sprintf (message, "Mismatch between begin_group(%s) and end_group(%s)",
	     group_names_[current_group_].c_str(),
	     group_names_[id_group].c_str());
    WARNING("Performance::end_group",message);
  }

  current_group_ = 0;

}

//----------------------------------------------------------------------

void Performance::new_region(unsigned    id_region, 
			     std::string region_name)
{
  new_item_ (region_names_,
	     id_region, 
	     region_name);
}

//----------------------------------------------------------------------

int Performance::region(unsigned id_region)
{
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_region(unsigned id_region)
{
}

//----------------------------------------------------------------------

void Performance::start_region(unsigned region_id)
{
}

//----------------------------------------------------------------------

void Performance::stop_region(unsigned region_id)
{
}

//----------------------------------------------------------------------

void Performance::new_counter(unsigned    id_counter,
			      std::string counter_name)
{
  new_item_ (counter_names_, 
	     id_counter, 
	     counter_name);
}

//----------------------------------------------------------------------

type_counter Performance::counter(unsigned id_counter)
{
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_counter(unsigned          id_counter,
			      type_counter value)
{
}

//----------------------------------------------------------------------

void Performance::increment_counter(unsigned          id_counter,
				    type_counter value)
{
}

//----------------------------------------------------------------------

void Performance::flush()
{
}


//----------------------------------------------------------------------

void Performance::new_item_ 
(
 std::vector<std::string> & item_names,
 unsigned       id_item, 
 std::string    item_name
)
{
  item_names[id_item] = item_name;
}

