/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

/** 
 *********************************************************************
 *
 * @file      schedule.cpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @bug       
 * @note      
 *
 * DESCRIPTION 
 * 
 *    
 *
 * PACKAGES
 *
 *    
 * 
 * INCLUDES
 *  
 *    
 *
 * PUBLIC FUNCTIONS
 *  
 *    
 *
 * PRIVATE FUCTIONS
 *  
 *    
 *
 * $Id$
 *
 *********************************************************************
 */

#include "schedule.hpp"
 
//====================================================================
// PUBLIC FUNCTIONS
//====================================================================

Schedule::Schedule(Monitor * monitor)
  : monitor_(monitor),
    parameters_(NULL),
    simulation_(NULL)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a Schedule object
 *
 *********************************************************************
 */
{
}

Schedule::~Schedule()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Delete a Schedule object
 *
 *********************************************************************
 */
{
}

/// 
void Schedule::read_parameters(FILE * fp)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Read parameter file
 *
 *********************************************************************
 */
{
  monitor_->print ("Reading parameters");

  if (parameters_ != NULL) {
    WARNING_MESSAGE("Schedule::read_parameters","Parameters object already exists");
    delete parameters_;
    parameters_ = 0;
  }

  parameters_ = new Parameters;

  parameters_->read(fp);
}

void Schedule::initialize_simulation()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Initialize the simulation
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Schedule::initialize_simulation","");

}

void Schedule::execute_simulation()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Run the simulation
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Schedule::execute_simulation","");
}

void Schedule::terminate_simulation()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Terminate the simulation
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Schedule::terminate_simulation","");
}

//====================================================================
// PRIVATE FUNCTIONS
//====================================================================

