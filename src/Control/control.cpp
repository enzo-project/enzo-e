//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      control.cpp
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

#include "control.hpp"
 
//====================================================================
// PUBLIC FUNCTIONS
//====================================================================

Control::Control(Monitor * monitor)
  : monitor_(monitor),
    parameters_(NULL),
    simulation_(NULL)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a Control object
 *
 *********************************************************************
 */
{
}

Control::~Control()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Delete a Control object
 *
 *********************************************************************
 */
{
}

/// 
void Control::read_parameters(FILE * fp)
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
    WARNING_MESSAGE("Control::read_parameters","Parameters object already exists");
    delete parameters_;
    parameters_ = 0;
  }

  parameters_ = new Parameters;

  parameters_->read(fp);
}

void Control::initialize_simulation()
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
  INCOMPLETE_MESSAGE("Control::initialize_simulation","");

}

void Control::execute_simulation()
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
  INCOMPLETE_MESSAGE("Control::execute_simulation","");
}

void Control::terminate_simulation()
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
  INCOMPLETE_MESSAGE("Control::terminate_simulation","");
}

//====================================================================
// PRIVATE FUNCTIONS
//====================================================================

