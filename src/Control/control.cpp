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

  
//   // Initialize domain

//   Domain { extent = [0.0, 0.3, 0.0, 0.3] } 

//   // Initialize grid

//   # The grid is AMR starting with array of 400 cells parallelized using MPI
 
//   Grid {
//     root      = [400,400]; # size of the root grid
//     levels    = 4;         # maximum effective levels (for r_factor = 2)
//     refine    = 4;         # refinement factor
//     type      = tree;      # AMR type: patch or tree
//     full_tree = false;     # whether tree is full (ala Flash) or not
//     backfill  = true;      # whether to backfill for refinement > 2
//     coalesce  = true;      # whether to coalesce small patches to one big one
//     patch_min = 4;         # minimum patch size
//     patch_max = 128;       # maximum patch size
//   }
//   // Define boundary conditions
                                
//   Boundary { type = "reflecting" }

//   # Define stopping criteria

//   // Initialize methods

//   # Only physics method is hydro using ppm

//   Method ppm {
//       diffusion = true;
//       flattening = true;
//       steepening = true;
//   }

//   // Initialize fields


//   # Define field properties with named groups 
                                
//   Field pressure { floor = 1.0e-6 }
//   Field density  { floor = 1.0e-6 }
                                   

//   Eos {
//       gamma     = 1.4;
//   }

//   Timestep {
//       courant   = 0.8;
//   }

//   # Define initial conditions

//   Initial {
//      density  = [1.0, x + y >= 0.15, 0.125];
//      pressure = [1.0, x + y >= 0.15, 0.14];
//      velocity_x = 0;
//      velocity_y = 0;
//   }

//   # Define boundary conditions


//   Stopping {
//      time  = 2.5;
//      cycle = 20000;
//   }

//   Output { 
//      period = 0.5;
//      time = [0.0, 0.3, 0.9, 1.6];
//      interval = [0.3, 0.5, 20]; # output t = 0.3 + k*0.5 for k = [0,20) 
//   }

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

