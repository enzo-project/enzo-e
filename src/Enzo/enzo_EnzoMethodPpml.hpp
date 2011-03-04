// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpml.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon May 17 14:16:01 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo PPML method

#ifndef ENZO_ENZO_METHOD_PPML_HPP
#define ENZO_ENZO_METHOD_PPML_HPP

class EnzoMethodPpml : public Hyperbolic {

/// @class    EnzoMethodPpml
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulate Enzo's PPML hydro method

public: // interface

/// Creae a new EnzoMethodPpml object
  EnzoMethodPpml(Parameters * parameters,
		 EnzoBlock * enzo)
    : Hyperbolic (parameters),
      enzo_(enzo)
  {};

  void compute_block( DataBlock * data_block,
		      double t, double dt ) throw(); 

private: // attributes

  EnzoBlock * enzo_;

};

#endif /* ENZO_ENZO_METHOD_PPML_HPP */
