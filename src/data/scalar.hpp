//345678901234567890123456789012345678901234567890123456789012345678901234567890

#ifndef SCALAR_HPP
#define SCALAR_HPP

/// Defines Scalar

// $Id$
/**
 * @file    scalar.hpp
 * @brief   Header file for the Array class
 * @author  James Bordner 
 * @version 1.0
 *
 * Definitions
 *
 * (*) Scalar
 *
 */
// $Log$

/// Define Scalar to be double.  Should be more flexible and easy to change.

#define Scalar double
#define SCALAR_SCANF "%lf"
#define SCALAR_PRINTF "%le "
#define MPI_SCALAR MPI_DOUBLE

#endif
