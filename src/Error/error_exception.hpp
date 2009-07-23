#ifndef ERROR_EXCEPTION_HPP
#define ERROR_EXCEPTION_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      error_exception.hpp
 * @brief     Exception class hierarchy
 * @author    James Bordner
 * @date      Sun Jul 12 13:12:09 PDT 2009
 *
 * Exception class hierarchy
 *
 * $Id$
 *
 *********************************************************************
 */

class Exception {

/** 
 *********************************************************************
 *
 * @class     Exception
 * @brief     Exception bass class
 * @ingroup   Error
 *
 * Exception bass class
 *
 *********************************************************************
 */

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------


};

class ExceptionBadPointer : public Exception {};

#endif /* EXCEPTION_HPP */

