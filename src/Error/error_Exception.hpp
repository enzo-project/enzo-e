// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     error_Exception.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Sun Jul 12 13:12:09 PDT 2009
/// @brief    Declaration of Exception classes

#ifndef ERROR_EXCEPTION_HPP
#define ERROR_EXCEPTION_HPP

/// @class    Exception
/// @ingroup  Error
/// @brief    Base class for the Exception class hierarchy
class Exception { };

/// @class    ExceptionBadPointer
/// @ingroup  Error
/// @brief    Global exception class for encountering a bad pointer
class ExceptionBadPointer           : public Exception {};

/// @class    ExceptionBadArrayDeallocation
/// @ingroup  Error
/// @brief    Array exception class for bad deallocation
class ExceptionBadArrayDeallocation : public Exception {};

/// @class    ExceptionBadArrayAllocation
/// @ingroup  Error
/// @brief    Array exception class for bad allocation
class ExceptionBadArrayAllocation   : public Exception {};

/// @class    ExceptionMemoryBadDeallocate
/// @ingroup  Error
/// @brief    Memory exception class for bad deallocation
class ExceptionMemoryBadDeallocate  : public Exception {};

/// @class    ExceptionMemoryBadAllocate
/// @ingroup  Error
/// @brief    Memory exception class for bad allocation
class ExceptionMemoryBadAllocate    : public Exception {};

/// @class    ExceptionMemoryBadGroupHandle
/// @ingroup  Error
/// @brief    Memory exception class for a bad group handle
class ExceptionMemoryBadGroupHandle : public Exception {};

/// @class    ExceptionParametersBadType
/// @ingroup  Error
/// @brief    Parameters exception class for bad parameter type
class ExceptionParametersBadType    : public Exception {};

#endif /* EXCEPTION_HPP */

