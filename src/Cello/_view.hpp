// See LICENSE_CELLO file for license and copyright information

/// @file     _view.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon Jun 3 2019
/// @brief    Private include file for the \ref View component

#ifndef _VIEW_HPP
#define _VIEW_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdio.h>
#include <cstddef>
#include <type_traits>
#include <limits>
#include <memory>

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "view_CelloView.hpp"

// for backwards compatability, we define CelloArray as an alias for CelloView
// (new code should avoid using this alias)
template<typename T, std::size_t D>
using CelloArray = CelloView<T, D>;

#include "view_ViewCollec.hpp"
#include "view_StringIndRdOnlyMap.hpp"

#endif /* _VIEW_HPP */
