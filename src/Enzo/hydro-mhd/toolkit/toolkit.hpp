// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/hydro-mhd/toolkit/toolkit.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-06-01
/// @brief    Include file for generic tools used to implement a hydro-mhd
///           solver
///
/// This is currently called toolkit because it provides generic functionality
/// that could be used for implementing hydro-mhd solvers. However, at the time
/// of writing this blurb, its currently just used by the VL+CT solver. For
/// that reason, it could make more sense to call this vlct_impl (Or after GH
/// PR #315 gets merged, multistage_impl).

#ifndef ENZO_HYDROMHD_TOOLKIT_TOOLKIT_HPP
#define ENZO_HYDROMHD_TOOLKIT_TOOLKIT_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <array>
#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include <charm++.h> // PUP::er

#include "Cello/error.hpp"
#include "Cello/view.hpp" // CelloView
#include "Enzo/enzo.hpp" // enzo_float, str_vec_t, EFltArrayMap

//----------------------------------------------------------------------
// Component Headers
//----------------------------------------------------------------------

#include "Enzo/hydro-mhd/toolkit/EnzoIntegrationQuanUpdate.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoLazyPassiveScalarFieldList.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoSourceGravity.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoSourceInternalEnergy.hpp"

// [order dependencies:]
#include "Enzo/hydro-mhd/toolkit/EnzoBfieldMethod.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoBfieldMethodCT.hpp"

// [order dependencies:]
#include "Enzo/hydro-mhd/toolkit/EnzoReconstructor.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoReconstructorNN.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoReconstructorPLM.hpp"

#endif /* ENZO_HYDROMHD_TOOLKIT_TOOLKIT_HPP */
