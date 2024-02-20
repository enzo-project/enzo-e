// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/enzo_driver_header.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    \ref Enzo component

// this is the only file in the Enzo-Layer that should include enzo.def.h
//
// enzo.def.h is generated for the enzo Charm++ module defined in enzo.ci
// - ASIDE: technically, the file's name has nothing to do with the name of the
//          .ci file (it's only related to the module name)
// - it's not actually a header file. It includes non-inlined definitions for
//   each generated entity.
// - Consequently, we need to be sure to only include it in one translation
//   unit so we don't duplicate the definitions
//
// Essentially, this file is the one-and-only place where we will include
// enzo.def.h



// first make sure to include the header for every (sub)component of the
// Enzo-layer that that holds a definition referenced inside of enzo.ci
// - in the future, it could make sense to declare separate charm++ modules
//   if we wanted to make the individual subcomponents more independent of
//   each other

#include "Enzo/enzo.hpp"
#include "Enzo/assorted/assorted.hpp"
#include "Enzo/gravity/gravity.hpp"
#include "Enzo/initial/initial.hpp"
#include "Enzo/hydro-mhd/hydro-mhd.hpp"
#include "Enzo/hydro-mhd/toolkit/toolkit.hpp"
#include "Enzo/io/io.hpp"
#include "Enzo/mesh/mesh.hpp"
#include "Enzo/particle/particle.hpp"
#include "Enzo/tests/tests.hpp"
#include "Enzo/utils/utils.hpp"

// next, make sure to include "charm_enzo.hpp" (which includes enzo.decl.h)
#include "Enzo/charm_enzo.hpp"

// finally, move on to include enzo.def.h

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

#include "enzo.def.h"
