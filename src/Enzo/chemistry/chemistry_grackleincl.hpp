// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/chemistry/chemistry_grackleincl.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2024-02-24
/// @brief    Include the <grackle.h> header
///
/// The whole point of this file is to include <grackle.h>.
///
/// Ideally, it would be nice for this to only be included in source files and
/// not in header files.
/// - Doing that would involve forward declaring the key grackle types.
/// - Unfortunately, that's non-trivial since Grackle employs the common C
///   idiom of typedefing an anonymous struct.
/// - This is problematic because you can only forward declare structs (and not
///   typedefs). Since the structs are anonymous, we can't forward declare them
/// - there are some ugly workarounds we can revisit in the future

#ifndef ENZO_CHEMISTRY_CHEMISTRY_PRIVATE_HPP
#define ENZO_CHEMISTRY_CHEMISTRY_PRIVATE_HPP

#ifdef CONFIG_USE_GRACKLE

#include <stdlib.h>
extern "C" {
  #define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
  #include <grackle.h>
}

#else

// declare the names of Grackle types to reduce the usage of ifdef statements
extern "C" { 
  struct chemistry_data;
  struct chemistry_data_storage;
  struct code_units;
  struct grackle_field_data;
}

#endif /* CONFIG_USE_GRACKLE */

#endif /* ENZO_CHEMISTRY_CHEMISTRY_PRIVATE_HPP */
