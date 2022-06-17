// See LICENSE_CELLO file for license and copyright information

/// @file     cello_defines.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs Dec 2 2021
/// @brief    Define macros that wrap language extensions that provide useful
///     optimization hints and are common to almost all mainstream compilers.

#ifndef CELLO_DEFINES_HPP
#define CELLO_DEFINES_HPP

// compilers define the following macros for identification purposes:
// - g++ defines __GNUC__
// - clang defines __clang__ and __GNUC__
// - classic intel compiler defines __INTEL_COMPILER and __GNUC__

/// @def      FORCE_INLINE
/// @brief    For functions/macros that you want to force the compiler to
///           inline, replace the `inline` specifier with `FORCE_INLINE`
#ifdef __GNUC__
// modern C++11 syntax: #define FORCE_INLINE [[gnu::always_inline]] inline
// (but the earliest compiler versions that supports this are unclear)
#define FORCE_INLINE inline __attribute__((always_inline))
#else // don't force inlining for unrecognized compilers
#define FORCE_INLINE inline
#endif

#endif /* CELLO_DEFINES_HPP */
