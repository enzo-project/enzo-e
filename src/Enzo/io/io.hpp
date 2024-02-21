// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/io/io.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for io subcomponent within the \ref Enzo layer

#ifndef ENZO_IO_IO_HPP
#define ENZO_IO_IO_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/mesh.hpp"    // Block
#include "Cello/io.hpp"      // Io
#include "Cello/problem.hpp" // Method, Initial

#include "Enzo/enzo.hpp" // enzo_float, EnzoBlock

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

#include "io/IoEnzoBlock.hpp"
#include "io/IoEnzoReader.hpp"
#include "io/IoEnzoWriter.hpp"

#include "io/EnzoMethodCheck.hpp"

#include "io/EnzoInitialHdf5.hpp"
#include "io/EnzoInitialMusic.hpp"

#endif /* ENZO_IO_IO_HPP */
