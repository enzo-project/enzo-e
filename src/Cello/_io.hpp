// See LICENSE_CELLO file for license and copyright information

/// @file     _io.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-11-01
/// @brief    Private include file for the \ref Io component 

#ifndef _IO_HPP
#define _IO_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <limits>
#include <boost/filesystem.hpp>
#include "pngwriter.h"

//----------------------------------------------------------------------
// Typedefs
//----------------------------------------------------------------------

enum meta_type {
  meta_type_file,
  meta_type_group
};

//----------------------------------------------------------------------
// Global functions
//----------------------------------------------------------------------

/// Namespace for global constants and functions
namespace cello {
  namespace io {

    /// Writes the version string to an output file
    inline void write_version_metadata(File* file)
    {
      ASSERT("cello::write_version_metadata()", "File object must not be null",
             file != nullptr);

      // define a variable to ensure CELLO_VERSION is actually defined and that
      // it is replaced with a string
      const char* version_str = CELLO_VERSION;

      // length intentionally excludes the terminating null character
      // (note: this is the default behvaior of strlen).
      int length = static_cast<int>(std::strlen(version_str));

      ASSERT("cello::write_version_metadata()", // sanity check!
             "version string must have at least 1 character", length > 0);
      file->file_write_meta(version_str, "version", type_char, length);
    }
  }
}

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "io_Colormap.hpp"
#include "io_ColormapRGB.hpp"

#include "io_Io.hpp"

#include "io_IoSimulation.hpp"
#include "io_IoBlock.hpp"

#include "io_IoFieldData.hpp"
#include "io_IoHierarchy.hpp"
#include "io_IoParticleData.hpp"
#include "io_IoReader.hpp"
#include "io_IoWriter.hpp"

#include "io_Input.hpp"

#include "io_Output.hpp"
#include "io_OutputCheckpoint.hpp"
#include "io_OutputData.hpp"
#include "io_OutputImage.hpp"

#include "io_Schedule.hpp"
#include "io_ScheduleInterval.hpp"
#include "io_ScheduleList.hpp"


#endif /* _IO_HPP */

