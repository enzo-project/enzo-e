// See LICENSE_CELLO file for license and copyright information

/// @file     disk_pngio.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-02-06
/// @brief    [\ref Disk] Declaration of functions in the pngio namespace

#ifndef DISK_PNGIO_HPP
#define DISK_PNGIO_HPP

namespace pngio{

  /// writes data stored in values to disk in a png file called fname
  ///
  /// @param[in] fname the output filename
  /// @param[in] values the array of data to write to disk. The values in each
  ///    row are stored contiguously.
  /// @param[in] width, height The dimensions of values (and the output image).
  /// @param[in] colormap Specifies how values are mapped to RGB values
  /// @param[in] vmin,vmax Defines the range of values the colormap covers
  void write(const std::string& fname, double* values, int width, int height,
             const std::vector<float> (&colormap)[3], double vmin, double vmax);

  /// reads data from a png file to create a mask
  ///
  /// @param[in] fname the output filename
  /// @param[out] width, height The dimensions of the mask. The values in each
  ///    each row of the mask are stored contiguously.
  bool* read_as_mask(const std::string& fname, int &width, int &height);

} // namespace pngio

#endif /* DISK_PNGIO_HPP */
