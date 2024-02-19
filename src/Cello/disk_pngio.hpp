// See LICENSE_CELLO file for license and copyright information

/// @file     disk_pngio.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-02-06
/// @brief    [\ref Disk] Declaration of functions in the pngio namespace

#ifndef DISK_PNGIO_HPP
#define DISK_PNGIO_HPP

namespace pngio{

  enum class ImgTransform{none, abs, log};

  /// writes data stored in values to disk in a png file called fname
  ///
  /// @param[in] fname the output filename
  /// @param[in] data the array of data to write to disk. The elements in each
  ///    row are assumed to be stored contiguously in memory.
  /// @param[in] width, height The dimensions of data (and the output image).
  /// @param[in] colormap Specifies how to map elements of data to RGB values
  /// @param[in] transform Specifies the transformation (if any) applied to
  ///    the elements of data before they are mapped to a color
  /// @param[in] min_max An optional array of 2 elements. The first element
  ///    and second element respectively specify the minimum and maximum of
  ///    the range of values that the colormap covers (after the transform has
  ///    been applied to entries of data). Values outside of this range are
  ///    clamped to the nearest value. Pass ``nullptr`` to indicate that the
  ///    full range of elements in data should be mapped to colors.
  void write(const std::string& fname, double* values, int width, int height,
             const std::vector<float> (&colormap)[3],
             ImgTransform transform,
             const std::array<double,2>* min_max) noexcept;

  /// reads data from a png file to create a mask. The mask is true if any RGB
  /// channel in the image is non-zero
  ///
  /// @param[in] fname the output filename
  /// @param[out] width, height The dimensions of the mask. The values in each
  ///    each row of the mask are stored contiguously.
  bool* read_as_mask(const std::string& fname,
                     int &width, int &height) noexcept;

} // namespace pngio

#endif /* DISK_PNGIO_HPP */
