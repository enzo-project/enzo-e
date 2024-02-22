// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructorNN.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Nearest Neighbor
///           Reconstruction

#ifndef ENZO_ENZO_RECONSTRUCTOR_NN_HPP
#define ENZO_ENZO_RECONSTRUCTOR_NN_HPP

class EnzoReconstructorNN : public EnzoReconstructor
{
  /// @class    EnzoReconstructorNN
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates nearest neighbor reconstruction of
  ///           primitives at cell interfaces

public: // interface

  /// Create a new EnzoReconstructorNN
  EnzoReconstructorNN(std::vector<std::string> active_key_names)
    : EnzoReconstructor(active_key_names)
  { }

  void reconstruct_interface
  (const EnzoEFltArrayMap &prim_map, EnzoEFltArrayMap &priml_map,
   EnzoEFltArrayMap &primr_map, const int dim, const int stale_depth,
   const str_vec_t& passive_list);

  int total_staling_rate()
  { return 1; }

  int immediate_staling_rate()
  { return 0; }

};

#endif /* ENZO_ENZO_RECONSTRUCTOR_NN_HPP */
