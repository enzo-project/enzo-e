// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Array.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-10-27
/// @brief    Implementation of the Adapt class for remeshing

#include "mesh.hpp"

//----------------------------------------------------------------------

void Adapt::allocate_level_bounds (int num_neighbors)
{
  level_curr_.resize(num_neighbors + 1);
  level_want_.resize(num_neighbors + 1);
  level_min_.resize(num_neighbors + 1);
  level_max_.resize(num_neighbors + 1);
  is_sibling_.resize(num_neighbors + 1);
  can_coarsen_.resize(num_neighbors + 1);
}

//----------------------------------------------------------------------

void Adapt::set_initial_level_bounds
(int i, Index index, int level_curr, int level_want, bool is_sibling)
{
  index_map_[index] = i;
  level_curr_[i] = level_curr;
  level_want_[i] = level_want;
  level_min_[i] = level_want;
  level_max_[i] = level_curr + 1;
  is_sibling_[i] = is_sibling;
  can_coarsen_[i] = false;
}

//----------------------------------------------------------------------

void Adapt::update_level_bounds
(Index index, int level_min, int level_max, bool can_coarsen)
{
  const int i = index_map_[index];
  level_min_[i] = std::max(level_min_[i],level_min);
  level_max_[i] = std::min(level_max_[i],level_max);
  can_coarsen_[i] = can_coarsen_[i] || can_coarsen;
}

//----------------------------------------------------------------------

bool Adapt::evaluate_level_bounds ()
{
  // Save values to test later if changed
  int level_min = level_min_[0];
  int level_max = level_max_[0];
  int can_coarsen = can_coarsen_[0];
  
  const int n = level_min_.size();
  for (int i=1; i<n; i++) {
    level_min_[0] = std::max(level_min_[0], (level_min_[i] - 1) );
  }
  int neighbor_max = 0;
  for (int i=1; i<n; i++) {
    neighbor_max = std::max(neighbor_max,level_max_[i]);
  }
  level_max_[0] = std::max(level_min_[0],neighbor_max-1);

  // adjust for coarsening: can only coarsen if all siblings can coarsen

  const bool want_to_coarsen = (level_min_[0] < level_curr_[0]);
  const bool converged = (level_min_[0] == level_max_[0]);
   if ( want_to_coarsen && converged ) {
    // block can coarsen, check that all neighbors can as well
    bool cant_coarsen = false;
    can_coarsen_[0] = true;
    int count_coarsen = 1;
    for (int i=1; i<n; i++) {
      if (is_sibling_[i]) {
        if (can_coarsen_[i]) ++count_coarsen;
        if (level_min_[i] >= level_curr_[0]) cant_coarsen = true;
      }
    }
    if (count_coarsen != cello::num_children()) {
      // if not known if can coarsen yet, reset max to curr
      level_max_[0] = level_curr_[0];
    }
    if (cant_coarsen) {
      // if cant coarsen, update level_min
      level_min_[0] = level_curr_[0];
    }
  }

  return ( (level_min_[0] != level_min) ||
           (level_max_[0] != level_max) ||
           (can_coarsen_[0] != can_coarsen) );
}

//----------------------------------------------------------------------

bool Adapt::is_committed(int i) const
{
  return (level_min_.at(i) == level_max_.at(i));
}

//----------------------------------------------------------------------

void Adapt::get_level_bounds
(int * level_min, int * level_max, bool * can_coarsen, int i) const
{
  (*level_min)   = level_min_[i];
  (*level_max)   = level_max_[i];
  (*can_coarsen) = can_coarsen_[i];
}

//======================================================================

