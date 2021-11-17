// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Array.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-10-27
/// @brief    Implementation of the Adapt class for remeshing

#include "mesh.hpp"

//----------------------------------------------------------------------

bool Adapt::insert_neighbor (Index index, bool is_sibling)
{
  const bool found = is_neighbor_(index);

  if (! found) {
    index_map_[index] = neighbor_list_.size();
    const int level = index.level();

    NeighborInfo neighbor { index, level, level-1, level+1, is_sibling, false };
    neighbor_list_.push_back(neighbor);
  }

  return (! found);
}

//----------------------------------------------------------------------

bool Adapt::delete_neighbor (Index index)
{
  int i;
  const bool found = is_neighbor_(index,&i);

  // If found...
  if (found) {
    delete_neighbor_(i);
  }
  return found;
}

//----------------------------------------------------------------------

bool Adapt::refine_neighbor (Index index)
{
  const bool success = delete_neighbor(index);
  if (success) {
    const int cxp = 2;
    const int cyp = (rank_ >= 2) ? 2 : 1;
    const int czp = (rank_ >= 3) ? 2 : 1;
    for (int icz=0; icz<czp; icz++) {
      for (int icy=0; icy<cyp; icy++) {
        for (int icx=0; icx<cxp; icx++) {
          Index index_child = index.index_child(icx,icy,icz);
          int adj = index_.adjacency(index_child,rank_);
          if (adj >= 0) {
            bool s = insert_neighbor(index_child);
          }
        }
      }
    }
  }
  return success;
}

//----------------------------------------------------------------------

bool Adapt::coarsen_neighbor (Index index)
{
  const bool success = delete_neighbor(index);
  if (success) {
    Index index_parent = index.index_parent();
    const int n = neighbor_list_.size();
    // delete all neighbors contained in index parent
    for (int i=0; i<n; i++) {
      const Index & index_sibling = neighbor_list_[i].index_;
      if (index_sibling.index_parent() == index_parent) {
        delete_neighbor (index_sibling);
      }
    }
    // insert parent (if not already inserted)
    insert_neighbor (index_parent);
  }
  return false;
}

//----------------------------------------------------------------------

void Adapt::refine (const Adapt & adapt_parent, int ic3[3])
{
  if (rank_ == 0) {
    copy_(adapt_parent);
    // update index_
    index_ = index_.index_child(ic3);
  }

  Index index_parent = adapt_parent.index();

  // delete non-adjacent neighbors
  const int n = adapt_parent.neighbor_list_.size();
  for (int i=0; i<n; i++) {
    const Index & index_neighbor = adapt_parent.neighbor_list_[i].index_;
    const int adj = index_.adjacency(index_neighbor,rank_);
    if (adj == -1) {
      delete_neighbor(index_neighbor);
    }
  }

  // add siblings
  const int level = index_parent.level();
  Index index_sibling = index_;
  
  index_sibling.set_level(level+1);
  const int icyp = (rank_ >= 2 ? 2:1);
  const int iczp = (rank_ >= 3 ? 2:1);
  for (int icz=0; icz<iczp; icz++) {
    for (int icy=0; icy<icyp; icy++) {
      for (int icx=0; icx<2; icx++) {
        index_sibling.set_child(level+1,icx,icy,icz);
        if ( ! (icx==ic3[0] && icy== ic3[1] && icz==ic3[2])) {
          insert_neighbor(index_sibling);
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void Adapt::coarsen (const Adapt & adapt_child)
{
  if (rank_ == 0) {
    // First time called
    copy_(adapt_child);
    // update index_
    index_ = index_.index_parent();

    // delete previous siblings (neighbors whose parent is self)
    const int n = adapt_child.neighbor_list_.size();
    for (int i=0; i<n; i++) {
      const Index & index_neighbor = adapt_child.neighbor_list_[i].index_;
      if (index_neighbor.index_parent() == index_) {
        delete_neighbor (index_neighbor);
      }
    }
  } else {
    // Subsequent time called

    // insert non-siblings
    const int n = adapt_child.neighbor_list_.size();
    for (int i=0; i<n; i++) {
      const Index & index_neighbor = adapt_child.neighbor_list_[i].index_;
      if (index_neighbor.index_parent() != index_) {
        insert_neighbor (index_neighbor);
      }
    }
  }
}

//----------------------------------------------------------------------

void Adapt::initialize_self
(Index index, int level_min, int level_now, int level_max)
{
  level_now_ = level_now;
  level_min_ = level_min;
  level_max_ = level_max;
  i_can_coarsen_ = false;
}

//----------------------------------------------------------------------

void Adapt::update_neighbor
(Index index, int level_min, int level_max, bool can_coarsen)
{
  const int i = index_map_[index];
  neighbor_list_[i].level_min_ = std::max(neighbor_list_[i].level_min_,level_min);
  neighbor_list_[i].level_max_ = std::min(neighbor_list_[i].level_max_,level_max);
  neighbor_list_[i].can_coarsen_ = neighbor_list_[i].can_coarsen_ || can_coarsen;
}

//----------------------------------------------------------------------

bool Adapt::update_bounds ()
{
  // Save values to test later if changed
  int level_min = level_min_;
  int level_max = level_max_;
  int can_coarsen = i_can_coarsen_;
  
  const int n = neighbor_list_.size();
  for (int i=0; i<n; i++) {
    level_min_ = std::max(level_min_, (neighbor_list_[i].level_min_ - 1) );
  }
  int neighbor_max = 0;
  for (int i=0; i<n; i++) {
    neighbor_max = std::max(neighbor_max,neighbor_list_[i].level_max_);
  }
  level_max_ = std::max(level_min_,neighbor_max-1);

  // adjust for coarsening: can only coarsen if all siblings can coarsen

  const bool want_to_coarsen = (level_min_ < level_now_);
  if ( want_to_coarsen && is_committed() ) {
    // block can coarsen, check that all neighbors can as well
    bool cant_coarsen = false;
    i_can_coarsen_ = true;
    int count_coarsen = 1;
    for (int i=0; i<n; i++) {
     if (neighbor_list_[i].is_sibling_) {
        if (neighbor_list_[i].can_coarsen_) ++count_coarsen;
        if (neighbor_list_[i].level_min_ >= level_now_) cant_coarsen = true;
      }
    }
    if (count_coarsen != cello::num_children(rank_)) {
      // if not known if can coarsen yet, reset max to curr
      level_max_ = level_now_;
    }
    if (cant_coarsen) {
      // if cant coarsen, update level_min
      level_min_ = level_now_;
    }
  }

  return ( (level_min_ != level_min) ||
           (level_max_ != level_max) ||
           (i_can_coarsen_ != can_coarsen) );
}

//----------------------------------------------------------------------

bool Adapt::is_committed() const
{
  return (level_min_ == level_max_);
}

//----------------------------------------------------------------------

void Adapt::get_level_bounds
(int * level_min, int * level_max, bool * can_coarsen) const
{
  (*level_min)   = level_min_;
  (*level_max)   = level_max_;
  (*can_coarsen) = i_can_coarsen_;
}

//----------------------------------------------------------------------

void Adapt::print(std::string message) const
  
{
  CkPrintf ("DEBUG_ADAPT adapt.print %s\n",message.c_str());
  const int n = num_neighbors();
  CkPrintf ("DEBUG_ADAPT    num_neighbors %d\n",n);
  int v3[3];
  index_.values(v3);
  CkPrintf ("DEBUG_ADAPT L%d [%8X %8X %8X]\n",
            index_.level(),v3[0],v3[1],v3[2]);
  for (int i=0; i<n; i++) {
    const NeighborInfo & info = neighbor_list_.at(i);
    info.index_.values(v3);
    int il3[3];
    const int max_level = 2;
    info.index_.index_level(il3,max_level);
    int it3[3],ia3[3];
    info.index_.array(ia3,ia3+1,ia3+2);
    info.index_.tree(it3,it3+1,it3+2);
    int level = info.index_.level();
    const int l = 1 << (max_level - level);
    CkPrintf ("DEBUG_ADAPT   %d L%d (%d <= %d <= %d) S%d C%d  %d %d %d - %d %d %d\n",
              i,level,
              info.level_min_,
              info.level_now_,
              info.level_max_,
              info.is_sibling_?1:0,
              info.can_coarsen_?1:0,
              il3[0],il3[1],il3[2],
              il3[0]+l,il3[1]+l,il3[2]+l);
  }
}

//======================================================================

bool Adapt::is_neighbor_ (Index index, int * ip) const
{
  // Search for index
  int i;
  const int n = neighbor_list_.size();
  for (i=0; i<n; i++) {
    if (neighbor_list_[i].index_==index) break;
  }
  // return index if used
  if (ip) (*ip) = i;
  // Save whether found for return value
  return (i < n);
}

//----------------------------------------------------------------------

void Adapt::delete_neighbor_ (int i)
{
  // ... shift indices following index back one
  const int n = neighbor_list_.size();
  for (; i<n-1; i++) {
    neighbor_list_[i] = neighbor_list_[i+1];
  }
  // ... and resize
  neighbor_list_.resize(n-1);
}

//----------------------------------------------------------------------

void Adapt::copy_ (const Adapt & adapt)
{
  rank_ = adapt.rank_;
  level_now_ = adapt.level_now_;
  level_min_ = adapt.level_min_;
  level_max_ = adapt.level_max_;
  i_can_coarsen_ = adapt.i_can_coarsen_;
  index_set_ = adapt.index_set_;
  index_map_ = adapt.index_map_;
  neighbor_list_ = adapt.neighbor_list_;
  index_ = adapt.index_;
}

//----------------------------------------------------------------------
