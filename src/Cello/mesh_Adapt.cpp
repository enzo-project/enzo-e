// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Array.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-10-27
/// @brief    Implementation of the Adapt class for remeshing

#include "mesh.hpp"

//----------------------------------------------------------------------

void Adapt::set_face_level (const int if3[3], LevelType level_type,
                            int level_min, int level_max, bool can_coarsen)
{
  set_face_level(if3,level_type,level_min);
}

//----------------------------------------------------------------------

void Adapt::set_face_level_last (const int ic3[3], const int if3[3],
                                 int level_min, int level_max, bool can_coarsen)
{
  set_face_level_last(ic3,if3,level_min);
}

//----------------------------------------------------------------------

bool Adapt::insert_neighbor (Index index, bool is_sibling)
{
  const bool found = is_neighbor(index);

  if (! found) {
    const int level = index.level();

    LevelInfo neighbor { index, level, level-1, level+1, is_sibling, false };
    neighbor_list_.push_back(neighbor);
  }

  return (! found);
}

//----------------------------------------------------------------------

bool Adapt::delete_neighbor (Index index)
{
  int i;
  const bool found = is_neighbor(index,&i);

  // If found...
  if (found) {
    delete_neighbor_(i);
  }
  return found;
}

//----------------------------------------------------------------------

void Adapt::reset_bounds()
{
  // reset self level bounds
  self_.level_min_ = std::max(self_.level_now_-1,min_level_);
  self_.level_max_ = std::min(self_.level_now_+1,max_level_);
  self_.can_coarsen_ = false;
  // reset neighbor level bounds
  const int n = num_neighbors();
  for (int i=0; i<n; i++) {
    LevelInfo & neighbor = neighbor_list_[i];
    neighbor.level_min_ = std::max(neighbor.level_now_-1, min_level_);
    neighbor.level_max_ = std::min(neighbor.level_now_+1, max_level_);
    neighbor.can_coarsen_ = false;
  }
}

//----------------------------------------------------------------------

void Adapt::refine_neighbor (Index index)
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
          int adj = self_.index_.adjacency(index_child,rank_,periodicity_);
          if (adj >= 0) {
            bool s = insert_neighbor(index_child);
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void Adapt::coarsen_neighbor (Index index)
{
  // can't coarsen if neighbor sibling
  if ( ! self_.index_.is_sibling(index)) {
    int n = num_neighbors();
    // delete all neighbors contained in index parent
    // (note: use copy of neighbor list to avoid issues with looping
    // over data-structure that is actively being modified)
    auto neighbors = index_neighbors_();
    for (int i=0; i<n; i++) {
      const Index index_sibling = neighbors[i];
      if (index.is_sibling(index_sibling)) {
        delete_neighbor (index_sibling);
      }
    }
    // insert parent (if not already inserted)
    Index index_parent = index.index_parent(min_level_);
    insert_neighbor (index_parent);
  }
}

//----------------------------------------------------------------------

void Adapt::refine (const Adapt & adapt_parent, int ic3[3])
{
  Index index_parent = adapt_parent.index();

  self_.index_ = index_parent.index_child(ic3);

  // insert adjacent neighbors
  const int np = adapt_parent.neighbor_list_.size();
  for (int i=0; i<np; i++) {
    const Index index_neighbor = adapt_parent.neighbor_list_[i].index_;
    const int adj = self_.index_.adjacency(index_neighbor,rank_,periodicity_);
    if (adj >= 0) {
      insert_neighbor(index_neighbor);
    }
  }

  // add siblings
  const int level = index_parent.level();
  Index index_sibling = self_.index_;

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
  // Delete parent if it is a neighbor (e.g. periodic
  // b.c. with blocking == 1 along an axis
  delete_neighbor(index_parent);
}

//----------------------------------------------------------------------

void Adapt::coarsen (const Adapt & adapt_child)
{
  Index index_child = adapt_child.index();

  if (! valid_) {
    // First time called copy all child neighbors
    copy_(adapt_child);
    valid_ = true;
    self_.index_ = self_.index_.index_parent(min_level_);

    // remove any blocks internal to this one (siblings to child block)
    const int n = adapt_child.neighbor_list_.size();
    for (int i=0; i<n; i++) {
      const Index index_neighbor = adapt_child.neighbor_list_[i].index_;
      if (index_child.is_sibling(index_neighbor)) {
        delete_neighbor (index_neighbor);
      }
    }
  } else {
    // not first time called--insert non-sibling neighbors of child block
    const int n = adapt_child.neighbor_list_.size();
    for (int i=0; i<n; i++) {
      const Index index_neighbor = adapt_child.neighbor_list_[i].index_;
      const bool not_sibling = ! index_child.is_sibling(index_neighbor);
      if (not_sibling) {
        insert_neighbor (index_neighbor);
      }
    }
  }
}


//----------------------------------------------------------------------

void Adapt::initialize_self
(Index index, int level_min, int level_now)
{
  self_.level_now_ = level_now;
  self_.level_min_ = std::max(level_min,min_level_);
  self_.level_max_ = std::min(level_now+1,max_level_);
  self_.can_coarsen_ = false;
}

//----------------------------------------------------------------------

void Adapt::update_neighbor
(Index index, int level_min, int level_max, bool can_coarsen)
{
  LevelInfo * neighbor = neighbor_(index);
  if (neighbor) {
    neighbor->level_min_ = std::max(neighbor->level_min_,level_min);
    neighbor->level_max_ = std::min(neighbor->level_max_,level_max);
    neighbor->can_coarsen_ = neighbor->can_coarsen_ || can_coarsen;
  }
}

//----------------------------------------------------------------------

bool Adapt::update_bounds()
{
  // Save values to test later if changed
  const int level_min = self_.level_min_;
  const int level_max = self_.level_max_;
  const int can_coarsen = self_.can_coarsen_;

  const int n = neighbor_list_.size();
  int level_min_new = self_.level_min_;
  int level_max_new = self_.level_max_;
  for (int i=0; i<n; i++) {
    level_min_new = std::max(level_min_new, (neighbor_list_[i].level_min_ - 1) );
  }
  int neighbor_max = min_level_;
  for (int i=0; i<n; i++) {
    neighbor_max = std::max(neighbor_max,neighbor_list_[i].level_max_);
  }

  level_max_new = std::max(level_min_new,neighbor_max-1);

  self_.can_coarsen_ = (level_max_new == self_.level_now_ - 1);

  // Reset level_max to >= level_now pending checking all sibling can
  // coarsen
  self_.level_min_ = level_min_new;
  self_.level_max_ = std::max(level_max_new,self_.level_now_);

  // adjust for coarsening: can only coarsen if all siblings can coarsen

  const bool want_to_coarsen = (self_.level_min_ < self_.level_now_);
  int count_coarsen = 0;
  bool cannot_coarsen = false;
  if ( want_to_coarsen && self_.can_coarsen_) {
    // count self
    ++count_coarsen;
    // block can coarsen, check that all neighbors can as well
    for (int i=0; i<n; i++) {
      const auto & neighbor = neighbor_list_[i];
      if (self_.index_.is_sibling(neighbor.index_)) {
        if (neighbor.can_coarsen_) {
          ++count_coarsen;
        }
        if (neighbor.level_min_ >= self_.level_now_) {
          cannot_coarsen = true;
        }
      } else if (self_.index_.is_nibling(neighbor.index_)) {
          // cannot coarsen if sibling is refined
          cannot_coarsen = true;
      }

    }
    if (count_coarsen == cello::num_children(rank_)) {
      // Can coarsen if all siblings can coarsen (count includes self)
      self_.level_max_ = self_.level_min_;
    }
    if (cannot_coarsen) {
      self_.level_min_ = self_.level_now_;
      self_.can_coarsen_ = false;
    }
  }

  const bool changed =
    (self_.level_min_ != level_min) ||
    (self_.level_max_ != level_max) ||
    (self_.can_coarsen_ != can_coarsen);

  return ( changed );
}

//----------------------------------------------------------------------

bool Adapt::is_converged() const
{
  return (self_.level_min_ == self_.level_max_);
}

//----------------------------------------------------------------------

bool Adapt::neighbors_converged() const
{
  const int n = neighbor_list_.size();
  bool converged = true;
  for (int i=0; i<n; i++) {
    bool is_self = self_.index_ == neighbor_list_[i].index_;
    const auto & neighbor = is_self ? self_ : neighbor_list_[i];
    converged = converged && 
       (neighbor.level_min_ == neighbor.level_max_);
  }
  return converged;
}

//----------------------------------------------------------------------

void Adapt::get_level_bounds
(int * level_min, int * level_max, bool * can_coarsen) const
{
  (*level_min)   = self_.level_min_;
  (*level_max)   = self_.level_max_;
  (*can_coarsen) = self_.can_coarsen_;
}

//----------------------------------------------------------------------

void Adapt::get_neighbor_level_bounds
(Index index, int * level_min, int * level_max, bool * can_coarsen) const
{
  for (int i=0; i<num_neighbors(); i++) {
    auto & neighbor = neighbor_list_[i];
    if (index == neighbor_list_[i].index_) {
      (*level_min)   = neighbor.level_min_;
      (*level_max)   = neighbor.level_max_;
      (*can_coarsen) = neighbor.can_coarsen_;
    }
  }
}

//----------------------------------------------------------------------

void Adapt::print(std::string message, Block * block) const
{
  char prefix[255];
  if (block) {
    sprintf (prefix,"%d PRINT_ADAPT %s %s %p",CkMyPe(),block->name().c_str(),message.c_str(),(void*)this);
  } else {
    sprintf (prefix,"PRINT_ADAPT %s",message.c_str());
  }
  CkPrintf ("DEBUG_ADAPT face_level curr: ");
  for (int i=0; i<face_level_[0].size(); i++) {
    CkPrintf ("%d ", face_level_[0].at(i));
  }
  CkPrintf ("\n");
  CkPrintf ("DEBUG_ADAPT face_level next: ");
  for (int i=0; i<face_level_[1].size(); i++) {
    CkPrintf ("%d ", face_level_[1].at(i));
  }
  CkPrintf ("\n");
  CkPrintf ("DEBUG_ADAPT periodicity %d %d %d\n",
            periodicity_[0],
            periodicity_[1],
            periodicity_[2]);
  for (int i=0; i<face_level_[2].size(); i++) {
    if (i%27==0) CkPrintf ("\nDEBUG_ADAPT face_level last: ");
    CkPrintf ("%2d ", face_level_[2].at(i));
  }
  CkPrintf ("\n");
  CkPrintf ("%s (%d <= %d <= %d) C%d\n",
            prefix,
            self_.level_min_,
            self_.level_now_,
            self_.level_max_,
            self_.can_coarsen_?1:0);

  const int n = num_neighbors();
  CkPrintf ("%s    num_neighbors %d\n",prefix,n);
  for (int i=0; i<n; i++) {
    const LevelInfo & info = neighbor_list_.at(i);
    int il3[3];
    info.index_.index_level(il3,max_level_);
    int level = info.index_.level();
    char neighbor_block[80];
    if (block) {
      sprintf (neighbor_block,"%s",block->name(info.index_).c_str());
    } else {
      int it3[3],ia3[3];
      info.index_.array(ia3,ia3+1,ia3+2);
      info.index_.tree(it3,it3+1,it3+2);
      sprintf (neighbor_block,"%X:%X %X:%X,%X:%X",
               ia3[0],it3[0],
               ia3[1],it3[1],
               ia3[2],it3[2]);
    }
    const int l = 1 << (max_level_ - level);
    CkPrintf ("%s   %d %s [%d %d] S%d C%d\n",
              prefix,
              i,neighbor_block,
              info.level_min_,
              info.level_max_,
              info.is_sibling_?1:0,
              info.can_coarsen_?1:0);
    fflush(stdout);
  }
}

//----------------------------------------------------------------------

void Adapt::write(std::string root, Block * block, int cycle_start)
{
  const int cycle = cello::simulation()->cycle();
  if (cycle >= cycle_start) {
    char filename[80];
    sprintf (filename,"%d-%s.%s",cycle,root.c_str(),block->name().c_str());

    CkPrintf ("%s\n",filename);

    FILE * fp = fopen (filename,"w");
    const int n = neighbor_list_.size();
    for (int i=0; i<n; i++) {
      fprintf (fp,"%s\n",
               block->name(neighbor_list_[i].index_).c_str());
    }
    fclose(fp);
  }
}

//======================================================================

int Adapt::data_size () const
{
  int size = 0;

  SIZE_VECTOR_TYPE(size,int,face_level_[0]);
  SIZE_VECTOR_TYPE(size,int,face_level_[1]);
  SIZE_VECTOR_TYPE(size,int,face_level_[2]);

  SIZE_SCALAR_TYPE(size,bool,valid_);
  SIZE_SCALAR_TYPE(size,int,rank_);
  SIZE_ARRAY_TYPE (size,int,periodicity_,3);
  SIZE_SCALAR_TYPE(size,int,min_level_);
  SIZE_SCALAR_TYPE(size,int,max_level_);
  SIZE_SCALAR_TYPE(size,LevelInfo,self_);
  SIZE_VECTOR_TYPE(size,LevelInfo,neighbor_list_);

  return size;
}

// ----------------------------------------------------------------------

char * Adapt::save_data (char * buffer) const
{
  char * pc = buffer;

  SAVE_VECTOR_TYPE(pc,int,face_level_[0]);
  SAVE_VECTOR_TYPE(pc,int,face_level_[1]);
  SAVE_VECTOR_TYPE(pc,int,face_level_[2]);

  SAVE_SCALAR_TYPE(pc,bool,valid_);
  SAVE_SCALAR_TYPE(pc,int,rank_);
  SAVE_ARRAY_TYPE (pc,int,periodicity_,3);
  SAVE_SCALAR_TYPE(pc,int,min_level_);
  SAVE_SCALAR_TYPE(pc,int,max_level_);
  SAVE_SCALAR_TYPE(pc,LevelInfo,self_);
  SAVE_VECTOR_TYPE(pc,LevelInfo,neighbor_list_);

  return pc;
}

// ----------------------------------------------------------------------

char * Adapt::load_data (char * buffer)
{
  char * pc = buffer;

  LOAD_VECTOR_TYPE(pc,int,face_level_[0]);
  LOAD_VECTOR_TYPE(pc,int,face_level_[1]);
  LOAD_VECTOR_TYPE(pc,int,face_level_[2]);

  LOAD_SCALAR_TYPE(pc,bool,valid_);
  LOAD_SCALAR_TYPE(pc,int,rank_);
  LOAD_ARRAY_TYPE (pc,int,periodicity_,3);
  LOAD_SCALAR_TYPE(pc,int,min_level_);
  LOAD_SCALAR_TYPE(pc,int,max_level_);
  LOAD_SCALAR_TYPE(pc,LevelInfo,self_);
  LOAD_VECTOR_TYPE(pc,LevelInfo,neighbor_list_);

  return pc;
}

//======================================================================

bool Adapt::is_neighbor (Index index, int * ip) const
{
  // Search for index
  int i;
  const int n = neighbor_list_.size();
  for (i=0; i<n; i++) {
    if (neighbor_list_[i].index_==index) break;
  }
  // return index if used
  if (ip) (*ip) = i;

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
  face_level_[0] = adapt.face_level_[0];
  face_level_[1] = adapt.face_level_[1];
  face_level_[2] = adapt.face_level_[2];

  valid_ =         adapt.valid_;
  rank_ =          adapt.rank_;
  periodicity_[0] = adapt.periodicity_[0];
  periodicity_[1] = adapt.periodicity_[1];
  periodicity_[2] = adapt.periodicity_[2];
  min_level_ =     adapt.min_level_;
  max_level_ =     adapt.max_level_;
  self_ =          adapt.self_;
  neighbor_list_ = adapt.neighbor_list_;
}

//----------------------------------------------------------------------

std::vector<Index> Adapt::index_neighbors_() const
{
  std::vector<Index> index_list;
  const int n = num_neighbors();
  index_list.resize(n);
  for (int i=0; i<n; i++) index_list[i] = index(i);
  return index_list;
}



