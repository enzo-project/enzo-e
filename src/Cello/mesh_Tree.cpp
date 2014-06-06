// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Tree.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "mesh.hpp"

//----------------------------------------------------------------------

Tree::Tree(int d, int r) throw ()
  : d_(d), r_(r), 
    root_(new Node), 
    num_nodes_(1),
    max_level_(0)
{
  c_ = 1;
  if (d_>=1) c_ *= r_;
  if (d_>=2) c_ *= r_;
  if (d_>=3) c_ *= r_;
}

//----------------------------------------------------------------------

Tree::~Tree() throw ()
{
  delete root_;
}

//----------------------------------------------------------------------

void Tree::refine_node (const NodeTrace & node_trace)
{
  int count =  node_trace.node()->refine(c_);
  
  int level = node_trace.level() + 1;
  num_nodes_ += count;
  if (level > max_level_) max_level_ = level;

}

//----------------------------------------------------------------------

void Tree::coarsen_node (const NodeTrace & node_trace)
{
  int count = node_trace.node()->coarsen(c_);

  num_nodes_ -= count;
}

//----------------------------------------------------------------------

bool Tree::node_neighbor 
(
 const NodeTrace & node_trace,
 NodeTrace * neighbor_trace,
 int ix,
 int iy,
 int iz) const
/// @param node_trace  Node trace of the node
/// @param neighbor_trace  Object defining the NodeTrace of the Node's neighbor
/// @param ix,iy,iz    Direction of the neighbor we wish to find
///
/// This function works by starting at the given node, traversing
/// back to a common ancestor, then tracing forward to the neighbor
/// node.  The function returns whether a neighbor is found (true) or
/// not (false).  (A neighbor may not be found it the node is on
/// the border.)
{

  *neighbor_trace = node_trace;

  if (ix==0 && iy==0 && iz==0) return neighbor_trace;

  if (d_ < 2) iy = 0;
  if (d_ < 3) iz = 0;

  ASSERT("Tree::node_neighbor",
	 "Can (currently) only find (d-1)-face neighbors",
	 (abs(ix) + abs(iy) + abs(iz)) <= 1);

  const int dx = 1;
  const int dy = r_;
  const int dz = r_*r_;

  int k,kx=0,ky=0,kz=0;
  k = neighbor_trace->index();
  index_(k,&kx,&ky,&kz);
  std::stack<int> index_stack;
  int level = neighbor_trace->level();
  int level0 = node_trace.level();
  //--------------------------------------------------
  if        (ix == -1) {
    if (kx > 0) {
      k -= dx;
      neighbor_trace->pop();
      neighbor_trace->push(k);
      return true;
    }
    index_stack.push(k);
    while ( (level-- > 0) && (kx == 0) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    if (level < 0) return false;
    k -= dx;
    neighbor_trace->pop();
    neighbor_trace->push(k);
    level++;
    index_stack.pop();
    while ( (level++ < level0) && 
	    ! neighbor_trace->node()->is_leaf())
      {
	k = index_stack.top();
	index_stack.pop();
	k += dx;
	neighbor_trace->push(k);
      }
  //--------------------------------------------------
  } else if (ix == +1) {
    if (kx < r_ - 1) {
      k += dx;
      neighbor_trace->pop();
      neighbor_trace->push(k);
      return true;
    }
    index_stack.push(k);
    while ( (level-- > 0) && (kx == r_-1) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    if (level < 0) return false;
    k += dx;
    neighbor_trace->pop();
    neighbor_trace->push(k);
    level++;
    index_stack.pop();
    while ( (level++ < level0) && 
	    ! neighbor_trace->node()->is_leaf())
      {
	k = index_stack.top();
	index_stack.pop();
	k -= dx;
	neighbor_trace->push(k);
      }
  //--------------------------------------------------
  } else if (iy == -1) {
    if (ky > 0) {
      k -= dy;
      neighbor_trace->pop();
      neighbor_trace->push(k);
      return true;
    }
    index_stack.push(k);
    while ( (level-- > 0) && (ky == 0) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    if (level < 0) return false;
    k -= dy;
    neighbor_trace->pop();
    neighbor_trace->push(k);
    level++;
    index_stack.pop();
    while ( (level++ < level0) && 
	    ! neighbor_trace->node()->is_leaf())
      {
	k = index_stack.top();
	index_stack.pop();
	k += dy;
	neighbor_trace->push(k);
      }
  //--------------------------------------------------
  } else if (iy == +1) {
    if (ky < r_ - 1) {
      k += dy;
      neighbor_trace->pop();
      neighbor_trace->push(k);
      return true;
    }
    index_stack.push(k);
    while ( (level-- > 0) && (ky == r_-1) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    if (level < 0) return false;
    k += dy;
    neighbor_trace->pop();
    neighbor_trace->push(k);
    level++;
    index_stack.pop();
    while ( (level++ < level0) && 
	    ! neighbor_trace->node()->is_leaf())
      {
	k = index_stack.top();
	index_stack.pop();
	k -= dy;
	neighbor_trace->push(k);
      }
  //--------------------------------------------------
  } else if (iz == -1) {
    if (kz > 0) {
      k -= dz;
      neighbor_trace->pop();
      neighbor_trace->push(k);
      return true;
    }
    index_stack.push(k);
    while ( (level-- > 0) && (kz == 0) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    if (level < 0) return false;
    k -= dz;
    neighbor_trace->pop();
    neighbor_trace->push(k);
    level++;
    index_stack.pop();
    while ( (level++ < level0) && 
	    ! neighbor_trace->node()->is_leaf())
      {
	k = index_stack.top();
	index_stack.pop();
	k += dz;
	neighbor_trace->push(k);
      }
  //--------------------------------------------------
  } else if (iz == +1) {
    if (kz < r_ - 1) {
      k += dz;
      neighbor_trace->pop();
      neighbor_trace->push(k);
      return true;
    }
    index_stack.push(k);
    while ( (level-- > 0) && (kz == r_-1) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    if (level < 0) return false;
    k += dz;
    neighbor_trace->pop();
    neighbor_trace->push(k);
    level++;
    index_stack.pop();
    while ( (level++ < level0) && 
	    ! neighbor_trace->node()->is_leaf())
      {
	k = index_stack.top();
	index_stack.pop();
	k -= dz;
	neighbor_trace->push(k);
      }
  }

  return true;
}

//----------------------------------------------------------------------

void Tree::index_(int k, int * kx, int *ky, int *kz) const
{
  if (kx) (*kx) = k % r_;
  if (d_ > 1) {
    k -= (*kx);
    k /= r_;
    if (ky) (*ky) = k % r_;
    if (d_ > 2) {
      k -= (*ky);
      k /= r_;
      if (kz) (*kz) = k % r_;
    }
  }
}

//----------------------------------------------------------------------

NodeTrace Tree::node_parent (const NodeTrace & node_trace) const
{
  NodeTrace parent_trace (node_trace);
  parent_trace.pop();
  return parent_trace;
}

//----------------------------------------------------------------------

NodeTrace Tree::node_child (const NodeTrace & node_trace, int index) const
{
  NodeTrace child_trace (node_trace);
  child_trace.push(index);
  return child_trace;
}

//----------------------------------------------------------------------

void Tree::balance ()
{
  int dx[] = {-1, +1,  0,  0,  0,  0};
  int dy[] = { 0,  0, -1, +1,  0,  0};
  int dz[] = { 0,  0,  0,  0, -1, +1};

  for (int level = max_level(); level >= 0; --level) {
    ItNode it_node (this,level);
    while (it_node.next_leaf()) {
      const NodeTrace * node_trace = it_node.node_trace();
      NodeTrace neighbor_trace (root_);
      bool is_neighbor;

      for (int i = 0; i < 2*dimension(); i++) {
	bool do_refine;
	do {
	  is_neighbor = 
	    node_neighbor(*node_trace,&neighbor_trace,dx[i],dy[i],dz[i]);
	  do_refine = is_neighbor &&
	    neighbor_trace.level() < node_trace->level() - 1;
	  if (do_refine) {
	    refine_node(neighbor_trace);
	  }
	} while (do_refine);
      }
      
    }
    
  }
}

//----------------------------------------------------------------------

void Tree::coalesce ()
{
  int * level_data = new int[num_nodes()];

  ItNode it_node (this);

  // Initialize node levels

  int index_data=0;
  while (it_node.next_leaf()) {

    const NodeTrace * node_trace = it_node.node_trace();

    level_data[index_data] = node_trace->level();

    node_trace->node()->set_data(&level_data[index_data++]);
  }

  // Repeatedly coarsen nodes with children but no grandchildren
  int passes = 0;
  int r2d = 1;
  for (int i=0; i<dimension(); i++) r2d *= refinement();

  int count;

  do {

    count = 0;
    for (int level = max_level(); level >= 0; --level) {
      ItNode it_node (this,level);
      while (it_node.next_leaf()) {

	const NodeTrace * node_trace = it_node.node_trace();
	// actually want leaf parents, not leaves themselves
	if (node_trace->index()==0) {
	  int level = node_trace->level();
	  -- level;
	  Node * node = 0;
	  if (level > 0) {
	    bool do_coarsen = true;
	    node = node_trace->node_level(level);
	    // loop over leaf parent's children
	    int * level0 = (int *) node->child(0)->data();
	    for (int i=0; i<r2d; i++) {
	      Node * child = node->child(i);
	      if (child->is_leaf()) {
		int * level1 = (int * ) child->data();
		if (*level1 != *level0) {
		  do_coarsen = false;
		}
	      } else do_coarsen = false;
	    }
	    if (do_coarsen) {
	      count ++;
	      // adjust mesh level 
	      coarsen_node(node);
	      int * level1 = (int * ) node->data();
	      if (level1 == NULL) {
		node->set_data(&level_data[index_data++]);
		level1 = (int * ) node->data();
	      }
	      *level1 = *level0;
	    }
	  }
	}
      }  
    }
    passes ++;
    printf ("Pass %d  Count %d\n",passes,count);
  } while (count > 0);
  
}
//----------------------------------------------------------------------
//======================================================================

