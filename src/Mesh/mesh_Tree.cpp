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
    num_nodes_(1)
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
  
  num_nodes_ += count;
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

  int k,kx,ky,kz;
  k = neighbor_trace->index();
  index_(k,&kx,&ky,&kz);
  std::stack<int> index_stack;
  int level = neighbor_trace->level();
  int level0 = node_trace.level();
  //--------------------------------------------------
  if        (ix == -1) {
    index_stack.push(k);
    while ( (level-- >= 0) && (kx == 0) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    neighbor_trace->pop();
    k -= dx;
    while ( (level++ < level0) && 
	    neighbor_trace->node()->child(k) != NULL) 
      {
	neighbor_trace->push(k);
	k = index_stack.top();
	k += dx;
	index_stack.pop();
      }
  //--------------------------------------------------
  } else if (ix == +1) {
    index_stack.push(k);
    while ( (level-- >= 0) && (kx == r_-1) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    neighbor_trace->pop();
    k += dx;
    while ( (level++ < level0) && 
	    neighbor_trace->node()->child(k) != NULL) 
      {
	neighbor_trace->push(k);
	k = index_stack.top();
	k -= dx;
	index_stack.pop();
      }
  //--------------------------------------------------
  } else if (iy == -1) {
    index_stack.push(k);

    while ( (level-- >= 0) && (ky == 0) )
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    neighbor_trace->pop();
    k -= dy;
    while ( (level++ < level0) && 
	    neighbor_trace->node()->child(k) != NULL) 
      {
	neighbor_trace->push(k);
	index_stack.pop();
	k = index_stack.top();
	k += dy;
      }
  //--------------------------------------------------
  } else if (iy == +1) {
    index_stack.push(k);
    while ( (level-- >= 0) && (ky == r_-1) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    neighbor_trace->pop();
    k += dy;
    while ( (level++ < level0) && 
	    neighbor_trace->node()->child(k) != NULL) 
      {
	neighbor_trace->push(k);
	k = index_stack.top();
	k -= dy;
	index_stack.pop();
      }
  //--------------------------------------------------
  } else if (iz == -1) {
    index_stack.push(k);
    while ( (level-- >= 0) && (kz == 0) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    neighbor_trace->pop();
    k -= dz;
    while ( (level++ < level0) && 
	    neighbor_trace->node()->child(k) != NULL) 
      {
	neighbor_trace->push(k);
	k = index_stack.top();
	k += dz;
	index_stack.pop();
      }
  //--------------------------------------------------
  } else if (iz == +1) {
    index_stack.push(k);
    while ( (level-- >= 0) && (kz == r_-1) ) 
      {
	neighbor_trace->pop();
	k = neighbor_trace->index();
	index_(k,&kx,&ky,&kz);
	index_stack.push(k);
      }
    neighbor_trace->pop();
    k += dz;
    while ( (level++ < level0) && 
	    neighbor_trace->node()->child(k) != NULL) 
      {
	neighbor_trace->push(k);
	k = index_stack.top();
	k -= dz;
	index_stack.pop();
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

void Tree::balance_node (const NodeTrace & node_trace)
{

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

//======================================================================

