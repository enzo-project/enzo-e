// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOrderMorton.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-01
/// @brief    [\ref Problem] Declaration of the MethodOrderMorton class for
///           generating the Morton (Z-)ordering of blocks in the hierarchy

#ifndef PROBLEM_METHOD_ORDER_MORTON_HPP
#define PROBLEM_METHOD_ORDER_MORTON_HPP

class MethodOrderMorton : public Method {

  /// @class    MethodOrderMorton
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  MethodOrderMorton(int min_level) throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodOrderMorton);
  
  /// Charm++ PUP::able migration constructor
  MethodOrderMorton (CkMigrateMessage *m)
    : Method (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    Method::pup(p);
    p | is_index_;
    p | is_count_;
    p | is_ratio_;
    p | is_next_;
    p | is_weight_;
    p | is_weight_child_;
    p | is_sync_index_;
    p | is_sync_weight_;
    p | min_level_;
  }

  void compute_continue( Block * block);
  void compute_complete( Block * block);
  void send_weight(Block * block, int weight, bool self);
  void recv_weight(Block * block, int ic3[3], int weight, bool self);
  void send_index(Block * block, int index, int count, bool self);
  void recv_index(Block * block, int index, int count, bool self);

public: // virtual methods
  
  /// Apply the method to determine the Morton ordering of blocks
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "order_morton"; }

private: // methods

  /// Return the pointer to the Block's Morton ordering index 
  long long * pindex_(Block * block);

  /// Return the pointer to the number of Block indices
  long long * pcount_(Block * block);

  /// Return the ratio of index to count
  double * pratio_(Block * block);
  
  /// Return the pointer to the Index of the "next" block
  Index * pnext_(Block * block);

  /// Return the pointer to the Block's weight (including self)
  long long * pweight_(Block * block);

  /// Return the pointer to the given Block's child weight
  long long * pweight_child_(Block * block, int index);

  /// Return the pointer to the Block's Morton ordering index 
  Sync * psync_index_(Block * block);

  /// Return the pointer to the Block's weight (including self)
  Sync * psync_weight_(Block * block);

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Block Scalar<int> index
  int is_index_;
  /// Block Scalar<int> count
  int is_count_;
  /// Block Scalar<double> ratio
  int is_ratio_;
  /// Block Scalar<Index> next
  int is_next_;
  /// Block Scalar<int> weight (#decendent blocks + self)
  int is_weight_;
  /// Block Scalar<int> child weight (array of size cello::num_children())
  int is_weight_child_;
  /// Block Scalar<Sync> sync counter for index (coarse->fine)
  int is_sync_index_;
  /// Block Scalar<sync> sync counter for weight (fine->coarse)
  int is_sync_weight_;

  /// Minimum refinement level for ordering; may be < 0
  int min_level_;
};

#endif /* PROBLEM_METHOD_ORDER_MORTON_HPP */

