// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodRefresh.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2021-03-09
/// @brief    [\ref Problem] Declaration for the MethodRefresh class

#ifndef PROBLEM_METHOD_REFRESH_HPP
#define PROBLEM_METHOD_REFRESH_HPP

class MethodRefresh : public Method
{
  /// @class    MethodRefresh
  /// @ingroup  MethodRefresh
  /// @brief    [\ref MethodRefresh] Declaration of MethodRefresh
  ///
  /// Method for refreshing data in ghost zones

public: // interface

  /// Create a new MethodRefresh
  MethodRefresh
  (std::vector< std::string > field_list,
   std::vector< std::string > particle_list,
   int ghost_depth,
   int min_face_rank,
   bool all_fields,
   bool all_particles);

  /// Destructor
  virtual ~MethodRefresh() throw()
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodRefresh);

  /// Charm++ PUP::able migration constructor
  MethodRefresh (CkMigrateMessage *m)
    : Method(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    Method::pup(p);
    p | field_list_;
    p | particle_list_;
    p | ghost_depth_;
    p | min_face_rank_;
  }

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute ( Block * block) throw();

  /// Return the name of this MethodRefresh
  virtual std::string name () throw ()
  { return "refresh"; }

protected: // functions

  void refresh_ (Block * block);
  
protected: // attributes

  /// List of id's of fields to refresh
  std::vector<int> field_list_;

  /// List of id's of particle types to refresh
  std::vector<int> particle_list_;

  /// Ghost layer depth
  int ghost_depth_;

  /// Minimum dimensional face to refresh (0 corners, 1 edges, 2 facets)
  int min_face_rank_;

  /// Whether to refresh all fields, ignoring field_list_
  bool all_fields_;

  /// Whether to refresh all particles, ignoring particle_list_
  bool all_particles_;
};


#endif /* PROBLEM_METHOD_REFRESH_HPP */
