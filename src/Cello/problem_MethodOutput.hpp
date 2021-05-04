// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOutput.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2021-03-09
/// @brief    [\ref Problem] Declaration for the MethodOutput class

#ifndef PROBLEM_METHOD_OUTPUT_HPP
#define PROBLEM_METHOD_OUTPUT_HPP

class IoBlock;
class IoFieldData;
class IoParticleData;
class FileHdf5;

class MethodOutput : public Method
{
  /// @class    MethodOutput
  /// @ingroup  MethodOutput
  /// @brief    [\ref MethodOutput] Declaration of MethodOutput
  ///
  /// Method for writing data to disk files.  Designed to be highly
  /// scalable.

public: // interface

  /// Create a new MethodOutput
  MethodOutput
  (std::vector< std::string > file_name,
   std::vector< std::string > path_name,
   std::vector< std::string > field_list,
   std::vector< std::string > particle_list,
   int ghost_depth,
   int min_face_rank,
   bool all_fields,
   bool all_particles,
   int blocking_x = 1,
   int blocking_y = 1,
   int blocking_z = 1);

  /// Destructor
  virtual ~MethodOutput() throw()
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodOutput);

  /// Charm++ PUP::able migration constructor
  MethodOutput (CkMigrateMessage *m)
    : Method(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    Method::pup(p);
    p | file_name_;
    p | path_name_;
    p | field_list_;
    p | particle_list_;
    p | ghost_depth_;
    p | min_face_rank_;
    p | all_fields_;
    p | all_particles_;
    PUParray(p,blocking_,3);
    p | is_count_;
  }

  void compute_continue (Block * block);
  void next (Block * block, MsgOutput *);
  void write (Block * block, MsgOutput *);

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute ( Block * block) throw();

  /// Return the name of this MethodOutput
  virtual std::string name () throw ()
  { return "output"; }

protected: // functions

  void output_ (Block * block);

  int is_writer_ (Index index);

  FileHdf5 * file_open_(Block * block, int a3[3]);
  void file_write_hierarchy_(FileHdf5 * file);
  void file_write_block_(FileHdf5 * , Block * , MsgOutput *);
  
  
protected: // attributes

  /// File name and format
  std::vector <std::string> file_name_;

  /// Path name and format
  std::vector <std::string> path_name_;
  
  /// List of id's of fields to output
  std::vector<int> field_list_;

  /// List of id's of particle types to output
  std::vector<int> particle_list_;

  /// Ghost layer depth
  int ghost_depth_;

  /// Minimum dimensional face to output (0 corners, 1 edges, 2 facets)
  int min_face_rank_;

  /// Whether to output all fields, ignoring field_list_
  bool all_fields_;

  /// Whether to output all particles, ignoring particle_list_
  bool all_particles_;

  /// Size of the root-level octree array partitioning.  Data in all
  /// blocks in each partition are written to a single file by the
  /// minimal root-level Block
  int blocking_[3];

  /// Block Scalar int for counting files
  int is_count_;
};


#endif /* PROBLEM_METHOD_OUTPUT_HPP */
