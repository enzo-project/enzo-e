// See LICENSE_CELLO file for license and copyright information

/// @file     problem_InitialTrace.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2015-11-06 23:00:30
/// @brief    [\ref Problem] Declaration for the InitialTrace component

#ifndef PROBLEM_INITIAL_TRACE_HPP
#define PROBLEM_INITIAL_TRACE_HPP

class Hierarchy;
class Field;
class Particle;

class InitialTrace : public Initial
{

  /// @class    InitialTrace
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Initialize trace particles

public: // interface

  /// empty constructor for charm++ pup()
  InitialTrace(std::string field,
	       double mpp,
	       int dx, int dy, int dz) throw() 
  : Initial(),
    mpp_(mpp),
    field_(field),
    dx_(dx),dy_(dy),dz_(dz)
  {}

  /// Destructor
  virtual ~InitialTrace() throw()
  {} ;

  /// CHARM++ PUP::able declaration
  PUPable_decl(InitialTrace);

  /// CHARM++ migration constructor for PUP::able
  InitialTrace (CkMigrateMessage *m)
    : Initial(m),
      mpp_(0.0),
      field_(""),
      dx_(0),dy_(0),dz_(0)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);


public: // virtual functions

  /// InitialTraceize a Block
  virtual void enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw();

protected: // functions

  /// Initial particle positions are a uniform array if mpp_ == 0
  void uniform_placement_ (Block * block, Field field, Particle particle);
  
  /// Initial particle positions are random based on local density
  void density_placement_ (Block * block, Field field, Particle particle);

public: // static attributes
  
  /// Next available particle ID, where id0_ = CkMyPe() + k * CkNumPes()
  static int id0_[CONFIG_NODE_SIZE];

protected: // attributes

  /// if mpp_ != 0, mass per particle--place on average one particle
  /// for each mpp grams
  double mpp_;

  /// Field to use when mpp_ != 0.
  std::string field_;

  /// if mpp_==0, put a particle in every cell ix, iy such that
  /// (ix%dx_==0) && (iy%dy_ == 0) && (iz%dz_ == 0)
  int dx_,dy_,dz_;

};

#endif /* PROBLEM_INITIAL_TRACE_HPP */
