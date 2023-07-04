// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoMethodMergeSinks.hpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   18 January 2022
/// @brief  Implementation of EnzoMethodMergeSinks class, a method for merging
///         sink particles separated by a distance less than the merging
///         radius.

#ifndef ENZO_ENZO_METHOD_MERGESINKS
#define ENZO_ENZO_METHOD_MERGESINKS

class EnzoMethodMergeSinks : public Method {

  /// @class   EnzoMethodMergeSinks
  /// @ingroup Enzo
  /// @brief   [\ref Enzo] Encapsulate Merge Sink Routines

public:

  // Create a new MergeSinks object
  EnzoMethodMergeSinks(double merging_radius_cells);

  /// Destructor
  virtual ~EnzoMethodMergeSinks() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodMergeSinks);

  /// Charm++ PUP::able migration constructor
  EnzoMethodMergeSinks (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method
  virtual void compute( Block * block) throw();

  /// Name
  virtual std::string name () throw()
  { return "merge_sinks"; }

  /// Not sure if this is needed
  virtual std::string particle_type () throw()
  { return "sink";}

  // Compute the maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // methods

  void compute_(Block * block);

  void get_particle_coordinates_(EnzoBlock * enzo_block, int it,
				enzo_float * particle_coordinates);

  bool particles_in_neighbouring_blocks_(EnzoBlock * enzo_block,
					 enzo_float * particle_coordinates,
					 int ** group_lists,int * group_sizes,
					 int i);

  // Checks to be performed at initial cycle
  void do_checks_(const Block* block) throw();

  
protected: // attributes

  enzo_float merging_radius_cells_;
};

#endif /* EnzoMethodMergeSinks */
