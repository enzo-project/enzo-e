// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodMergeStars.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date
/// @brief

#ifndef ENZO_ENZO_METHOD_MERGESTARS
#define ENZO_ENZO_METHOD_MERGESTARS

class EnzoMethodMergeStars : public Method {

  /// @class   EnzoMethodMergeStars
  /// @ingroup Enzo
  /// @brief   [\ref Enzo] Encapsulate Merge Star Routines

public:

  // Create a new MergeStars object
  EnzoMethodMergeStars();

  /// Destructor
  virtual ~EnzoMethodMergeStars() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodMergeStars);

  /// Charm++ PUP::able migration constructor
  EnzoMethodMergeStars (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method
  virtual void compute( Block * block) throw();

  /// Name
  virtual std::string name () throw()
  { return "merge_stars"; }

  /// Not sure if this is needed
  virtual std::string particle_type () throw()
  { return "star";}

  // Compute the maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // methods

  void compute_(Block * block);

  void get_particle_coordinates_block_units_
                                (EnzoBlock * enzo_block, int it,
                                double * particle_coordinates_block_units,
                                double * merging_radius_block_units);

  bool particles_in_neighbouring_blocks_(double * particle_coordinates,
					int ** group_lists,int * group_sizes,
					int i);
  
protected: // attributes

  double merging_radius_cells_;
};

#endif /* EnzoMethodMergeStars */
