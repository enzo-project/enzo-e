// See LICENSE_CELLO file for license and copyright information

/// @file     particle_ParticleBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 17:37:19 PDT 2014
/// @brief    Implementation of the ParticleBlock class
///
/// A ParticleBlock object represents collections of particles of
/// multiple types, where the attributes associated with the types
/// are defined using a separate ParticleDescr class.
///
/// To conserve memory, particle attributes for all particles of a given 
/// type are "packed" into a single data_[index_type][] array.  The format
/// of this array is
///
/// valid  attribute-1       attribute-2   ...   attribute-k
///
///   V_1 A1_1 A2_1 A3_1 ... Ak_1  
///   V_2 A1_2 A2_2 A3_2 ... Ak_2  
///     ...
///   V_N A1_N A2_N A3_N ... Ak_N  
///
///  where Aj_k are bytes associated with the j'th attribute of particle k,
///  and V_j is a single byte that is non-zero if the particle is valid.
///  The number of bytes used by Aj_k is given by 
///  ParticleDescr::attribute_size(j)

#include "particle.hpp"

