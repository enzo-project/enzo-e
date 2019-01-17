// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMaker.hpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
/// @date
/// @brief

#ifndef ENZO_ENZO_METHOD_FEEDBACK
#define ENZO_ENZO_METHOD_FEEDBACK

class EnzoMethodFeedback : public Method {

  /// @class   EnzoMethodFeedback
  /// @ingroup Enzo
  /// @btief   [\ref Enzo] Encapsulate Feedback Routines

public:

  EnzoMethodFeedback();

  /// Destructor
  virtual ~EnzoMethodFeedback() throw() {};

  /// CHarm++ Pup::able declarations
  PUPable_decl(EnzoMethodFeedback);

  /// Charm++ Pup::able migration Constructor
  EnzoMethodFeedback (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// Charm++ Pack / Unpack function
  void pup(PUP::er &p);

  /// Apply the method
  virtual void compute (Block * block) throw();

  /// name
  virtual std::string name() throw()
  { return "feedback"; }

  // Compute the maximum timestep for this method
  virtual double timestep (Block * block) const throw();

// protected::

};


#endif /* EnzoMethodFeedback */
