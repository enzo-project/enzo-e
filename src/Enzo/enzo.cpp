#include "enzo.hpp"

namespace enzo {

  EnzoSimulation * simulation()
  {
    return proxy_enzo_simulation.ckLocalBranch();
  }

  EnzoProblem * problem()
  {
    return (EnzoProblem *) simulation()->problem();
  }

  const EnzoConfig * config()
  {
    return (const EnzoConfig *) simulation()->config();
  }

  EnzoPhysicsCosmology * cosmology()
  {
    return (EnzoPhysicsCosmology *) problem()->physics("cosmology");
  }

  EnzoUnits * units()
  {
    return (EnzoUnits *) problem()->units();
  }

  CProxy_EnzoBlock block_array()
  {
    return (CProxy_EnzoBlock) enzo::simulation()->hierarchy()->block_array();
  }

  EnzoBlock * block ( Block * block)
  {
    return static_cast<EnzoBlock*> (block);
  }

}
