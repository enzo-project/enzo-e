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

}
