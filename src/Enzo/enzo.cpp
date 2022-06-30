#include "enzo.hpp"

namespace enzo {

  EnzoSimulation * simulation()
  {
    return proxy_enzo_simulation.ckLocalBranch();
  }

  const EnzoFactory * factory()
  {
    return (const EnzoFactory *) simulation()->factory();
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

  const EnzoMethodGrackle * grackle_method()
  {
    if (!enzo::config()->method_grackle_use_grackle) {return NULL;}
    return (const EnzoMethodGrackle *) problem()->method("grackle");
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

  bool uses_dual_energy_formalism(bool default_ret /* false */)
  {
    // TODO(mabruzzo): this is meant to be a short-term solution. My immediate
    // priority is to create a Physics object for storing EOS properties
    // (including dual energy formalism parameters). This function will be
    // modified or deleted at that time.

    EnzoProblem* prob = problem();
    if (prob->method("ppm")){
      return config()->ppm_dual_energy;
    } else if (prob->method("mhd_vlct")){
      return config()->method_vlct_dual_energy;
    } else if (prob->method("ppml")){
      return false;
    } else {
      return default_ret;
    }
  }

}
