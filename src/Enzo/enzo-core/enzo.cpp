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

  EnzoPhysicsFluidProps * fluid_props()
  {
    Physics* out = problem()->physics("fluid_props");
    // handling in EnzoProblem::initialize_physics_coda_ should ensure that
    // this is never a nullptr
    ASSERT("enzo::fluid_props", "Something went wrong", out != nullptr);
    return (EnzoPhysicsFluidProps *) out;
  }

  const EnzoMethodGrackle * grackle_method()
  {
    // previously, this explicitly checked the value of
    // enzo::config()->method_grackle_use_grackle, but this was redundant
    return (const EnzoMethodGrackle *) problem()->method("grackle");
  }

  const GrackleChemistryData * grackle_chemistry()
  {
    const EnzoMethodGrackle* ptr = grackle_method();
    return (ptr != nullptr) ? ptr->try_get_chemistry() : nullptr;
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
