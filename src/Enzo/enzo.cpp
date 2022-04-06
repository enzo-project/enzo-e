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

  void check_particle_attribute(std::string type,
				std::string attribute){
    
  ParticleDescr * particle_descr = cello::particle_descr();
  ASSERT1("enzo::check_particle_attribute",
	  "Error: This problem requires %s particle_type.",
	  type, particle_descr->type_exists(type));

  int it = particle_descr->type_index(type);
  ASSERT2("enzo::check_particle_attribute",
	  "Error: This problem requires %s particles to have an "
	  "attribute called %s.",
	  type, attribute,
	  particle_descr->has_attribute(it,attribute));
  }
}
