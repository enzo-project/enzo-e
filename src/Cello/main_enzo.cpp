#define CHARM_ENZO

#include "enzo.hpp"

//----------------------------------------------------------------------

CkReduction::reducerType r_method_turbulence_type;

extern CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs);

//--------------------------------------------------

void register_method_turbulence(void)
{
  r_method_turbulence_type = CkReduction::addReducer(r_method_turbulence); 
}

//--------------------------------------------------

// SEE enzo_EnzoMethodTurbulence.cpp for context
CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs)
{
  double accum[MAX_TURBULENCE_ARRAY] = { 0.0 };
  accum[INDEX_TURBULENCE_minD] = std::numeric_limits<double>::max();
  accum[INDEX_TURBULENCE_maxD] = std::numeric_limits<double>::min();

  for (int i=0; i<n; i++) {
    double * values = (double *) msgs[i]->getData();
    for (int ig=0; ig<MAX_TURBULENCE_ARRAY-2; ig++) {
      accum [ig] += values[ig];
    }
    accum [INDEX_TURBULENCE_minD] = 
      std::min(accum[INDEX_TURBULENCE_minD],values[INDEX_TURBULENCE_minD]);
    accum [INDEX_TURBULENCE_maxD] = 
      std::max(accum[INDEX_TURBULENCE_maxD],values[INDEX_TURBULENCE_maxD]);
  }
  return CkReductionMsg::buildNew(MAX_TURBULENCE_ARRAY*sizeof(double),accum);
}

//======================================================================

CkReduction::reducerType r_method_gravity_cg_type;

extern CkReductionMsg * r_method_gravity_cg(int n, CkReductionMsg ** msgs);

//--------------------------------------------------

void register_method_gravity_cg(void)
{
  r_method_gravity_cg_type = CkReduction::addReducer(r_method_gravity_cg); 
}

//--------------------------------------------------

// SEE enzo_EnzoMethodGravityCg.cpp for context
CkReductionMsg * r_method_gravity_cg(int n, CkReductionMsg ** msgs)
{
  long double accum[3] = { 0.0, 0.0, 0.0 };

  for (int i=0; i<n; i++) {
    long double * values = (long double *) msgs[i]->getData();
    for (int ig=0; ig<3; ig++) {
      accum [ig] += values[ig];
    }
  }
  return CkReductionMsg::buildNew(3*sizeof(long double),accum);
}

//======================================================================

#include "main.cpp"
