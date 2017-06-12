#include "enzo.hpp"

CkReduction::reducerType r_method_turbulence_type;

CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs)
{
  double accum[max_turbulence_array] = { 0.0 };
  accum[index_turbulence_mind] = std::numeric_limits<double>::max();
  accum[index_turbulence_maxd] = - std::numeric_limits<double>::max();

  for (int i=0; i<n; i++) {
    double * values = (double *) msgs[i]->getData();
    for (int ig=0; ig<max_turbulence_array-2; ig++) {
      accum [ig] += values[ig];
    }
    accum [index_turbulence_mind] = 
      std::min(accum[index_turbulence_mind],values[index_turbulence_mind]);
    accum [index_turbulence_maxd] = 
      std::max(accum[index_turbulence_maxd],values[index_turbulence_maxd]);
  }
  return CkReductionMsg::buildNew(max_turbulence_array*sizeof(double),accum);
}

void register_method_turbulence(void)
{ r_method_turbulence_type = CkReduction::addReducer(r_method_turbulence); }

