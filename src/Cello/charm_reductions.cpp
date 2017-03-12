#include "charm.hpp"

//======================================================================

CkReduction::reducerType r_reduce_performance_type;

CkReductionMsg * r_reduce_performance(int n, CkReductionMsg ** msgs)
{
  return NULL;
  // NOT IMPLEMENTED YET
}

void register_reduce_performance(void)
{ r_reduce_performance_type = CkReduction::addReducer(r_reduce_performance); }


//======================================================================

CkReduction::reducerType sum_long_double_type;

void register_sum_long_double(void)
{ sum_long_double_type = CkReduction::addReducer(sum_long_double); }

CkReductionMsg * sum_long_double(int n, CkReductionMsg ** msgs)
{
  long double accum = 0.0;

  for (int i=0; i<n; i++) {
    long double * values = (long double *) msgs[i]->getData();
    accum += values[0];
  }
  return CkReductionMsg::buildNew(sizeof(long double),&accum);
}

//------------------------------------------------------------------------


CkReduction::reducerType sum_long_double_2_type;

void register_sum_long_double_2(void)
{ sum_long_double_2_type = CkReduction::addReducer(sum_long_double_2); }

CkReductionMsg * sum_long_double_2(int n, CkReductionMsg ** msgs)
{
  long double accum[2] = { 0.0, 0.0 };

  for (int i=0; i<n; i++) {
    long double * values = (long double *) msgs[i]->getData();
    accum [0] += values[0];
    accum [1] += values[1];
  }
  return CkReductionMsg::buildNew(2*sizeof(long double),accum);
}

//------------------------------------------------------------------------


CkReduction::reducerType sum_long_double_3_type;

void register_sum_long_double_3(void)
{ sum_long_double_3_type = CkReduction::addReducer(sum_long_double_3); }

CkReductionMsg * sum_long_double_3(int n, CkReductionMsg ** msgs)
{
  long double accum[3] = { 0.0, 0.0, 0.0 };
  for (int i=0; i<n; i++) {
    long double * values = (long double *) msgs[i]->getData();
    accum [0] += values[0];
    accum [1] += values[1];
    accum [2] += values[2];
  }
  return CkReductionMsg::buildNew(3*sizeof(long double),accum);
}

//======================================================================

CkReduction::reducerType sum_long_double_4_type;

void register_sum_long_double_4(void)
{ sum_long_double_4_type = CkReduction::addReducer(sum_long_double_4); }

CkReductionMsg * sum_long_double_4(int n, CkReductionMsg ** msgs)
{
  long double accum[4] = { 0.0, 0.0, 0.0, 0.0 };

  for (int i=0; i<n; i++) {
    long double * values = (long double *) msgs[i]->getData();
    accum [0] += values[0];
    accum [1] += values[1];
    accum [2] += values[2];
    accum [3] += values[3];
  }
  return CkReductionMsg::buildNew(4*sizeof(long double),accum);
}

//======================================================================

