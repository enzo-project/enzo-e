#include "charm.hpp"

//======================================================================

CkReduction::reducerType r_reduce_performance_type;

void register_reduce_performance(void)
{ r_reduce_performance_type = CkReduction::addReducer(r_reduce_performance); }

CkReductionMsg * r_reduce_performance(int n, CkReductionMsg ** msgs)
{
  if (n <= 0) return NULL;

  long long length = ((long long*) (msgs[0]->getData()))[0];

  long long * accum = new long long [length];
  for (int i=0; i<length; i++) accum[i] = 0.0;

  // save length
  accum [0] = length;

  // sum remaining values
  for (int i=0; i<n; i++) {
    ASSERT2("r_reduce_performance()",
	    "CkReductionMsg actual size %d is different from expected %d",
	    msgs[i]->getSize(),length*sizeof(long long),
	    (msgs[i]->getSize() == length*sizeof(long long)));
      
    long long * values = (long long *) msgs[i]->getData();
    for (int j=1; j<length; j++) {
      accum [j] += values[j];
    }
  }
  return CkReductionMsg::buildNew(length*sizeof(long long),accum);
}


//======================================================================

CkReduction::reducerType sum_long_double_type;

void register_sum_long_double(void)
{ sum_long_double_type = CkReduction::addReducer(sum_long_double); }

CkReductionMsg * sum_long_double(int n, CkReductionMsg ** msgs)
{
  long double accum = 0.0;

  for (int i=0; i<n; i++) {
    ASSERT2("sum_long_double()",
	    "CkReductionMsg actual size %d is different from expected %d",
	    msgs[i]->getSize(),sizeof(long double),
	    (msgs[i]->getSize() == sizeof(long double)));
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
    ASSERT2("sum_long_double_2()",
	    "CkReductionMsg actual size %d is different from expected %d",
	    msgs[i]->getSize(),2*sizeof(long double),
	    (msgs[i]->getSize() == 2*sizeof(long double)));
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
    ASSERT2("sum_long_double_3()",
	    "CkReductionMsg actual size %d is different from expected %d",
	    msgs[i]->getSize(),3*sizeof(long double),
	    (msgs[i]->getSize() == 3*sizeof(long double)));
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

    ASSERT2("sum_long_double_4()",
	    "CkReductionMsg actual size %d is different from expected %d",
	    msgs[i]->getSize(),4*sizeof(long double),
	    (msgs[i]->getSize() == 4*sizeof(long double)));
    
    long double * values = (long double *) msgs[i]->getData();
    
    accum [0] += values[0];
    accum [1] += values[1];
    accum [2] += values[2];
    accum [3] += values[3];
  }
  return CkReductionMsg::buildNew(4*sizeof(long double),accum);
}

//======================================================================

