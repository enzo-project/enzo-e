#include "charm.hpp"

//======================================================================

CkReduction::reducerType r_reduce_performance_type;

void register_reduce_performance(void)
{ r_reduce_performance_type = CkReduction::addReducer(r_reduce_performance); }

CkReductionMsg * r_reduce_performance(int n, CkReductionMsg ** msgs)
{
  if (n <= 0) return NULL;

  long long num_sum = ((long long*) (msgs[0]->getData()))[0];
  long long num_max = ((long long*) (msgs[0]->getData()))[1];

  const int length = 2 + num_sum + num_max;
  std::vector<long long> accum;
  ASSERT1 ("r_reduce_performance",
	   "Sanity check failed on expected accumulator array %d",
	   length, (length < 500));
  accum.resize(length);
  accum.clear();

  // save length
  accum [0] = num_sum;
  accum [1] = num_max;

  // sum remaining values
  for (int i=0; i<n; i++) {
    ASSERT4("r_reduce_performance()",
	    "CkReductionMsg actual size %d is different from expected %lu num_sum %d num_max %d",
	    msgs[i]->getSize(),length*sizeof(long long),num_sum,num_max,
	    (((long unsigned)(msgs[i]->getSize()) ==
              length*sizeof(long long))));
      
    long long * values = (long long *) msgs[i]->getData();
    int j = 2;
    for (int count=1; count <= num_sum; count++) {
      accum [j] += values[j];
      ++j;
    }
    for (int count=1; count <= num_max; count++) {
      accum [j] = std::max(accum[j],values[j]);
      ++j;
    }
  }

  return CkReductionMsg::buildNew(length*sizeof(long long),&accum[0]);
}

//======================================================================

CkReduction::reducerType r_reduce_method_debug_type;

void register_reduce_method_debug(void)
{ r_reduce_method_debug_type = CkReduction::addReducer(r_reduce_method_debug); }

CkReductionMsg * r_reduce_method_debug(int n, CkReductionMsg ** msgs)
{
  if (n <= 0) return NULL;

  long double length = ((long double*) (msgs[0]->getData()))[0];
  int num_fields = (length - 1) / 4;
  std::vector<long double> accum;
  accum.resize(length);
  accum.clear();

  // save length
  accum [0] = num_fields;

  // initialize reductions min max sum count
  int j=1;
  for (int i_f=0; i_f < num_fields; i_f++) {
    accum [j] = std::numeric_limits<long double>::max();
    ++j;
    accum [j] = -std::numeric_limits<long double>::max();
    ++j;
    accum [j] = 0.0;
    ++j;
    accum [j] = 0.0;
    ++j;
  }
  for (int i=0; i<n; i++) {
    long double * values = (long double *) msgs[i]->getData();
    j = 1;
    for (int i_f=0; i_f < num_fields; i_f++) {
      accum [j] = std::min(accum[j],values[j]); ++j;
      accum [j] = std::max(accum[j],values[j]); ++j;
      accum [j] += values[j]; ++j;
      accum [j] += values[j]; ++j;
    }
  }

  return CkReductionMsg::buildNew(length*sizeof(long double),&accum[0]);
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
	    "CkReductionMsg actual size %d is different from expected %lu",
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
  long double accum[2];
  std::fill_n(accum,2,0.0);

  for (int i=0; i<n; i++) {

    ASSERT2("sum_long_double_2()",
	    "CkReductionMsg actual size %d is different from expected %lu",
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
  long double accum[3];
  std::fill_n(accum,3,0.0);
  
  for (int i=0; i<n; i++) {

    ASSERT2("sum_long_double_3()",
	    "CkReductionMsg actual size %d is different from expected %lu",
	    msgs[i]->getSize(),3*sizeof(long double),
	    (msgs[i]->getSize() == 3*sizeof(long double)));

    long double * values = (long double *) msgs[i]->getData();

    accum [0] += values[0];
    accum [1] += values[1];
    accum [2] += values[2];
  }
  return CkReductionMsg::buildNew(3*sizeof(long double),accum);
}

//----------------------------------------------------------------------

CkReduction::reducerType sum_long_double_n_type;

void register_sum_long_double_n(void)
{ sum_long_double_n_type = CkReduction::addReducer(sum_long_double_n); }

CkReductionMsg * sum_long_double_n(int n, CkReductionMsg ** msgs)
{

  const int N = int(*((long double *) msgs[0]->getData()));
  
  long double * accum = new long double [N+1];
  
  std::fill_n(accum,N+1,0.0);

  for (int i=0; i<n; i++) {

    ASSERT2("sum_long_double_n()",
	    "CkReductionMsg actual size %d is different from expected %lu",
	    msgs[i]->getSize(),(N+1)*sizeof(long double),
            ((long unsigned)msgs[i]->getSize() == (N+1)*sizeof(long double)));
    
    long double * values = (long double *) msgs[i]->getData();

    accum[0] = values[0];
    for (int i=1; i<=N; i++) {
      accum[i] += values[i];
    }
  }
  CkReductionMsg * msg = CkReductionMsg::buildNew
    ((N+1)*sizeof(long double),accum);
  delete [] accum;
  return msg;
}

//======================================================================

