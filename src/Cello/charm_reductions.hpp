extern CkReduction::reducerType r_reduce_performance_type;
extern CkReductionMsg * r_reduce_performance(int n, CkReductionMsg ** msgs);
extern void register_reduce_performance(void);

extern CkReduction::reducerType sum_long_double_type;
extern CkReductionMsg * sum_long_double_2(int n, CkReductionMsg ** msgs);
extern void register_sum_long_double_2(void);

extern CkReduction::reducerType sum_long_double_2_type;
extern CkReductionMsg * sum_long_double_3(int n, CkReductionMsg ** msgs);
extern void register_sum_long_double_3(void);

extern CkReduction::reducerType sum_long_double_3_type;
extern CkReductionMsg * sum_long_double_4(int n, CkReductionMsg ** msgs);
extern void register_sum_long_double_4(void);

extern CkReduction::reducerType sum_long_double_4_type;
extern CkReductionMsg * sum_long_double(int n, CkReductionMsg ** msgs);
extern void register_sum_long_double(void);
  
