extern CkReductionMsg * r_reduce_performance(int n, CkReductionMsg ** msgs);
extern CkReduction::reducerType r_reduce_performance_type;
extern void register_reduce_performance(void);

extern CkReductionMsg * sum_long_double(int n, CkReductionMsg ** msgs);
extern CkReduction::reducerType sum_long_double_type;
extern void register_sum_long_double(void);

extern CkReductionMsg * sum_long_double_2(int n, CkReductionMsg ** msgs);
extern CkReduction::reducerType sum_long_double_2_type;
extern void register_sum_long_double_2(void);

extern CkReductionMsg * sum_long_double_3(int n, CkReductionMsg ** msgs);
extern CkReduction::reducerType sum_long_double_3_type;
extern void register_sum_long_double_3(void);

extern CkReductionMsg * sum_long_double_n(int n, CkReductionMsg ** msgs);
extern CkReduction::reducerType sum_long_double_n_type;
extern void register_sum_long_double_n(void);

extern CkReductionMsg * r_reduce_method_debug(int n, CkReductionMsg ** msgs);
extern CkReduction::reducerType r_reduce_method_debug_type;
extern void register_reduce_method_debug(void);

