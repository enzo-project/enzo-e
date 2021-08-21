extern CkReduction::reducerType r_method_mhd_turbulence_it_type;
extern CkReduction::reducerType r_method_turbulence_type;

extern CkReductionMsg * r_method_mhd_turbulence_it(int n, CkReductionMsg ** msgs);
extern CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs);

extern void register_method_mhd_turbulence_it(void);
extern void register_method_turbulence(void);

