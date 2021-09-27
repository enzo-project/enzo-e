//extern CkReduction::reducerType r_method_turbulence_mhd_ig_type;
extern CkReduction::reducerType r_method_turbulence_mhd_it_type;
extern CkReduction::reducerType r_method_turbulence_type;

//extern CkReductionMsg * r_method_turbulence_mhd_ig(int n, CkReductionMsg ** msgs);
extern CkReductionMsg * r_method_turbulence_mhd_it(int n, CkReductionMsg ** msgs);
extern CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs);

//extern void register_method_turbulence_mhd_ig(void);
extern void register_method_turbulence_mhd_it(void);
extern void register_method_turbulence(void);

