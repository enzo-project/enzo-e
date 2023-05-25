class StarFind {

public:
  StarFind();

#ifdef CONFIG_USE_TORCH
  // store deserialized checkpoint files as attributes
  torch::jit::script::Module stage1() const throw() {return stage1_checkpoint_;}

  torch::jit::script::Module stage2() const throw() {return stage2_checkpoint_;}

private:
  torch::jit::script::Module load_checkpoint(std::string file) throw() {
    return torch::jit::load(file);
  }

  torch::jit::script::Module stage1_checkpoint_, stage2_checkpoint_;

#endif
};


class FBNet {

public:
  FBNet();

  double mass_CDF     (int i) const throw() {return mass_CDF_[i];}
  double mass_CDF_bins(int i) const throw() {return mass_CDF_bins_[i];}

  double Nstar_CDF     (int i) const throw() {return Nstar_CDF_[i];}
  double Nstar_CDF_bins(int i) const throw() {return Nstar_CDF_bins_[i];}

  double creationtime_CDF     (int i) const throw() {return creationtime_CDF_[i];}
  double creationtime_CDF_bins(int i) const throw() {return creationtime_CDF_bins_[i];}

  double M0(int i) const throw() {return M0_[i];}
  double M1(int i) const throw() {return M1_[i];}
  double M2(int i) const throw() {return M2_[i];}
  double M3(int i) const throw() {return M3_[i];}


  int get_Nstars() throw();

  double get_mass() throw(); 

  double get_creationtime() throw();

  double get_metal_yield(double mass) throw();

  double metal_yield_SNe(double mass) throw();

  double metal_yield_HNe(double mass) throw();

  double metal_yield_PISNe(double mass) throw();

  void get_bincounts(std::vector<int> * bin_counts,
                       std::vector<double> * masses   , std::vector<double> * creationtimes,
                       std::vector<double> * mass_bins, std::vector<double> * time_bins) throw(); 

private:
  void read_file(std::string file, std::vector<std::vector<double>*> * vars) throw();

  void tokenize_sample() throw();

  std::vector<double> mass_CDF_, mass_CDF_bins_; 
  std::vector<double> Nstar_CDF_, Nstar_CDF_bins_;
  std::vector<double> creationtime_CDF_, creationtime_CDF_bins_;

  std::vector<double> M0_, M1_, M2_, M3_; 

};
