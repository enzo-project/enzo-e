#include <sstream>

#include "cello.hpp"
#include "enzo.hpp"

StarFind::StarFind ()
{
#ifdef CONFIG_USE_TORCH
  // load S1 checkpoint
  const std::string stage1_filepath = "/home1/07320/whick002/StarNetRuntime/model_checkpoints/smalldense.jtpt";
  stage1_checkpoint_ = this->load_checkpoint(stage1_filepath); 

  // load S2 checkpoint
  const std::string stage2_filepath = "/home1/07320/whick002/StarNetRuntime/model_checkpoints/incepunet.jtpt";
  stage2_checkpoint_ = this->load_checkpoint(stage2_filepath);
#endif
}

//---------------------------

FBNet::FBNet ()
{
  const EnzoConfig * enzo_config = enzo::config();

  // read all input file
  // TODO: make input filenames parameters to EnzoMethodInference

  std::string mass_CDF_file = "/home1/07320/whick002/enzo-e_inference/input/FBNet_inputs/CDF_mass.txt";
  std::vector<std::vector<double>*> mass_vars = {&mass_CDF_, &mass_CDF_bins_};
  read_file(mass_CDF_file, &mass_vars);

  std::string Nstar_CDF_file = "/home1/07320/whick002/enzo-e_inference/input/FBNet_inputs/CDF_Nstar.txt";
  std::vector<std::vector<double>*> Nstar_vars = {&Nstar_CDF_, &Nstar_CDF_bins_};
  read_file(Nstar_CDF_file, &Nstar_vars );

  std::string creationtime_CDF_file = "/home1/07320/whick002/enzo-e_inference/input/FBNet_inputs/CDF_creationtime.txt";
  std::vector<std::vector<double>*> creationtime_vars = {&creationtime_CDF_, &creationtime_CDF_bins_};
  read_file(creationtime_CDF_file, &creationtime_vars );


  // Only load in one regression weights file. Which one to load is an input parameter
  //       numbers 7-20 correspond to "model_time", which is an input parameter. Seems like
  //       IMF_sampler.py is hard-coded to use regression_weights_12.txt though.
  //NOTE: See Tables 3 and 4 of https://iopscience.iop.org/article/10.3847/1538-4357/ac6c87/pdf
  
  int model_time = 12; // length of time Pop III FB is active for??
  std::string weights_file = "/home1/07320/whick002/enzo-e_inference/input/FBNet_inputs/regression_weights_" + std::to_string(model_time) + ".txt";
  std::vector<std::vector<double>*> regression_vars = {&M0_, &M1_, &M2_, &M3_};
  read_file(weights_file, &regression_vars);
}

//---------------------------

void FBNet::read_file(std::string file, std::vector< std::vector<double> * > * vars) throw()
{
  // reads in text file and stores each row in a separate array
  std::ifstream inFile;
  inFile.open(file);

  ASSERT1("FBNet::read_file()", "file %s failed to open!", file.c_str(), inFile.is_open());

  int line_counter = 0;
  for (std::string line; std::getline(inFile, line, '\n'); ) {
    std::stringstream s(line);
    for (std::string val; std::getline(s, val, ' '); ) {
      std::cout << file << ' ' << line_counter << std::endl;
      std::cout << *((*vars)[line_counter]) << std::endl;
      (*((*vars)[line_counter])).push_back(std::stod(val));
    }
    line_counter++;
  }
}

//---------------------------

int FBNet::get_Nstars() throw() 
{
  srand(time(NULL));
  double random = (double) rand()/RAND_MAX;
  
  int Nstars = 0;
  // find into which bin the random number falls
  for (int i=0; i<Nstar_CDF_.size()-1; i++) {
    if ((Nstar_CDF_[i] <= random) && (random < Nstar_CDF_[i+1])) {
      Nstars = (int) Nstar_CDF_bins_[i];
      break;
    }
  }

  return Nstars;
}

double FBNet::get_mass() throw() 
{
  srand(time(NULL));
  double random = (double) rand()/RAND_MAX;
  
  double Mstar = 0;
  // find into which bin the random number falls
  for (int i=0; i<mass_CDF_.size()-1; i++) {
    if ((mass_CDF_[i] <= random) && (random < mass_CDF_[i+1])) {
      Mstar = (double) mass_CDF_bins_[i];
      break;
    }
  }

  return Mstar;
}

double FBNet::get_creationtime() throw() 
{
  srand(time(NULL));
  double random = (double) rand()/RAND_MAX;
  
  int creationtime = 0;
  // find into which bin the random number falls
  for (int i=0; i<creationtime_CDF_.size()-1; i++) {
    if ((creationtime_CDF_[i] <= random) && (random < creationtime_CDF_[i+1])) {
      creationtime = (int) creationtime_CDF_bins_[i];
      break;
    }
  }

  return creationtime;
}

//---------------------------

double FBNet::get_metal_yield(double mass) throw() 
{
  if ((11.0 <= mass) && (mass < 20.0)) {
    return metal_yield_SNe(mass);
  } 
  else if ((20.0 <= mass) && (mass < 40.0)) {
    return metal_yield_HNe(mass);
  }
  else if ((140.0 <= mass) && (mass < 260.0)) {
    return metal_yield_PISNe(mass);
  }
  else {
    return 0.0;
  }
}

double FBNet::metal_yield_SNe(double mass) throw() 
{
  return 0.1077 + 0.3383 * (mass-11.0);
}

double FBNet::metal_yield_HNe(double mass) throw() 
{
  std::vector<double> HypernovaMetals = {3.36 , 3.53 , 5.48 , 7.03 , 8.59}; // Msun
  std::vector<double> HypernovaMass   = {19.99, 25.0 , 30.0 , 35.0 , 40.01};

  double yield = 0.0;
  for (int i=0; i<HypernovaMetals.size()-1; i++) {
    if ((HypernovaMass[i] <= mass) && (mass < HypernovaMass[i+1])) {
      double frac = (mass - HypernovaMass[i]) / (HypernovaMass[i+1]-HypernovaMass[i]);
      yield = HypernovaMetals[i] + frac*(HypernovaMetals[i+1]-HypernovaMetals[i]);
    }
  }
  return yield; 
}

double FBNet::metal_yield_PISNe(double mass) throw() 
{
  double He_core_mass = (13.0/24.0) * (mass-20.0);
  return 5.0 + 1.304*(He_core_mass - 64.0);
}

//-----------------------

void FBNet::get_bincounts(std::vector<int> * bin_counts,
                          std::vector<double> * masses   , std::vector<double> * creationtimes,
                          std::vector<double> * mass_bins, std::vector<double> * time_bins) throw() 
{    
  
}
