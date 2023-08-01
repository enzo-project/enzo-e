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
  read_file(mass_CDF_file, mass_vars);

  std::string Nstar_CDF_file = "/home1/07320/whick002/enzo-e_inference/input/FBNet_inputs/CDF_Nstar.txt";
  std::vector<std::vector<double>*> Nstar_vars = {&Nstar_CDF_, &Nstar_CDF_bins_};
  read_file(Nstar_CDF_file, Nstar_vars );

  std::string creationtime_CDF_file = "/home1/07320/whick002/enzo-e_inference/input/FBNet_inputs/CDF_creationtime.txt";
  std::vector<std::vector<double>*> creationtime_vars = {&creationtime_CDF_, &creationtime_CDF_bins_};
  read_file(creationtime_CDF_file, creationtime_vars );


  // Only load in one regression weights file. Which one to load is an input parameter
  //       numbers 7-20 correspond to "model_time", which is an input parameter. Seems like
  //       IMF_sampler.py is hard-coded to use regression_weights_12.txt though.
  //NOTE: See Tables 3 and 4 of https://iopscience.iop.org/article/10.3847/1538-4357/ac6c87/pdf
  
  model_time_ = 12; // length of time in Myr that we're modelling
  model_dt_ = 6; // determines number of time bins per call to model (Nbins = model_time/model_dt)

  Nstar_mean_ = 1.311; // log-mean number of stars per population (used for tokenization)
  Nstar_std_ = 0.352;
  
  std::string weights_file = "/home1/07320/whick002/enzo-e_inference/input/FBNet_inputs/regression_weights_" + std::to_string(model_time_) + ".txt";
  std::vector<std::vector<double>*> regression_vars = {&M0_, &M1_, &M2_, &M3_};
  read_file(weights_file, regression_vars);

}

//---------------------------

void FBNet::read_file(std::string file, std::vector< std::vector<double> * > vars) throw()
{
  // reads in text file and stores each row in a separate array
  std::ifstream inFile;
  inFile.open(file);

  ASSERT1("FBNet::read_file()", "file %s failed to open!", file.c_str(), inFile.is_open());

  int line_counter = 0;
  for (std::string line; std::getline(inFile, line, '\n'); ) {
    std::stringstream s(line);
    for (std::string val; std::getline(s, val, ' '); ) {
      (*(vars[line_counter])).push_back(std::stod(val));
    }
    line_counter++;
  }
}

//---------------------------

int FBNet::get_Nstars() throw() 
{
  double random = (double) rand()/RAND_MAX;
  
  int Nstars = 0;
  // find into which bin the random number falls
  for (int i=0; i<Nstar_CDF_.size()-1; i++) {
    if ((Nstar_CDF_[i] <= random) && (random < Nstar_CDF_[i+1])) {
      Nstars = (int) Nstar_CDF_bins_[i];
      break;
    }
  }

  return std::max(Nstars, 1); // always return at least 1 star
}

double FBNet::get_mass() throw() 
{
  double random = (double) rand()/RAND_MAX;
  
  double Mstar = 0;
  // find into which bin the random number falls
  for (int i=0; i<mass_CDF_.size()-1; i++) {
    if ((mass_CDF_[i] <= random) && (random < mass_CDF_[i+1])) {
      Mstar = (double) mass_CDF_bins_[i];
      break;
    }
  }

  return std::max(Mstar, 1.0); // always return at least 1 Msun
}

double FBNet::get_creationtime() throw() 
{
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

int FBNet::get_binindex(double val, std::vector<double> bins) throw()
{
  for (int i=0; i < bins.size()-1; i++) {
    if ( (bins[i] <= val) && (val < bins[i+1]) ) {
      return i;
    }
  }

  return INT_MAX; // return huge number if val doesn't fall within any of the bins
}

//--------------------

double FBNet::get_radius(std::vector<double> masses, std::vector<double> creationtimes) throw()
{
  // Start by tokenizing sample

  // initialize mass bins (TODO: link to parameter?)
  std::vector<double> massbins = {1.0, 11.0, 20.0, 40.0, 100.0, 140.0, 200.0, 260.0, 300.0};
  std::vector<double> mass_bincounts;  
  mass_bincounts.resize(massbins.size()-1);
  std::fill(mass_bincounts.begin(), mass_bincounts.end(), 0);

  // initialize time bins
  std::vector<double> timebins;
  for (int i=0; i < std::floor(model_time_/model_dt_); i++) {
    timebins.push_back(i * (double) model_dt_);  
  }
  std::vector<double> time_bincounts;
  time_bincounts.resize(timebins.size()-1);
  std::fill(time_bincounts.begin(), time_bincounts.end(), 0);
  
  // get bincounts
  int Nstars = masses.size();
  for (int i=0; i < Nstars; i++) {
    mass_bincounts[get_binindex( masses[i], massbins )] += 1;
    time_bincounts[get_binindex( creationtimes[i], timebins )] += 1;
  }
  int Nbins_mass = mass_bincounts.size();
  int Nbins_time = time_bincounts.size();
  int Nbins = Nbins_mass + Nbins_time;

  // determine if we're in a sparsely unpopulated, averagely populated, or densely populated system
  double N1 = std::pow(10, Nstar_mean_ - Nstar_std_);
  double N2 = std::pow(10, Nstar_mean_);
  double N3 = std::pow(10, Nstar_mean_ + Nstar_std_);

  std::vector<double> * M = NULL;
  if (Nstars < N1) {
    M = &M0_;
  }
  else if ( (N1 <= Nstars) && (Nstars < N2) ) {
    M = &M1_;
  }
  else if ( (N2 <= Nstars) && (Nstars < N3) ) {
    M = &M2_;
  }
  else { // N3 <= Nstars
    M = &M3_;
  }

  ASSERT("FBNet::get_radius()", "Nbins+1 != M.size()", 
          Nbins+1 == (*M).size()); // plus 1 to include bias feature

  ASSERT("FBNet::get_radius()", "&M is NULL pointer!",
          M != NULL);


  bool include_timebins = true; // TODO: Make this a parameter!

  double radius = 1 * (*M)[0]; // index 0 in M reserved for bias feature
  for (int i=0; i < Nbins_mass; i++) {
    radius += mass_bincounts[i] * (*M)[i+1];
  }
  if (include_timebins) {
    for (int i=0; i < Nbins_time; i++) {
      radius += time_bincounts[i] * (*M)[Nbins_mass+1 + i];
    }
  }
  // NOTE: this matrix multiplication gives us log(r) in kpc
  radius = std::pow(10, std::max(radius, std::log10(0.25)));

  return radius;
}

//---------------------------------------------------------

void FBNet::update_mesh(EnzoBlock * enzo_block, EnzoObjectFeedbackSphere sphere) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();
  double lunit = enzo_units->length();
  double munit = enzo_units->mass();

  Field field = enzo_block->data()->field();

  int mx,my,mz;  
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);

  double xm,ym,zm;
  double xp,yp,zp;
  double hx,hy,hz;
  enzo_block->data()->lower(&xm,&ym,&zm);
  enzo_block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  enzo_float * density = (enzo_float *) field.values("density");
  enzo_float * metal_density = (enzo_float *) field.values("metal_density");
  enzo_float * PopIII_metal_density  = (enzo_float *) field.values("PopIII_metal_density");
  enzo_float * PopIII_SNe_metal_density = (enzo_float *) field.values("PopIII_SNe_metal_density");
  enzo_float * PopIII_HNe_metal_density = (enzo_float *) field.values("PopIII_HNe_metal_density");
  enzo_float * PopIII_PISNe_metal_density = (enzo_float *) field.values("PopIII_PISNe_metal_density");

  double metal_mass_tot = sphere.metal_mass_SNe() + sphere.metal_mass_HNe() + sphere.metal_mass_PISNe();
  if (metal_mass_tot == 0.0) { 
    // do nothing if no yield
    return;
  }
  double x = sphere.pos(0), y = sphere.pos(1), z = sphere.pos(2);

  // get 3D grid index for sphere center - account for ghost zones!!
  int ix = (int) std::floor((x - xm) / hx) + gx;
  int iy = (int) std::floor((y - ym) / hy) + gy;
  int iz = (int) std::floor((z - zm) / hz) + gz;

  double r_code = sphere.r();
  double inv_vol = 3 / (4*cello::pi * r_code*r_code*r_code);

  // compute average density for deposited metals
  double drho_SNe   = sphere.metal_mass_SNe()  * inv_vol;
  double drho_HNe   = sphere.metal_mass_HNe()  * inv_vol;
  double drho_PISNe = sphere.metal_mass_PISNe()* inv_vol;
  double drho = drho_SNe + drho_HNe + drho_PISNe;


  #ifdef DEBUG_METHOD_FBNET
    CkPrintf("[%d] FBNet::update_mesh -- rho = %1.2e, rho_SNe = %1.2e, rho_HNe = %1.2e, rho_PISNe = %1.2e\n", CkMyPe(), drho, drho_SNe, drho_HNe, drho_PISNe);
  #endif

  int r_hx = std::ceil(r_code/hx), r_hy = std::ceil(r_code/hy), r_hz = std::ceil(r_code/hz);
  // CiC deposit with cloud radius of r_code
  // Feedback spheres may intersect neighboring blocks in such a way that
  // the volume of the sphere will extend past the available ghost zones.
  // Account for this by looping through ALL identified feedback spheres in the domain,
  // and depositing portions of external spheres that intersect the block.
  // Since spheres have uniform density, this is as simple as adding a check in the 
  // CiC routine that the current cell being indexed actually lies within the block.

  for (int iz_ = iz-r_hz; iz_ <= iz+r_hz; iz_++) {
    // if out of bounds, go to next iteration
    if ((iz_ < 0) || (mz <= iz_)) continue;

    for (int iy_ = iy-r_hy; iy_ <= iy+r_hy; iy_++) {
      if ((iy_ < 0) || (my <= iy_)) continue;

      for (int ix_ = ix-r_hx; ix_ <= ix+r_hx; ix_++) {
        if ((ix_ < 0) || (mx <= ix_)) continue;

        int i_ = INDEX(ix_,iy_,iz_,mx,my);

        // if cell is within the deposition radius
        bool contained = ( (iz-iz_)*(iz-iz_) + (iy-iy_)*(iy-iy_) + (ix-ix_)*(ix-ix_) <= r_hx*r_hx );
        density[i_]                    += contained*drho;
        metal_density[i_]              += contained*drho;
        PopIII_metal_density[i_]       += contained*drho;
        PopIII_SNe_metal_density[i_]   += contained*drho_SNe;
        PopIII_HNe_metal_density[i_]   += contained*drho_HNe;
        PopIII_PISNe_metal_density[i_] += contained*drho_PISNe;
  
        //if (contained) {    
        //  CkPrintf("FBNet::update_mesh -- metal density = %1.2e; radius = %1.2e kpc\n", drho, r_code * lunit/enzo_constants::kpc_cm);
        //}


      }
    }
  }

}
