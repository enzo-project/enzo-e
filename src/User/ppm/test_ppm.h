const int num_problems = 4;
const int num_ghosts   = 3;

enum type_problem {  
  problem_ppm_unknown, 
  problem_ppm_image,  
  problem_ppm_implosion,  
  problem_ppm_implosion3
};

const char * problem_name[] = {
  "",
  "ppm-image",
  "ppm-implosion",
  "ppm-implosion3"
};

const int problem_size [] = {
  0,
  512,
  400,
  32 
};

const int problem_cycles [] = {
  0,
  10000,
  10000,
  10000 
};
