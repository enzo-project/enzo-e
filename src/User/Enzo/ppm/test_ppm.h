const int num_problems = 4;
const int num_ghosts   = 3;

enum type_problem {  
  problem_unknown, 
  problem_image,  
  problem_implosion,  
  problem_implosion3
};

const char * problem_name[] = {
  "",
  "image",
  "implosion",
  "implosion3"
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
