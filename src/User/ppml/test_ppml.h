const int num_ghosts   = 3;

enum problem_ppml_type {
  problem_ppml_unknown,
  problem_ppml_blast,
  problem_ppml_implosion3,
  num_problems
};

const char * problem_name[] = {
  "",
  "ppml-blast",
  "ppml-implosion3"
};

const int problem_size [] = {
  0,
  32,
  32
};

const int problem_cycles [] = {
  0,
  10000,
  10000
};
