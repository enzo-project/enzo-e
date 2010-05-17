const int num_problems = 3;
const int num_ghosts   = 3;

enum type_problem {
  problem_unknown,
  problem_blast,
  problem_implosion3
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
