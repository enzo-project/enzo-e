namespace test {
  const int max_level = 3;
  struct face_test_type {
    int child[3];
    int size[3];
    int centered[3];
    double cell_width[3];
    double time_step_1;
    double time_step_2;
  };

  const int num_face_tests = 1;
  face_test_type face_test[] =
    {
      {
       .child = { 0, 0, 0}, .size   = {10, 6, 8 },
       .centered = { 0, 0, 0}, .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125, .time_step_2 = 0.125
      },
    };
};
