namespace test {
  const int max_level = 3;
  struct face_test_type {
    int test;
    int rank;
    int face[3];
    int normal[3];
    int level_1;
    int level_2;
    int child[3];
    int size[3];
    int centered[3];
    double cell_width[3];
    double time_step_1;
    double time_step_2;
  };

  
  const int num_face_tests = 36;
  face_test_type face_test[] =
    {
      // Non-boundary
      // level 0 x-axis
      {
        .test = 0,
        .rank = 3,
        .face   = {-1, 0, 0},
        .normal = {-1, 0, 0},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = {10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // level 0 y-axis
      {
        .test = 1,
        .rank = 3,
        .face   = {0, +1, 0},
        .normal = {0, +1, 0},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // level 0 z-axis
      {
        .test = 2,
        .rank = 3,
        .face   = {0, 0, +1},
        .normal = {0, 0, +1},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root same parent same level x-axis
      {
        .test = 3,
        .rank = 3,
        .face   = {-1, 0, 0},
        .normal = {+1, 0, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root same parent same level y-axis
      {
        .test = 4,
        .rank = 3,
        .face   = {0, -1, 0},
        .normal = {0, -1, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root same parent same level z-axis
      {
        .test = 5,
        .rank = 3,
        .face   = {0, 0, +1},
        .normal = {0, 0, +1},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root different parent same level x-axis
      {
        .test = 6,
        .rank = 3,
        .face   = {+1, 0, 0},
        .normal = {+1, 0, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root different parent same level y-axis
      {
        .test = 7,
        .rank = 3,
        .face   = {0, -1, 0},
        .normal = {0, -1, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root different parent same level z-axis
      {
        .test = 8,
        .rank = 3,
        .face   = {0, 0, +1},
        .normal = {0, 0, +1},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // different root same level x-axis
      {
        .test = 9,
        .rank = 3,
        .face   = {+1, 0, 0},
        .normal = {-1, 0, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // different root same level y-axis
      {
        .test = 10,
        .rank = 3,
        .face   = {0, -1, 0},
        .normal = {0, +1, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // different root same level z-axis
      {
        .test = 11,
        .rank = 3,
        .face   = {0, 0, +1},
        .normal = {0, 0, -1},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // coarser level x-axis
      {
        .test = 12,
        .rank = 3,
        .face   = {+1, 0, 0},
        .normal = {-1, 0, 0},
        .level_1 = 1,
        .level_2 = 2,
        .child = {1, 0, 1},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // finer level x-axis
      {
        .test = 13,
        .rank = 3,
        .face   = {-1, 0, 0},
        .normal = {-1, 0, 0},
        .level_1 = 2,
        .level_2 = 1,
        .child = {1, 0, 1},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // coarser level y-axis
      {
        .test = 14,
        .rank = 3,
        .face   = {0, -1, 0},
        .normal = {0, +1, 0},
        .level_1 = 1,
        .level_2 = 2,
        .child = {1, 0, 0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // finer level y-axis
      {
        .test = 15,
        .rank = 3,
        .face   = {0, +1, 0},
        .normal = {0, +1, 0},
        .level_1 = 2,
        .level_2 = 1,
        .child = {1, 0, 0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // coarser level z-axis
      {
        .test = 16,
        .rank = 3,
        .face   = {0, 0, +1},
        .normal = {0, 0, -1},
        .level_1 = 1,
        .level_2 = 2,
        .child = {1, 0, 1},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // finer level z-axis
      {
        .test = 17,
        .rank = 3,
        .face   = {0, 0, -1},
        .normal = {0, 0, -1},
        .level_1 = 2,
        .level_2 = 1,
        .child = {1, 0, 1},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // P3 B.C.'S
      // p3 level 0 x-axis p3 boundary
      {
        .test = 18,
        .rank = 3,
        .face   = {+1, 0, 0},
        .normal = {+1, 0, 0},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // p3 level 0 y-axis
      {
        .test = 19,
        .rank = 3,
        .face   = {0, -1, 0},
        .normal = {0, -1, 0},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // p3 level 0 z-axis
      {
        .test = 20,
        .rank = 3,
        .face   = {0, 0, -1},
        .normal = {0, 0, -1},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // p3 same level x-axis
      {
        .test = 21,
        .rank = 3,
        .face   = {-1, 0, 0},
        .normal = {-1, 0, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // p3 same level y-axis
      {
        .test = 22,
        .rank = 3,
        .face   = {0, +1, 0},
        .normal = {0, +1, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // p3 same level z-axis
      {
        .test = 23,
        .rank = 3,
        .face   = {0, 0, -1},
        .normal = {0, 0, -1},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 8 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      //======================================================================
      // 2D tests
      //======================================================================
      // level 0 x-axis
      {
        .test = 24,
        .rank = 2,
        .face   = {-1, 0},
        .normal = {-1, 0},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // level 0 y-axis
      {
        .test = 25,
        .rank = 2,
        .face   = {0, +1},
        .normal = {0, +1},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root same parent same level x-axis
      {
        .test = 26,
        .rank = 2,
        .face   = {-1, 0},
        .normal = {+1, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root same parent same level y-axis
      {
        .test = 27,
        .rank = 2,
        .face   = {0, -1},
        .normal = {0, -1},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root different parent same level x-axis
      {
        .test = 28,
        .rank = 2,
        .face   = {+1, 0},
        .normal = {+1, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // same root different parent same level y-axis
      {
        .test = 29,
        .rank = 2,
        .face   = {0, -1},
        .normal = {0, -1},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // different root same level x-axis
      {
        .test = 30,
        .rank = 2,
        .face   = {+1, 0},
        .normal = {-1, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // different root same level y-axis
      {
        .test = 31,
        .rank = 2,
        .face   = {0, -1},
        .normal = {0, +1},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // P3 b.c.'s
      // p3 level 0 x-axis p3 boundary
      {
        .test = 32,
        .rank = 2,
        .face   = {+1, 0},
        .normal = {+1, 0},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // p3 level 0 y-axis
      {
        .test = 33,
        .rank = 2,
        .face   = {0, -1},
        .normal = {0, -1},
        .level_1 = 0,
        .level_2 = 0,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // p3 same level x-axis
      {
        .test = 34,
        .rank = 2,
        .face   = {-1, 0},
        .normal = {-1, 0},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      },
      // p3 same level y-axis
      {
        .test = 35,
        .rank = 2,
        .face   = {0, +1},
        .normal = {0, +1},
        .level_1 = 2,
        .level_2 = 2,
        .child = {0,0,0},
        .size = { 10, 6, 1 },
        .centered = { 0, 0, 0},
        .cell_width = { 0.25, 0.5, 1.0},
        .time_step_1 = 0.125,
        .time_step_2 = 0.125
      }
    };
};
