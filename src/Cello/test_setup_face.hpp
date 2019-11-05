  const int max_level = 3;
  struct test_type {
    int rank;  int array[3];   bool periodic[3];
    int levels_1;  int array_1[3];  int tree_1[max_level][3];
    int levels_2;  int array_2[3];  int tree_2[max_level][3];
    float size[3];  int normal[3];  int NX,NY,NZ;
    int axis; int face;
    bool l_degenerate;
  };

  Face * prev_A = nullptr;
  Face * prev_B = nullptr;
  
  const int num_tests = 49;
  test_type test[] =
    {
      // Non-boundary
      // test 0: level 0 x-axis
      {
        .rank = 3, .array = {4,5,1}, .periodic = {false,false,false},
        .levels_1 = 0, .array_1 = {3,4,0}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {2,4,0}, .tree_2 = { },
        .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = -1,
        .l_degenerate = false
      },
      // test 1: level 0 y-axis
      {
        .rank = 3, .array = {3,5,4}, .periodic = {false,false,false},
        .levels_1 = 0, .array_1 = {0,0,1}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {0,1,1}, .tree_2 = { },
        .size = {1.0, 0.0, 1.0}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = +1,
        .l_degenerate = false
      },
      // test 2: level 0 z-axis
      {
        .rank = 3, .array = {4,5,7}, .periodic = {false,false,false},
        .levels_1 = 0, .array_1 = {2,1,5}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {2,1,6}, .tree_2 = { },
        .size = {1.0, 1.0, 0.0}, .normal = {0,0,1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = +1,
        .l_degenerate = false
      },
      // test 3: same root same parent same level x-axis
      {
        .rank = 3, .array = {4,5,3}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {3,1,2}, .tree_1 = {{ 0, 1, 0}, {1, 0, 1}},
        .levels_2 = 2, .array_2 = {3,1,2}, .tree_2 = {{ 0, 1, 0}, {0, 0, 1}},
        .size = {0.0, 1.0, 1.0}, .normal = {1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = -1,
        .l_degenerate = false
      },
      // test 4: same root same parent same level y-axis
      {
        .rank = 3, .array = {4,5,4}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {3,1,2}, .tree_1 = {{ 1, 1, 0}, {0, 1, 1}},
        .levels_2 = 2, .array_2 = {3,1,2}, .tree_2 = {{ 1, 1, 0}, {0, 0, 1}},
        .size = {1.0, 0.0, 1.0}, .normal = {0,-1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 5: same root same parent same level z-axis
      {
        .rank = 3, .array = {4,5,5}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {3,1,2}, .tree_1 = {{ 1, 0, 0}, {1, 1, 0}},
        .levels_2 = 2, .array_2 = {3,1,2}, .tree_2 = {{ 1, 0, 0}, {1, 1, 1}},
        .size = {1.0, 1.0, 0.0}, .normal = {0,0,1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = +1,
        .l_degenerate = false
      },
      // test 6: same root different parent same level x-axis
      {
        .rank = 3, .array = {3,1,3}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 1}, {1, 1, 0}},
        .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {0, 1, 0}},
        .size = {0.0, 1.0, 1.0}, .normal = {1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 7: same root different parent same level y-axis
      {
        .rank = 3, .array = {3,1,3}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 1, 1, 1}, {1, 0, 0}},
        .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 0, 1}, {1, 1, 0}},
        .size = {1.0, 0.0, 1.0}, .normal = {0,-1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 8: same root different parent same level z-axis
      {
        .rank = 3, .array = {2,2,4}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 1, 1, 0}, {1, 0, 1}},
        .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 0}},
        .size = {1.0, 1.0, 0.0}, .normal = {0,0,1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = +1,
        .l_degenerate = false
      },
      // test 9: different root same level x-axis
      {
        .rank = 3, .array = {2,5,3}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 1}, {0, 0, 1}},
        .levels_2 = 2, .array_2 = {0,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
        .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 10: different root same level y-axis
      {
        .rank = 3, .array = {3,5,7}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 0}, {1, 1, 0}},
        .levels_2 = 2, .array_2 = {1,1,2}, .tree_2 = {{ 0, 0, 0}, {1, 0, 0}},
        .size = {1.0, 0.0, 1.0}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 11: different root same level z-axis
      {
        .rank = 3, .array = {3,3,5}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {1,0,3}, .tree_1 = {{ 1, 1, 0}, {1, 0, 0}},
        .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
        .size = {1.0, 1.0, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = +1,
        .l_degenerate = false
      },
      // test 12: coarser level x-axis
      {
        .rank = 3, .array = {2,5,3}, .periodic = {false,false,false},
        .levels_1 = 1, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 1}},
        .levels_2 = 2, .array_2 = {0,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
        .size = {0.0, 0.5, 0.5}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 13: finer level x-axis
      {
        .rank = 3, .array = {2,5,3}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {0,0,2}, .tree_1 = {{ 1, 1, 1}, {1, 0, 1}},
        .levels_2 = 1, .array_2 = {1,0,2}, .tree_2 = {{ 0, 1, 1}},
        .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 14: coarser level y-axis
      {
        .rank = 3, .array = {3,5,7}, .periodic = {false,false,false},
        .levels_1 = 1, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 0}},
        .levels_2 = 2, .array_2 = {1,1,2}, .tree_2 = {{ 0, 0, 0}, {1, 0, 0}},
        .size = {0.5, 0.0, 0.5}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = -1,
        .l_degenerate = false
      },
      // test 15: finer level y-axis
      {
        .rank = 3, .array = {3,5,7}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {1,1,2}, .tree_1 = {{ 0, 0, 0}, {1, 0, 0}},
        .levels_2 = 1, .array_2 = {1,0,2}, .tree_2 = {{ 0, 1, 0}},
        .size = {1.0, 0.0, 1.0}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 16: coarser level z-axis
      {
        .rank = 3, .array = {3,3,5}, .periodic = {false,false,false},
        .levels_1 = 1, .array_1 = {1,0,3}, .tree_1 = {{ 1, 1, 0}},
        .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
        .size = {0.5, 0.5, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = +1,
        .l_degenerate = false
      },
      // test 17: finer level z-axis
      {
        .rank = 3, .array = {3,3,5}, .periodic = {false,false,false},
        .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 1, 1, 1}, {1, 0, 1}},
        .levels_2 = 1, .array_2 = {1,0,3}, .tree_2 = {{ 1, 1, 0}},
        .size = {1.0, 1.0, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = -1,
        .l_degenerate = false
      },
      // Periodic b.c.'s
      // test 18: periodic level 0 x-axis periodic boundary
      {
        .rank = 3, .array = {4,2,2}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {3,1,0}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
        .levels_2 = 0, .array_2 = {0,1,0}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
        .size = {0.0, 1.0, 1.0}, .normal = {1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 19: periodic level 0 y-axis
      {
        .rank = 3, .array = {3,5,4}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {0,0,1}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
        .levels_2 = 0, .array_2 = {0,4,1}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
        .size = {1.0, 0.0, 1.0}, .normal = {0,-1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 20: periodic level 0 z-axis
      {
        .rank = 3, .array = {4,5,7}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {2,1,0}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
        .levels_2 = 0, .array_2 = {2,1,6}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
        .size = {1.0, 1.0, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = -1,
        .l_degenerate = false
      },
      // test 21: periodic same level x-axis
      {
        .rank = 3, .array = {2,5,3}, .periodic = {true,true,true},
        .levels_1 = 2, .array_1 = {0,0,2}, .tree_1 = {{ 0, 1, 1}, {0, 0, 1}},
        .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
        .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = -1,
        .l_degenerate = false
      },
      // test 22: periodic same level y-axis
      {
        .rank = 3, .array = {3,5,7}, .periodic = {true,true,true},
        .levels_1 = 2, .array_1 = {1,4,2}, .tree_1 = {{ 0, 1, 0}, {1, 1, 0}},
        .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 0, 0, 0}, {1, 0, 0}},
        .size = {1.0, 0.0, 1.0}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = +1,
        .l_degenerate = false
      },
      // test 23: periodic same level z-axis
      {
        .rank = 3, .array = {3,3,5}, .periodic = {true,true,true},
        .levels_1 = 2, .array_1 = {1,0,0}, .tree_1 = {{ 1, 1, 0}, {1, 0, 0}},
        .levels_2 = 2, .array_2 = {1,0,4}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
        .size = {1.0, 1.0, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = -1,
        .l_degenerate = false
      },
      // test 24: non-periodic same level x-axis
      {
        .rank = 3, .array = {4,2,2}, .periodic = {false,false,false},
        .levels_1 = 0, .array_1 = {3,1,0}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
        .levels_2 = 0, .array_2 = {0,1,0}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
        .size = {-1.0, 1.0, 1.0}, .normal = {1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 25: non-periodic level 0 y-axis
      {
        .rank = 3, .array = {3,5,4}, .periodic = {false,false,false},
        .levels_1 = 0, .array_1 = {0,0,1}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
        .levels_2 = 0, .array_2 = {0,4,1}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
        .size = {1.0, -1.0, 1.0}, .normal = {0,-1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 26: periodic level 0 z-axis
      {
        .rank = 3, .array = {4,5,7}, .periodic = {false,false,false},
        .levels_1 = 0, .array_1 = {2,1,0}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
        .levels_2 = 0, .array_2 = {2,1,6}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
        .size = {1.0, 1.0, -1.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = -1,
        .l_degenerate = false
      },
      //======================================================================
      // 2D tests
      //======================================================================
      // test 27: level 0 x-axis
      {
        .rank = 2, .array = {4,5}, .periodic = {false,false},
        .levels_1 = 0, .array_1 = {3,4}, .tree_1 = { {-1,-1}, {-1,-1} },
        .levels_2 = 0, .array_2 = {2,4}, .tree_2 = { {-1,-1}, {-1,-1} },
        .size = {0.0, 1.0}, .normal = {-1,0}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 0, .face = -1,
        .l_degenerate = false
      },
      // test 28: level 0 y-axis
      {
        .rank = 2, .array = {3,5}, .periodic = {false,false},
        .levels_1 = 0, .array_1 = {0,0}, .tree_1 = {{-1,-1}, {-1,-1}},
        .levels_2 = 0, .array_2 = {0,1}, .tree_2 = {{-1,-1}, {-1,-1}},
        .size = {1.0, 0.0}, .normal = {0,1}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 1, .face = +1,
        .l_degenerate = false
      },
      // test 29: same root same parent same level x-axis
      {
        .rank = 2, .array = {4,5}, .periodic = {false,false},
        .levels_1 = 2, .array_1 = {3,1}, .tree_1 = {{ 0, 1}, {1, 0}},
        .levels_2 = 2, .array_2 = {3,1}, .tree_2 = {{ 0, 1}, {0, 0}},
        .size = {0.0, 1.0}, .normal = {1,0}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 0, .face = -1,
        .l_degenerate = false
      },
      // test 30: same root same parent same level y-axis
      {
        .rank = 2, .array = {4,5}, .periodic = {false,false},
        .levels_1 = 2, .array_1 = {3,1}, .tree_1 = {{ 1, 1}, {0, 1}},
        .levels_2 = 2, .array_2 = {3,1}, .tree_2 = {{ 1, 1}, {0, 0}},
        .size = {1.0, 0.0}, .normal = {0,-1}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 31 same root different parent same level x-axis
      {
        .rank = 2, .array = {3,1}, .periodic = {false,false},
        .levels_1 = 2, .array_1 = {1,0}, .tree_1 = {{ 0, 1}, {1, 1}},
        .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 1, 1}, {0, 1}},
        .size = {0.0, 1.0}, .normal = {1,0}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 32: same root different parent same level y-axis
      {
        .rank = 2, .array = {3,1}, .periodic = {false,false},
        .levels_1 = 2, .array_1 = {1,0}, .tree_1 = {{ 1, 1}, {1, 0}},
        .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 1, 0}, {1, 1}},
        .size = {1.0, 0.0}, .normal = {0,-1}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 33: different root same level x-axis
      {
        .rank = 2, .array = {2,5}, .periodic = {false,false},
        .levels_1 = 2, .array_1 = {1,0}, .tree_1 = {{ 0, 1}, {0, 0}},
        .levels_2 = 2, .array_2 = {0,0}, .tree_2 = {{ 1, 1}, {1, 0}},
        .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 34: different root same level y-axis
      {
        .rank = 2, .array = {3,5}, .periodic = {false,false},
        .levels_1 = 2, .array_1 = {1,0}, .tree_1 = {{ 0, 1}, {1, 1}},
        .levels_2 = 2, .array_2 = {1,1}, .tree_2 = {{ 0, 0}, {1, 0}},
        .size = {1.0, 0.0}, .normal = {0,1}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // Periodic b.c.'s
      // test 35: periodic level 0 x-axis periodic boundary
      {
        .rank = 2, .array = {4,2}, .periodic = {true,true},
        .levels_1 = 0, .array_1 = {3,1}, .tree_1 = {{-1,-1}, {-1,-1}},
        .levels_2 = 0, .array_2 = {0,1}, .tree_2 = {{-1,-1}, {-1,-1}},
        .size = {0.0, 1.0}, .normal = {1,0}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 36: periodic level 0 y-axis
      {
        .rank = 2, .array = {3,5}, .periodic = {true,true},
        .levels_1 = 0, .array_1 = {0,0}, .tree_1 = {{-1,-1}, {-1,-1}},
        .levels_2 = 0, .array_2 = {0,4}, .tree_2 = {{-1,-1}, {-1,-1}},
        .size = {1.0, 0.0}, .normal = {0,-1}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 37: periodic same level x-axis
      {
        .rank = 2, .array = {2,5}, .periodic = {true,true},
        .levels_1 = 2, .array_1 = {0,0}, .tree_1 = {{ 0, 1}, {0, 0}},
        .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 1, 1}, {1, 0}},
        .size = {0.0, 1.0}, .normal = {-1,0}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 0, .face = -1,
        .l_degenerate = false
      },
      // test 38: periodic same level y-axis
      {
        .rank = 2, .array = {3,5}, .periodic = {true,true},
        .levels_1 = 2, .array_1 = {1,4}, .tree_1 = {{ 0, 1}, {1, 1}},
        .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 0, 0}, {1, 0}},
        .size = {1.0, 0.0}, .normal = {0,1}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 1, .face = +1,
        .l_degenerate = false
      },
      // Non-periodic b.c.'s
      // test 39: periodic level 0 x-axis non-periodic boundary
      {
        .rank = 2, .array = {4,2}, .periodic = {false,false},
        .levels_1 = 0, .array_1 = {3,1}, .tree_1 = {{-1,-1}, {-1,-1}},
        .levels_2 = 0, .array_2 = {0,1}, .tree_2 = {{-1,-1}, {-1,-1}},
        .size = {-1.0, 1.0}, .normal = {1,0}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 0, .face = +1,
        .l_degenerate = false
      },
      // test 40: periodic level 0 y-axis non-periodic boundary
      {
        .rank = 2, .array = {3,5}, .periodic = {false,false},
        .levels_1 = 0, .array_1 = {0,0}, .tree_1 = {{-1,-1}, {-1,-1}},
        .levels_2 = 0, .array_2 = {0,4}, .tree_2 = {{-1,-1}, {-1,-1}},
        .size = {1.0,-1.0}, .normal = {0,-1}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 1, .face = -1,
        .l_degenerate = false
      },
      // test 41: periodic same level x-axis non-periodic boundary
      {
        .rank = 2, .array = {2,5}, .periodic = {false,false},
        .levels_1 = 2, .array_1 = {0,0}, .tree_1 = {{ 0, 1}, {0, 0}},
        .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 1, 1}, {1, 0}},
        .size = {-1.0, 1.0}, .normal = {-1,0}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 0, .face = -1,
        .l_degenerate = false
      },
      // test 42: periodic same level y-axis non-periodic boundary
      {
        .rank = 2, .array = {3,5}, .periodic = {false,false},
        .levels_1 = 2, .array_1 = {1,4}, .tree_1 = {{ 0, 1}, {1, 1}},
        .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 0, 0}, {1, 0}},
        .size = {1.0, -1.0}, .normal = {0,1}, .NX = 10, .NY = 6, .NZ = 0,
        .axis = 1, .face = +1,
        .l_degenerate = false
      },
      // test 43: same block, periodic boundaries
      {
        .rank = 3, .array = {1,5,3}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {0,2,1}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {0,2,1}, .tree_2 = { },
        .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = -1,
        .l_degenerate = true
      },
      // test 44: same block, periodic boundaries
      {
        .rank = 3, .array = {5,1,3}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {3,0,2}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {3,0,2}, .tree_2 = { },
        .size = {1.0, 0.0, 1.0}, .normal = {0,+1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = +1,
        .l_degenerate = true
      },
      // test 45: same block, periodic boundaries
      {
        .rank = 3, .array = {3,5,1}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {2,4,0}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {2,4,0}, .tree_2 = { },
        .size = {1.0, 1.0, 0.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = -1,
        .l_degenerate = true
      },
      // test 46: same block, periodic boundaries
      {
        .rank = 3, .array = {1,5,3}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {0,2,1}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {0,2,1}, .tree_2 = { },
        .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 0, .face = -1,
        .l_degenerate = true
      },
      // test 47: same block, periodic boundaries
      {
        .rank = 3, .array = {5,1,3}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {3,0,2}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {3,0,2}, .tree_2 = { },
        .size = {1.0, 0.0, 1.0}, .normal = {0,+1,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 1, .face = +1,
        .l_degenerate = true
      },
      // test 48: same block, periodic boundaries
      {
        .rank = 3, .array = {3,5,1}, .periodic = {true,true,true},
        .levels_1 = 0, .array_1 = {2,4,0}, .tree_1 = { },
        .levels_2 = 0, .array_2 = {2,4,0}, .tree_2 = { },
        .size = {1.0, 1.0, 0.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
        .axis = 2, .face = -1,
        .l_degenerate = true
      },
      
    };

