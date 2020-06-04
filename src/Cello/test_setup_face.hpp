namespace test {
  Index create_index (const int levels, int array[3], int tree[][3])
{
  Index index;
  index.clear();
  index.set_array (array[0],array[1],array[2]);
  for (int i=0; i<levels; i++) {
    index.set_child(i+1,tree[i][0],tree[i][1],tree[i][2]);
  }
  index.set_level(levels);
  return index;
};

const int max_level = 3;
struct face_test_type {
  int test;
  int rank;  int array[3];   bool periodic[3];
  int levels_1;  int array_1[3];  int tree_1[max_level][3];
  int levels_2;  int array_2[3];  int tree_2[max_level][3];
  float size[3];  int normal[3];  int NX,NY,NZ;
  int ghost[3];
  int offset[3];
  int cx,cy,cz;
  double hx,hy,hz,dt;
  int axis; int face;
  bool l_degenerate;
};

  
  const int num_face_tests = 49;
  face_test_type face_test[] =
  {
    // Non-boundary
    // level 0 x-axis
    {
      .test = 0,
      .rank = 3, .array = {4,5,1}, .periodic = {false,false,false},
      .levels_1 = 0, .array_1 = {3,4,0}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {2,4,0}, .tree_2 = { },
      .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = false
    },
    // level 0 y-axis
    {
      .test = 1,
      .rank = 3, .array = {3,5,4}, .periodic = {false,false,false},
      .levels_1 = 0, .array_1 = {0,0,1}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {0,1,1}, .tree_2 = { },
      .size = {1.0, 0.0, 1.0}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = +1,
      .l_degenerate = false
    },
    // level 0 z-axis
    {
      .test = 2,
      .rank = 3, .array = {4,5,7}, .periodic = {false,false,false},
      .levels_1 = 0, .array_1 = {2,1,5}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {2,1,6}, .tree_2 = { },
      .size = {1.0, 1.0, 0.0}, .normal = {0,0,1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = +1,
      .l_degenerate = false
    },
    // same root same parent same level x-axis
    {
      .test = 3,
      .rank = 3, .array = {4,5,3}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {3,1,2}, .tree_1 = {{ 0, 1, 0}, {1, 0, 1}},
      .levels_2 = 2, .array_2 = {3,1,2}, .tree_2 = {{ 0, 1, 0}, {0, 0, 1}},
      .size = {0.0, 1.0, 1.0}, .normal = {1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = false
    },
    // same root same parent same level y-axis
    {
      .test = 4,
      .rank = 3, .array = {4,5,4}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {3,1,2}, .tree_1 = {{ 1, 1, 0}, {0, 1, 1}},
      .levels_2 = 2, .array_2 = {3,1,2}, .tree_2 = {{ 1, 1, 0}, {0, 0, 1}},
      .size = {1.0, 0.0, 1.0}, .normal = {0,-1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // same root same parent same level z-axis
    {
      .test = 5,
      .rank = 3, .array = {4,5,5}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {3,1,2}, .tree_1 = {{ 1, 0, 0}, {1, 1, 0}},
      .levels_2 = 2, .array_2 = {3,1,2}, .tree_2 = {{ 1, 0, 0}, {1, 1, 1}},
      .size = {1.0, 1.0, 0.0}, .normal = {0,0,1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = +1,
      .l_degenerate = false
    },
    // same root different parent same level x-axis
    {
      .test = 6,
      .rank = 3, .array = {3,1,3}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 1}, {1, 1, 0}},
      .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {0, 1, 0}},
      .size = {0.0, 1.0, 1.0}, .normal = {1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // same root different parent same level y-axis
    {
      .test = 7,
      .rank = 3, .array = {3,1,3}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 1, 1, 1}, {1, 0, 0}},
      .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 0, 1}, {1, 1, 0}},
      .size = {1.0, 0.0, 1.0}, .normal = {0,-1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // same root different parent same level z-axis
    {
      .test = 8,
      .rank = 3, .array = {2,2,4}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 1, 1, 0}, {1, 0, 1}},
      .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 0}},
      .size = {1.0, 1.0, 0.0}, .normal = {0,0,1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = +1,
      .l_degenerate = false
    },
    // different root same level x-axis
    {
      .test = 9,
      .rank = 3, .array = {2,5,3}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 1}, {0, 0, 1}},
      .levels_2 = 2, .array_2 = {0,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
      .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // different root same level y-axis
    {
      .test = 10,
      .rank = 3, .array = {3,5,7}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 0}, {1, 1, 0}},
      .levels_2 = 2, .array_2 = {1,1,2}, .tree_2 = {{ 0, 0, 0}, {1, 0, 0}},
      .size = {1.0, 0.0, 1.0}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // different root same level z-axis
    {
      .test = 11,
      .rank = 3, .array = {3,3,5}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {1,0,3}, .tree_1 = {{ 1, 1, 0}, {1, 0, 0}},
      .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
      .size = {1.0, 1.0, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = +1,
      .l_degenerate = false
    },
    // coarser level x-axis
    {
      .test = 12,
      .rank = 3, .array = {2,5,3}, .periodic = {false,false,false},
      .levels_1 = 1, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 1}},
      .levels_2 = 2, .array_2 = {0,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
      .size = {0.0, 0.5, 0.5}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,1},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // finer level x-axis
    {
      .test = 13,
      .rank = 3, .array = {2,5,3}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {0,2,2}, .tree_1 = {{ 1, 1, 1}, {1, 0, 1}},
      .levels_2 = 1, .array_2 = {1,2,2}, .tree_2 = {{ 0, 1, 1}},
      .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = false
    },
    // coarser level y-axis
    {
      .test = 14,
      .rank = 3, .array = {3,5,7}, .periodic = {false,false,false},
      .levels_1 = 1, .array_1 = {1,0,2}, .tree_1 = {{ 0, 1, 0}},
      .levels_2 = 2, .array_2 = {1,1,2}, .tree_2 = {{ 0, 0, 0}, {1, 0, 0}},
      .size = {0.5, 0.0, 0.5}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {1,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // finer level y-axis
    {
      .test = 15,
      .rank = 3, .array = {3,5,7}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {1,1,3}, .tree_1 = {{ 0, 0, 0}, {1, 0, 0}},
      .levels_2 = 1, .array_2 = {1,0,3}, .tree_2 = {{ 0, 1, 0}},
      .size = {1.0, 0.0, 1.0}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = +1,
      .l_degenerate = false
    },
    // coarser level z-axis
    {
      .test = 16,
      .rank = 3, .array = {3,3,5}, .periodic = {false,false,false},
      .levels_1 = 1, .array_1 = {1,0,3}, .tree_1 = {{ 1, 1, 0}},
      .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
      .size = {0.5, 0.5, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {1,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = +1,
      .l_degenerate = false
    },
    // finer level z-axis
    {
      .test = 17,
      .rank = 3, .array = {3,3,5}, .periodic = {false,false,false},
      .levels_1 = 2, .array_1 = {1,0,3}, .tree_1 = {{ 1, 1, 1}, {1, 0, 1}},
      .levels_2 = 1, .array_2 = {1,0,4}, .tree_2 = {{ 1, 1, 0}},
      .size = {1.0, 1.0, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = -1,
      .l_degenerate = false
    },
    // PERIODIC B.C.'S
    // periodic level 0 x-axis periodic boundary
    {
      .test = 18,
      .rank = 3, .array = {4,2,2}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {3,1,0}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
      .levels_2 = 0, .array_2 = {0,1,0}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
      .size = {0.0, 1.0, 1.0}, .normal = {1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // periodic level 0 y-axis
    {
      .test = 19,
      .rank = 3, .array = {3,5,4}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {0,0,1}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
      .levels_2 = 0, .array_2 = {0,4,1}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
      .size = {1.0, 0.0, 1.0}, .normal = {0,-1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // periodic level 0 z-axis
    {
      .test = 20,
      .rank = 3, .array = {4,5,7}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {2,1,0}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
      .levels_2 = 0, .array_2 = {2,1,6}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
      .size = {1.0, 1.0, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = -1,
      .l_degenerate = false
    },
    // periodic same level x-axis
    {
      .test = 21,
      .rank = 3, .array = {2,5,3}, .periodic = {true,true,true},
      .levels_1 = 2, .array_1 = {0,0,2}, .tree_1 = {{ 0, 1, 1}, {0, 0, 1}},
      .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
      .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = false
    },
    // periodic same level y-axis
    {
      .test = 22,
      .rank = 3, .array = {3,5,7}, .periodic = {true,true,true},
      .levels_1 = 2, .array_1 = {1,4,2}, .tree_1 = {{ 0, 1, 0}, {1, 1, 0}},
      .levels_2 = 2, .array_2 = {1,0,2}, .tree_2 = {{ 0, 0, 0}, {1, 0, 0}},
      .size = {1.0, 0.0, 1.0}, .normal = {0,1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = +1,
      .l_degenerate = false
    },
    // periodic same level z-axis
    {
      .test = 23,
      .rank = 3, .array = {3,3,5}, .periodic = {true,true,true},
      .levels_1 = 2, .array_1 = {1,0,0}, .tree_1 = {{ 1, 1, 0}, {1, 0, 0}},
      .levels_2 = 2, .array_2 = {1,0,4}, .tree_2 = {{ 1, 1, 1}, {1, 0, 1}},
      .size = {1.0, 1.0, 0.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = -1,
      .l_degenerate = false
    },
    // non-periodic same level x-axis
    {
      .test = 24,
      .rank = 3, .array = {4,2,2}, .periodic = {false,false,false},
      .levels_1 = 0, .array_1 = {3,1,0}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
      .levels_2 = 0, .array_2 = {0,1,0}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
      .size = {-1.0, 1.0, 1.0}, .normal = {1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // non-periodic level 0 y-axis
    {
      .test = 25,
      .rank = 3, .array = {3,5,4}, .periodic = {false,false,false},
      .levels_1 = 0, .array_1 = {0,0,1}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
      .levels_2 = 0, .array_2 = {0,4,1}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
      .size = {1.0, -1.0, 1.0}, .normal = {0,-1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // periodic level 0 z-axis
    {
      .test = 26,
      .rank = 3, .array = {4,5,7}, .periodic = {false,false,false},
      .levels_1 = 0, .array_1 = {2,1,0}, .tree_1 = {{-1,-1,-1}, {-1,-1,-1}},
      .levels_2 = 0, .array_2 = {2,1,6}, .tree_2 = {{-1,-1,-1}, {-1,-1,-1}},
      .size = {1.0, 1.0, -1.0}, .normal = {0,0,-1}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = -1,
      .l_degenerate = false
    },
    //======================================================================
    // 2D tests
    //======================================================================
    // level 0 x-axis
    {
      .test = 27,
      .rank = 2, .array = {4,5}, .periodic = {false,false},
      .levels_1 = 0, .array_1 = {3,4}, .tree_1 = { {-1,-1}, {-1,-1} },
      .levels_2 = 0, .array_2 = {2,4}, .tree_2 = { {-1,-1}, {-1,-1} },
      .size = {0.0, 1.0}, .normal = {-1,0}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = false
    },
    // level 0 y-axis
    {
      .test = 28,
      .rank = 2, .array = {3,5}, .periodic = {false,false},
      .levels_1 = 0, .array_1 = {0,0}, .tree_1 = {{-1,-1}, {-1,-1}},
      .levels_2 = 0, .array_2 = {0,1}, .tree_2 = {{-1,-1}, {-1,-1}},
      .size = {1.0, 0.0}, .normal = {0,1}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = +1,
      .l_degenerate = false
    },
    // same root same parent same level x-axis
    {
      .test = 29,
      .rank = 2, .array = {4,5}, .periodic = {false,false},
      .levels_1 = 2, .array_1 = {3,1}, .tree_1 = {{ 0, 1}, {1, 0}},
      .levels_2 = 2, .array_2 = {3,1}, .tree_2 = {{ 0, 1}, {0, 0}},
      .size = {0.0, 1.0}, .normal = {1,0}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = false
    },
    // same root same parent same level y-axis
    {
      .test = 30,
      .rank = 2, .array = {4,5}, .periodic = {false,false},
      .levels_1 = 2, .array_1 = {3,1}, .tree_1 = {{ 1, 1}, {0, 1}},
      .levels_2 = 2, .array_2 = {3,1}, .tree_2 = {{ 1, 1}, {0, 0}},
      .size = {1.0, 0.0}, .normal = {0,-1}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // same root different parent same level x-axis
    {
      .test = 31,
      .rank = 2, .array = {3,1}, .periodic = {false,false},
      .levels_1 = 2, .array_1 = {1,0}, .tree_1 = {{ 0, 1}, {1, 1}},
      .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 1, 1}, {0, 1}},
      .size = {0.0, 1.0}, .normal = {1,0}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // same root different parent same level y-axis
    {
      .test = 32,
      .rank = 2, .array = {3,1}, .periodic = {false,false},
      .levels_1 = 2, .array_1 = {1,0}, .tree_1 = {{ 1, 1}, {1, 0}},
      .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 1, 0}, {1, 1}},
      .size = {1.0, 0.0}, .normal = {0,-1}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // different root same level x-axis
    {
      .test = 33,
      .rank = 2, .array = {2,5}, .periodic = {false,false},
      .levels_1 = 2, .array_1 = {1,0}, .tree_1 = {{ 0, 1}, {0, 0}},
      .levels_2 = 2, .array_2 = {0,0}, .tree_2 = {{ 1, 1}, {1, 0}},
      .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // different root same level y-axis
    {
      .test = 34,
      .rank = 2, .array = {3,5}, .periodic = {false,false},
      .levels_1 = 2, .array_1 = {1,0}, .tree_1 = {{ 0, 1}, {1, 1}},
      .levels_2 = 2, .array_2 = {1,1}, .tree_2 = {{ 0, 0}, {1, 0}},
      .size = {1.0, 0.0}, .normal = {0,1}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // Periodic b.c.'s
    // periodic level 0 x-axis periodic boundary
    {
      .test = 35,
      .rank = 2, .array = {4,2}, .periodic = {true,true},
      .levels_1 = 0, .array_1 = {3,1}, .tree_1 = {{-1,-1}, {-1,-1}},
      .levels_2 = 0, .array_2 = {0,1}, .tree_2 = {{-1,-1}, {-1,-1}},
      .size = {0.0, 1.0}, .normal = {1,0}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // periodic level 0 y-axis
    {
      .test = 36,
      .rank = 2, .array = {3,5}, .periodic = {true,true},
      .levels_1 = 0, .array_1 = {0,0}, .tree_1 = {{-1,-1}, {-1,-1}},
      .levels_2 = 0, .array_2 = {0,4}, .tree_2 = {{-1,-1}, {-1,-1}},
      .size = {1.0, 0.0}, .normal = {0,-1}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // periodic same level x-axis
    {
      .test = 37,
      .rank = 2, .array = {2,5}, .periodic = {true,true},
      .levels_1 = 2, .array_1 = {0,0}, .tree_1 = {{ 0, 1}, {0, 0}},
      .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 1, 1}, {1, 0}},
      .size = {0.0, 1.0}, .normal = {-1,0}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = false
    },
    // periodic same level y-axis
    {
      .test = 38,
      .rank = 2, .array = {3,5}, .periodic = {true,true},
      .levels_1 = 2, .array_1 = {1,4}, .tree_1 = {{ 0, 1}, {1, 1}},
      .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 0, 0}, {1, 0}},
      .size = {1.0, 0.0}, .normal = {0,1}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = +1,
      .l_degenerate = false
    },
    // Non-periodic b.c.'s
    // periodic level 0 x-axis non-periodic boundary
    {
      .test = 39,
      .rank = 2, .array = {4,2}, .periodic = {false,false},
      .levels_1 = 0, .array_1 = {3,1}, .tree_1 = {{-1,-1}, {-1,-1}},
      .levels_2 = 0, .array_2 = {0,1}, .tree_2 = {{-1,-1}, {-1,-1}},
      .size = {-1.0, 1.0}, .normal = {1,0}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = +1,
      .l_degenerate = false
    },
    // periodic level 0 y-axis non-periodic boundary
    {
      .test = 40,
      .rank = 2, .array = {3,5}, .periodic = {false,false},
      .levels_1 = 0, .array_1 = {0,0}, .tree_1 = {{-1,-1}, {-1,-1}},
      .levels_2 = 0, .array_2 = {0,4}, .tree_2 = {{-1,-1}, {-1,-1}},
      .size = {1.0,-1.0}, .normal = {0,-1}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = -1,
      .l_degenerate = false
    },
    // periodic same level x-axis non-periodic boundary
    {
      .test = 41,
      .rank = 2, .array = {2,5}, .periodic = {false,false},
      .levels_1 = 2, .array_1 = {0,0}, .tree_1 = {{ 0, 1}, {0, 0}},
      .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 1, 1}, {1, 0}},
      .size = {-1.0, 1.0}, .normal = {-1,0}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = false
    },
    // periodic same level y-axis non-periodic boundary
    {
      .test = 42,
      .rank = 2, .array = {3,5}, .periodic = {false,false},
      .levels_1 = 2, .array_1 = {1,4}, .tree_1 = {{ 0, 1}, {1, 1}},
      .levels_2 = 2, .array_2 = {1,0}, .tree_2 = {{ 0, 0}, {1, 0}},
      .size = {1.0, -1.0}, .normal = {0,1}, .NX = 10, .NY = 6, .NZ = 1,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = +1,
      .l_degenerate = false
    },
    // same block, periodic boundaries
    {
      .test = 43,
      .rank = 3, .array = {1,5,3}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {0,2,1}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {0,2,1}, .tree_2 = { },
      .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = true
    },
    // same block, periodic boundaries
    {
      .test = 44,
      .rank = 3, .array = {5,1,3}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {3,0,2}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {3,0,2}, .tree_2 = { },
      .size = {1.0, 0.0, 1.0}, .normal = {0,+1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = +1,
      .l_degenerate = true
    },
    // same block, periodic boundaries
    {
      .test = 45,
      .rank = 3, .array = {3,5,1}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {2,4,0}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {2,4,0}, .tree_2 = { },
      .size = {1.0, 1.0, 0.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = -1,
      .l_degenerate = true
    },
    // same block, periodic boundaries
    {
      .test = 46,
      .rank = 3, .array = {1,5,3}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {0,2,1}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {0,2,1}, .tree_2 = { },
      .size = {0.0, 1.0, 1.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 0, .face = -1,
      .l_degenerate = true
    },
    // same block, periodic boundaries
    {
      .test = 47,
      .rank = 3, .array = {5,1,3}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {3,0,2}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {3,0,2}, .tree_2 = { },
      .size = {1.0, 0.0, 1.0}, .normal = {0,+1,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 1, .face = +1,
      .l_degenerate = true
    },
    // same block, periodic boundaries
    {
      .test = 48,
      .rank = 3, .array = {3,5,1}, .periodic = {true,true,true},
      .levels_1 = 0, .array_1 = {2,4,0}, .tree_1 = { },
      .levels_2 = 0, .array_2 = {2,4,0}, .tree_2 = { },
      .size = {1.0, 1.0, 0.0}, .normal = {-1,0,0}, .NX = 10, .NY = 6, .NZ = 8,
      .ghost = {0,0,0},
      .offset = {0,0,0},
      .cx = 0, .cy = 0, .cz = 0,
      .hx = 0.25, .hy = 0.5, .hz = 1.0, .dt = 0.125,
      .axis = 2, .face = -1,
      .l_degenerate = true
    },
      
  };

};
