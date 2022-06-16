// See LICENSE_CELLO file for license and copyright information

/// @file     test_Box.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-01-28
/// @brief    Test program for the Box class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Box");

  int n[] = {8,8,1};
  int g4[] = {4,4,0};
  int g0[] = {0,0,0};
  int g1[] = {1,1,0};
  // linear / enzo
  // refine parent-child (0,0)

  {
    int n3[3] = {2,2,2};
    int g3[3] = {1,1,1};
    Box adapt(2,n3,g3);
    BoxType box_type = BoxType_receive;
    int level = 0;
    int if3[3] = {0,-1,0};
    int ic3[3] = {0,0,0};
    adapt.set_block(box_type,level,if3,ic3);
    adapt.compute_region();
    adapt.print("adapt");
    int im3[3],ip3[3];
    bool lpad=false;
    adapt.get_start_stop(im3,ip3,BlockType::none,BlockType::receive,lpad);
    CkPrintf ("im3 %d %d %d ip3 %d %d %d\n",
              im3[0],im3[1],im3[2],
              ip3[0],ip3[1],ip3[2]);
    exit_();
  }
  { 
    Box * box000_l0 = new Box (2,n,g4);
    Box * box000_l1 = new Box (2,n,g4);
    Box * box000_l4 = new Box (2,n,g4);
    Box * box000_e0 = new Box (2,n,g4);
    Box * box000_e1 = new Box (2,n,g4);
    Box * box000_e4 = new Box (2,n,g4);

    unit_assert (box000_l0 != NULL);
    unit_assert (box000_l1 != NULL);
    unit_assert (box000_l4 != NULL);
    unit_assert (box000_e0 != NULL);
    unit_assert (box000_e1 != NULL);
    unit_assert (box000_e4 != NULL);

    int f3_000[3] = {0,0,0};
    int g3_000[3] = {0,0,0};

    box000_l0->set_block(BoxType_receive,1,f3_000,g3_000);
    box000_l1->set_block(BoxType_receive,1,f3_000,g3_000);
    box000_l4->set_block(BoxType_receive,1,f3_000,g3_000);
    box000_e0->set_block(BoxType_receive,1,f3_000,g3_000);
    box000_e1->set_block(BoxType_receive,1,f3_000,g3_000);
    box000_e4->set_block(BoxType_receive,1,f3_000,g3_000);

    box000_l0->set_recv_ghosts(g0);
    box000_l1->set_recv_ghosts(g1);
    box000_l4->set_recv_ghosts(g4);
    box000_e0->set_recv_ghosts(g0);
    box000_e1->set_recv_ghosts(g1);
    box000_e4->set_recv_ghosts(g4);

    // box000_l0->set_send_ghosts(g0);
    // box000_l1->set_send_ghosts(g1);
    // box000_l4->set_send_ghosts(g4);
    // box000_e0->set_send_ghosts(g0);
    // box000_e1->set_send_ghosts(g1);
    // box000_e4->set_send_ghosts(g4);

    box000_e0->set_padding(1);
    box000_e1->set_padding(1);
    box000_e4->set_padding(1);

    box000_l0->compute_block_start(BoxType_receive);
    box000_l1->compute_block_start(BoxType_receive);
    box000_l4->compute_block_start(BoxType_receive);
    box000_e0->compute_block_start(BoxType_receive);
    box000_e1->compute_block_start(BoxType_receive);
    box000_e4->compute_block_start(BoxType_receive);

    box000_l0->compute_region();
    box000_l1->compute_region();
    box000_l4->compute_region();
    box000_e0->compute_region();
    box000_e1->compute_region();
    box000_e4->compute_region();

    int if3_l0[3],nf3_l0[3];
    int ic3_l0[3],nc3_l0[3];
    int if3_l1[3],nf3_l1[3];
    int ic3_l1[3],nc3_l1[3];
    int if3_l4[3],nf3_l4[3];
    int ic3_l4[3],nc3_l4[3];
    int if3_e0[3],nf3_e0[3];
    int ic3_e0[3],nc3_e0[3];
    int if3_e1[3],nf3_e1[3];
    int ic3_e1[3],nc3_e1[3];
    int if3_e4[3],nf3_e4[3];
    int ic3_e4[3],nc3_e4[3];

    // box000_l0->print("box000_l0");
    // box000_l1->print("box000_l1");
    // box000_l4->print("box000_l4");
    // box000_e0->print("box000_e0");
    // box000_e1->print("box000_e1");
    // box000_e4->print("box000_e4");

    bool lpad;
    box000_l0->get_start_size
      (if3_l0,nf3_l0,BlockType::none,BlockType::receive,lpad=false);
    box000_l0->get_start_size
      (ic3_l0,nc3_l0,BlockType::none,   BlockType::send,lpad=true);
    box000_l1->get_start_size
      (if3_l1,nf3_l1,BlockType::none,BlockType::receive,lpad=false);
    box000_l1->get_start_size
      (ic3_l1,nc3_l1,BlockType::none,   BlockType::send,lpad=true);
    box000_l4->get_start_size
      (if3_l4,nf3_l4,BlockType::none,BlockType::receive,lpad=false);
    box000_l4->get_start_size
      (ic3_l4,nc3_l4,BlockType::none,   BlockType::send,lpad=true);
    box000_e0->get_start_size
      (if3_e0,nf3_e0,BlockType::none,BlockType::receive,lpad=false);
    box000_e0->get_start_size
      (ic3_e0,nc3_e0,BlockType::none,   BlockType::send,lpad=true);
    box000_e1->get_start_size
      (if3_e1,nf3_e1,BlockType::none,BlockType::receive,lpad=false);
    box000_e1->get_start_size
      (ic3_e1,nc3_e1,BlockType::none,   BlockType::send,lpad=true);
    box000_e4->get_start_size
      (if3_e4,nf3_e4,BlockType::none,BlockType::receive,lpad=false);
    box000_e4->get_start_size
      (ic3_e4,nc3_e4,BlockType::none,   BlockType::send,lpad=true);

    CkPrintf ("Linear      ghost = 0  face %d %d %d:",0,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_l0[0],ic3_l0[1], nc3_l0[0],nc3_l0[1],
              if3_l0[0],if3_l0[1], nf3_l0[0],nf3_l0[1]);
    CkPrintf ("Linear      ghost = 1  face %d %d %d:",0,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_l1[0],ic3_l1[1], nc3_l1[0],nc3_l1[1],
              if3_l1[0],if3_l1[1], nf3_l1[0],nf3_l1[1]);
    CkPrintf ("Linear      ghost = 4  face %d %d %d:",0,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_l4[0],ic3_l4[1], nc3_l4[0],nc3_l4[1],
              if3_l4[0],if3_l4[1], nf3_l4[0],nf3_l4[1]);
    CkPrintf ("EnzoProlong ghost = 0  face %d %d %d:",0,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_e0[0],ic3_e0[1], nc3_e0[0],nc3_e0[1],
              if3_e0[0],if3_e0[1], nf3_e0[0],nf3_e0[1]);
    CkPrintf ("EnzoProlong ghost = 1  face %d %d %d:",0,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_e1[0],ic3_e1[1], nc3_e1[0],nc3_e1[1],
              if3_e1[0],if3_e1[1], nf3_e1[0],nf3_e1[1]);
    CkPrintf ("EnzoProlong ghost = 4  face %d %d %d:",0,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_e4[0],ic3_e4[1], nc3_e4[0],nc3_e4[1],
              if3_e4[0],if3_e4[1], nf3_e4[0],nf3_e4[1]);
            
    //--------------------------------------------------
    // Linear      ghost = 0: ic  4  4 nc  8  8 if  4  4 nf  8  8
    // Linear      ghost = 1: ic  4  4 nc 10 10 if  4  4 nf  8  8
    // Linear      ghost = 4: ic  4  4 nc 12 12 if  4  4 nf 12 12
    // EnzoProlong ghost = 0: ic  4  4 nc 10 10 if  4  4 nf  8  8
    // EnzoProlong ghost = 1: ic  4  4 nc 12 12 if  4  4 nf  8  8
    // EnzoProlong ghost = 4: ic  4  4 nc 14 14 if  2  2 nf 14 14
    //--------------------------------------------------
                                                  
    delete box000_l0;
    delete box000_l1;
    delete box000_l4;
    delete box000_e0;
    delete box000_e1;
    delete box000_e4;

  }
  {

    Box * box100_l0 = new Box (2,n,g4);
    Box * box100_l1 = new Box (2,n,g4);
    Box * box100_l4 = new Box (2,n,g4);
    Box * box100_e0 = new Box (2,n,g4);
    Box * box100_e1 = new Box (2,n,g4);
    Box * box100_e4 = new Box (2,n,g4);

    unit_assert (box100_l0 != NULL);
    unit_assert (box100_l1 != NULL);
    unit_assert (box100_l4 != NULL);
    unit_assert (box100_e0 != NULL);
    unit_assert (box100_e1 != NULL);
    unit_assert (box100_e4 != NULL);

    int f3_100[3] = {1,0,0};
    int g3_100[3] = {0,0,0};
    
    box100_l0->set_block(BoxType_receive,1,f3_100,g3_100);
    box100_l1->set_block(BoxType_receive,1,f3_100,g3_100);
    box100_l4->set_block(BoxType_receive,1,f3_100,g3_100);
    box100_e0->set_block(BoxType_receive,1,f3_100,g3_100);
    box100_e1->set_block(BoxType_receive,1,f3_100,g3_100);
    box100_e4->set_block(BoxType_receive,1,f3_100,g3_100);

    box100_l0->set_recv_ghosts(g0);
    box100_l1->set_recv_ghosts(g1);
    box100_l4->set_recv_ghosts(g4);
    box100_e0->set_recv_ghosts(g0);
    box100_e1->set_recv_ghosts(g1);
    box100_e4->set_recv_ghosts(g4);

    // box100_l0->set_send_ghosts(g0);
    // box100_l1->set_send_ghosts(g1);
    // box100_l4->set_send_ghosts(g4);
    // box100_e0->set_send_ghosts(g0);
    // box100_e1->set_send_ghosts(g1);
    // box100_e4->set_send_ghosts(g4);

    box100_e0->set_padding(1);
    box100_e1->set_padding(1);
    box100_e4->set_padding(1);

    box100_l0->compute_block_start(BoxType_receive);
    box100_l1->compute_block_start(BoxType_receive);
    box100_l4->compute_block_start(BoxType_receive);
    box100_e0->compute_block_start(BoxType_receive);
    box100_e1->compute_block_start(BoxType_receive);
    box100_e4->compute_block_start(BoxType_receive);

    box100_l0->compute_region();
    box100_l1->compute_region();
    box100_l4->compute_region();
    box100_e0->compute_region();
    box100_e1->compute_region();
    box100_e4->compute_region();

    int if3_l0[3],nf3_l0[3];
    int ic3_l0[3],nc3_l0[3];
    int if3_l1[3],nf3_l1[3];
    int ic3_l1[3],nc3_l1[3];
    int if3_l4[3],nf3_l4[3];
    int ic3_l4[3],nc3_l4[3];
    int if3_e0[3],nf3_e0[3];
    int ic3_e0[3],nc3_e0[3];
    int if3_e1[3],nf3_e1[3];
    int ic3_e1[3],nc3_e1[3];
    int if3_e4[3],nf3_e4[3];
    int ic3_e4[3],nc3_e4[3];

    // box100_l0->print("box100_l0");
    // box100_l1->print("box100_l1");
    // box100_l4->print("box100_l4");
    // box100_e0->print("box100_e0");
    // box100_e1->print("box100_e1");
    // box100_e4->print("box100_e4");

    bool lpad;

    box100_l0->get_start_size
      (if3_l0,nf3_l0,BlockType::receive,BlockType::receive,lpad=false);
    box100_l0->get_start_size
      (ic3_l0,nc3_l0,BlockType::send,   BlockType::send,lpad=true);

    box100_l1->get_start_size
      (if3_l1,nf3_l1,BlockType::receive,BlockType::receive,lpad=false);
    box100_l1->get_start_size
      (ic3_l1,nc3_l1,BlockType::send,   BlockType::send,lpad=true);

    box100_l4->get_start_size
      (if3_l4,nf3_l4,BlockType::receive,BlockType::receive,lpad=false);
    box100_l4->get_start_size
      (ic3_l4,nc3_l4,BlockType::send,   BlockType::send,lpad=true);

    box100_e0->get_start_size
      (if3_e0,nf3_e0,BlockType::receive,BlockType::receive,lpad=false);
    box100_e0->get_start_size
      (ic3_e0,nc3_e0,BlockType::send,   BlockType::send,lpad=true);

    box100_e1->get_start_size
      (if3_e1,nf3_e1,BlockType::receive,BlockType::receive,lpad=false);
    box100_e1->get_start_size
      (ic3_e1,nc3_e1,BlockType::send,   BlockType::send,lpad=true);

    box100_e4->get_start_size
      (if3_e4,nf3_e4,BlockType::receive,BlockType::receive,lpad=false);
    box100_e4->get_start_size
      (ic3_e4,nc3_e4,BlockType::send,   BlockType::send,lpad=true);
    box100_e4->print("box100_e4");

    CkPrintf ("Linear      ghost = 0  face %d %d %d:",1,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_l0[0],ic3_l0[1], nc3_l0[0],nc3_l0[1],
              if3_l0[0],if3_l0[1], nf3_l0[0],nf3_l0[1]);
    CkPrintf ("Linear      ghost = 1  face %d %d %d:",1,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_l1[0],ic3_l1[1], nc3_l1[0],nc3_l1[1],
              if3_l1[0],if3_l1[1], nf3_l1[0],nf3_l1[1]);
    CkPrintf ("Linear      ghost = 4  face %d %d %d:",1,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_l4[0],ic3_l4[1], nc3_l4[0],nc3_l4[1],
              if3_l4[0],if3_l4[1], nf3_l4[0],nf3_l4[1]);
    CkPrintf ("EnzoProlong ghost = 0  face %d %d %d:",1,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_e0[0],ic3_e0[1], nc3_e0[0],nc3_e0[1],
              if3_e0[0],if3_e0[1], nf3_e0[0],nf3_e0[1]);
    CkPrintf ("EnzoProlong ghost = 1  face %d %d %d:",1,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_e1[0],ic3_e1[1], nc3_e1[0],nc3_e1[1],
              if3_e1[0],if3_e1[1], nf3_e1[0],nf3_e1[1]);
    CkPrintf ("EnzoProlong ghost = 4  face %d %d %d:",1,0,0);
    CkPrintf ("   ic %2d %2d nc %2d %2d if %2d %2d nf %2d %2d\n",
              ic3_e4[0],ic3_e4[1], nc3_e4[0],nc3_e4[1],
              if3_e4[0],if3_e4[1], nf3_e4[0],nf3_e4[1]);
            
    //--------------------------------------------------
    // Linear      ghost = 0: ic  4  4 nc  8  8 if  4  4 nf  8  8
    // Linear      ghost = 1: ic  4  4 nc 10 10 if  4  4 nf  8  8
    // Linear      ghost = 4: ic  4  4 nc 12 12 if  4  4 nf 12 12
    // EnzoProlong ghost = 0: ic  4  4 nc 10 10 if  4  4 nf  8  8
    // EnzoProlong ghost = 1: ic  4  4 nc 12 12 if  4  4 nf  8  8
    // EnzoProlong ghost = 4: ic  4  4 nc 14 14 if  2  2 nf 14 14
    //--------------------------------------------------
                                                  
    delete box100_l0;
    delete box100_l1;
    delete box100_l4;
    delete box100_e0;
    delete box100_e1;
    delete box100_e4;

  }

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

