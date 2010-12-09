// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-19
/// @brief    Unit tests for the Layout class

#include <math.h>
#include "cello.hpp"

#include "error.hpp"
#include "test.hpp"
#include "parallel.hpp"

#define TOL 2e-16

#include "parallel.def"
#include PARALLEL_CHARM_INCLUDE(test_Layout.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  unit_init (parallel->rank(), parallel->size());

  unit_class("Layout");

  //----------------------------------------------------------------------
  // index conversions
  //----------------------------------------------------------------------

  {
    unit_func("index_*");

    int nx=5, ny=3, nz=4;
    bool passed = true;
    for (int iz=0; iz<nz && passed; iz++) {
      for (int iy=0; iy<ny && passed; iy++) {
	for (int ix=0; ix<nx && passed; ix++) {

	  int i,j,jx,jy,jz,kx,ky,kz;

	  // ix,iy,iz -> i
	  j = index_3_to_1(ix,iy,iz,nx,ny,nz);
	  index_1_to_3(j,jx,jy,jz,nx,ny,nz);
	  passed = passed && (ix==jx && iy==jy && iz==jz);

	  // i -> ix, i -> iy, i -> iz
	  kx = index_1_to_x(j,nx,ny,nz);
	  ky = index_1_to_y(j,nx,ny,nz);
	  kz = index_1_to_z(j,nx,ny,nz);
	  passed = passed && (ix==kx && iy==ky && iz==kz);

	  i = ix + nx*(iy + ny*iz);

	  // i -> ix,iy,iz
	  index_1_to_3(i,jx,jy,jz,nx,ny,nz);
	  j = index_3_to_1(jx,jy,jz,nx,ny,nz);
	  passed = passed && (i == j);
	
	}
      }
    }
    unit_assert (passed);
  }

  //----------------------------------------------------------------------
  // serial layout: (processes,threads,data blocks) = (1,1,1)
  //----------------------------------------------------------------------

  {
    
    unit_func("Layout");
    Layout layout;
    unit_assert (true);

    layout.set_periodic(axis_x,true);
    layout.set_periodic(axis_y,true);
    layout.set_periodic(axis_z,true);

    unit_func("processes");

    int p0,p1,p2;
    unit_assert (layout.processes(&p0,&p1,&p2) == 1);
    unit_assert (p0*p1*p2==1);

    unit_func("threads");

    int t0,t1,t2;
    unit_assert (layout.threads(&t0,&t1,&t2)  == 1);
    unit_assert (t0*t1*t2==1);

    unit_func("blocks_per_process");

    int bp0,bp1,bp2;

    unit_assert (layout.blocks_per_process(&bp0,&bp1,&bp2)  == 1);

    unit_func("blocks_per_thread");

    int bt0,bt1,bt2;
    unit_assert (layout.blocks_per_thread(&bt0,&bt1,&bt2)  == 1);

    unit_func("is_periodic");
    unit_assert (layout.is_periodic(axis_x) == true);
    unit_assert (layout.is_periodic(axis_y) == true);
    unit_assert (layout.is_periodic(axis_z) == true);

    unit_func("neighbor_is_internal");
    // periodic, so all neighbors should be internal
    unit_assert (layout.neighbor_is_internal(0,0,0,axis_x,+1));
    unit_assert (layout.neighbor_is_internal(0,0,0,axis_x,-1));
    unit_assert (layout.neighbor_is_internal(0,0,0,axis_y,+1));
    unit_assert (layout.neighbor_is_internal(0,0,0,axis_y,-1));
    unit_assert (layout.neighbor_is_internal(0,0,0,axis_z,+1));
    unit_assert (layout.neighbor_is_internal(0,0,0,axis_z,-1));

    unit_func("neighbor_process");
    unit_assert(layout.neighbor_process(0,0,0,axis_x,-1) == 0);
    unit_assert(layout.neighbor_process(0,0,0,axis_x,+1) == 0);
    unit_assert(layout.neighbor_process(0,0,0,axis_y,-1) == 0);
    unit_assert(layout.neighbor_process(0,0,0,axis_y,+1) == 0);
    unit_assert(layout.neighbor_process(0,0,0,axis_z,-1) == 0);
    unit_assert(layout.neighbor_process(0,0,0,axis_z,+1) == 0);

    unit_func("neighbor_thread");
    unit_assert(layout.neighbor_thread(0,0,0,axis_x,-1) == 0);
    unit_assert(layout.neighbor_thread(0,0,0,axis_x,+1) == 0);
    unit_assert(layout.neighbor_thread(0,0,0,axis_y,-1) == 0);
    unit_assert(layout.neighbor_thread(0,0,0,axis_y,+1) == 0);
    unit_assert(layout.neighbor_thread(0,0,0,axis_z,-1) == 0);
    unit_assert(layout.neighbor_thread(0,0,0,axis_z,+1) == 0);

    unit_func("extent");
    double lower_extent[3],upper_extent[3];
    layout.extent(0,0,0,lower_extent,upper_extent);
    unit_assert(lower_extent[axis_x] == 0.0);
    unit_assert(lower_extent[axis_y] == 0.0);
    unit_assert(lower_extent[axis_z] == 0.0);
    unit_assert(upper_extent[axis_x] == 1.0);
    unit_assert(upper_extent[axis_y] == 1.0);
    unit_assert(upper_extent[axis_z] == 1.0);

    unit_func("array_indices");
    int index_lower[3],index_upper[3];
    int nx = 13;
    int ny = 7;
    int nz = 15;
    layout.array_indices(0,0,0,nx,ny,nz,index_lower,index_upper);
    unit_assert(index_lower[axis_x] == 0);
    unit_assert(index_lower[axis_y] == 0);
    unit_assert(index_lower[axis_z] == 0);
    unit_assert(index_upper[axis_x] == nx);
    unit_assert(index_upper[axis_y] == ny);
    unit_assert(index_upper[axis_z] == nz);
  }  

  //----------------------------------------------------------------------
  // parallel layout: (processes,threads,data blocks) = (NP,1,1)
  //----------------------------------------------------------------------

  {

    int pb3[3] = {5,3,7};  // processor blocks
    int npb = pb3[0]*pb3[1]*pb3[2];

    unit_func("Layout");
    Layout layout;
    unit_assert (true);

    layout.set_periodic(axis_x,false);
    layout.set_periodic(axis_y,false);
    layout.set_periodic(axis_z,false);

    layout.set_processes(pb3[0],pb3[1],pb3[2]);

    unit_func("processes");

    int p0,p1,p2;
    unit_assert (layout.processes(&p0,&p1,&p2) == npb);
    unit_assert (p0*p1*p2==npb);

    unit_func("threads");

    int t0,t1,t2;
    unit_assert (layout.threads(&t0,&t1,&t2)  == 1);
    unit_assert (t0*t1*t2==1);

    unit_func("blocks_per_process");

    int bp0,bp1,bp2;
    unit_assert (layout.blocks_per_process(&bp0,&bp1,&bp2)  == 1);

    unit_func("blocks_per_thread");

    int bt0,bt1,bt2;
    unit_assert (layout.blocks_per_thread(&bt0,&bt1,&bt2)  == 1);

    unit_func("is_periodic");
    unit_assert (layout.is_periodic(axis_x) == false);
    unit_assert (layout.is_periodic(axis_y) == false);
    unit_assert (layout.is_periodic(axis_z) == false);

    //----------------------------------------------------------------------

    unit_func("neighbor_is_internal");

    // periodic, so all neighbors should be internal

    bool passed = true;

    for (int ipz=0; ipz<pb3[2]; ipz++) {
      for (int ipy=0; ipy<pb3[1]; ipy++) {
	for (int ipx=0; ipx<pb3[0]; ipx++) {
	  int ip = ipx + pb3[0]*(ipy + pb3[1]*ipz);
	  bool lxm = layout.neighbor_is_internal(ip,0,0,axis_x,-1);
	  bool lxp = layout.neighbor_is_internal(ip,0,0,axis_x,+1);
	  bool lym = layout.neighbor_is_internal(ip,0,0,axis_y,-1);
	  bool lyp = layout.neighbor_is_internal(ip,0,0,axis_y,+1);
	  bool lzm = layout.neighbor_is_internal(ip,0,0,axis_z,-1);
	  bool lzp = layout.neighbor_is_internal(ip,0,0,axis_z,+1);
	  passed = passed && (lxp == (ipx != pb3[0]-1));
 	  passed = passed && (lxm == (ipx != 0));
 	  passed = passed && (lyp == (ipy != pb3[1]-1));
 	  passed = passed && (lym == (ipy != 0));
 	  passed = passed && (lzp == (ipz != pb3[2]-1));
 	  passed = passed && (lzm == (ipz != 0));
	}
      }
    }

    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("neighbor_process");

    passed = true;

    for (int ipz=0; ipz<pb3[2]; ipz++) {
      for (int ipy=0; ipy<pb3[1]; ipy++) {
	for (int ipx=0; ipx<pb3[0]; ipx++) {
	  int ip = ipx + pb3[0]*(ipy + pb3[1]*ipz);
	  bool lxm = layout.neighbor_is_internal(ip,0,0,axis_x,-1);
	  bool lxp = layout.neighbor_is_internal(ip,0,0,axis_x,+1);
	  bool lym = layout.neighbor_is_internal(ip,0,0,axis_y,-1);
	  bool lyp = layout.neighbor_is_internal(ip,0,0,axis_y,+1);
	  bool lzm = layout.neighbor_is_internal(ip,0,0,axis_z,-1);
	  bool lzp = layout.neighbor_is_internal(ip,0,0,axis_z,+1);
	  int ipxm = layout.neighbor_process(ip,0,0,axis_x,-1);
	  int ipxp = layout.neighbor_process(ip,0,0,axis_x,+1);
	  int ipym = layout.neighbor_process(ip,0,0,axis_y,-1);
	  int ipyp = layout.neighbor_process(ip,0,0,axis_y,+1);
	  int ipzm = layout.neighbor_process(ip,0,0,axis_z,-1);
	  int ipzp = layout.neighbor_process(ip,0,0,axis_z,+1);
	  if (lxm) passed = passed && (ip+ipxm == (pb3[0]+ipx-1)%pb3[0]);
 	  if (lxp) passed = passed && (ip+ipxp == (pb3[0]+ipx+1)%pb3[0]);
 	  if (lym) passed = passed && (ip+ipym == (pb3[1]+ipy-1)%pb3[1]);
 	  if (lyp) passed = passed && (ip+ipyp == (pb3[1]+ipy+1)%pb3[1]);
 	  if (lzm) passed = passed && (ip+ipzm == (pb3[2]+ipz-1)%pb3[2]);
 	  if (lzp) passed = passed && (ip+ipzp == (pb3[2]+ipz+1)%pb3[2]);
	}
      }
    }
    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("neighbor_thread");

    passed = true;

    for (int ipz=0; ipz<pb3[2]; ipz++) {
      for (int ipy=0; ipy<pb3[1]; ipy++) {
	for (int ipx=0; ipx<pb3[0]; ipx++) {
	  int ip = ipx + pb3[0]*(ipy + pb3[1]*ipz);
	  int itxm = layout.neighbor_thread(ip,0,0,axis_x,-1);
	  int itxp = layout.neighbor_thread(ip,0,0,axis_x,+1);
	  int itym = layout.neighbor_thread(ip,0,0,axis_y,-1);
	  int ityp = layout.neighbor_thread(ip,0,0,axis_y,+1);
	  int itzm = layout.neighbor_thread(ip,0,0,axis_z,-1);
	  int itzp = layout.neighbor_thread(ip,0,0,axis_z,+1);
	  passed = passed && (itxm == 0);
	  passed = passed && (itxp == 0);
	  passed = passed && (itym == 0);
	  passed = passed && (ityp == 0);
	  passed = passed && (itzm == 0);
	  passed = passed && (itzp == 0);
	}
      }
    }

    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("extent");
    passed = true;

    for (int ipz=0; ipz<pb3[2]; ipz++) {
      for (int ipy=0; ipy<pb3[1]; ipy++) {
	for (int ipx=0; ipx<pb3[0]; ipx++) {
	  int ip = ipx + pb3[0]*(ipy + pb3[1]*ipz);
	  double lower_extent[3],upper_extent[3];
	  layout.extent(ip,0,0,lower_extent,upper_extent);
 	  passed = passed && (fabs(lower_extent[axis_x] - 1.0*ipx/pb3[0])<TOL);
 	  passed = passed && (fabs(lower_extent[axis_y] - 1.0*ipy/pb3[1])<TOL);
  	  passed = passed && (fabs(lower_extent[axis_z] - 1.0*ipz/pb3[2])<TOL);
	  passed = passed && (fabs(upper_extent[axis_x] - 1.0*(ipx+1)/pb3[0])<TOL);
	  passed = passed && (fabs(upper_extent[axis_y] - 1.0*(ipy+1)/pb3[1])<TOL);
	  passed = passed && (fabs(upper_extent[axis_z] - 1.0*(ipz+1)/pb3[2])<TOL);
	}
      }
    }
    // NOTE: results sensitive to roundoff
    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("array_indices");

    // array size

    int nx = 158;
    int ny = 720;
    int nz = 48;

    passed = true;

    for (int ipz=0; ipz<pb3[2]; ipz++) {
      for (int ipy=0; ipy<pb3[1]; ipy++) {
	for (int ipx=0; ipx<pb3[0]; ipx++) {
	  int ip = ipx + pb3[0]*(ipy + pb3[1]*ipz);
	  int index_lower[3],index_upper[3];
	  layout.array_indices(ip,0,0,nx,ny,nz,index_lower,index_upper);
	  passed = passed && (index_lower[axis_x] == nx * ipx/pb3[0]);
	  passed = passed && (index_lower[axis_y] == ny * ipy/pb3[1]);
	  passed = passed && (index_lower[axis_z] == nz * ipz/pb3[2]);
	  passed = passed && (index_upper[axis_x] == nx * (ipx+1)/pb3[0]);
	  passed = passed && (index_upper[axis_y] == ny * (ipy+1)/pb3[1]);
	  passed = passed && (index_upper[axis_z] == nz * (ipz+1)/pb3[2]);
	}
      }
    }

    unit_assert (passed);

  }  

  //----------------------------------------------------------------------
  // blocked layout: (processes,threads,data blocks) = (1,1,NB)
  //----------------------------------------------------------------------

  {

    int db3[3] = {4,1,9};  // data blocks
    int ndb = db3[0]*db3[1]*db3[2];

    unit_func("Layout");
    Layout layout;
    unit_assert (true);

    layout.set_periodic(axis_x,false);
    layout.set_periodic(axis_y,false);
    layout.set_periodic(axis_z,false);

    layout.set_blocks(db3[0],db3[1],db3[2]);

    unit_func("processes");
    int p0,p1,p2;
    unit_assert (layout.processes(&p0,&p1,&p2) == 1);
    unit_assert (p0*p1*p2==1);

    unit_func("threads");

    int t0,t1,t2;
    unit_assert (layout.threads(&t0,&t1,&t2)  == 1);
    unit_assert (t0*t1*t2==1);

    unit_func("blocks_per_process");

    int bp0,bp1,bp2;
    unit_assert (layout.blocks_per_process(&bp0,&bp1,&bp2)  == ndb);

    unit_func("blocks_per_thread");

    int bt0,bt1,bt2;
    unit_assert (layout.blocks_per_thread(&bt0,&bt1,&bt2)  == ndb);

    unit_func("is_periodic");
    unit_assert (layout.is_periodic(axis_x) == false);
    unit_assert (layout.is_periodic(axis_y) == false);
    unit_assert (layout.is_periodic(axis_z) == false);

    //----------------------------------------------------------------------

    unit_func("neighbor_is_internal");

    // periodic, so all neighbors should be internal

    bool passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
	  bool lxm = layout.neighbor_is_internal(0,0,idb,axis_x,-1);
	  bool lxp = layout.neighbor_is_internal(0,0,idb,axis_x,+1);
	  bool lym = layout.neighbor_is_internal(0,0,idb,axis_y,-1);
	  bool lyp = layout.neighbor_is_internal(0,0,idb,axis_y,+1);
	  bool lzm = layout.neighbor_is_internal(0,0,idb,axis_z,-1);
	  bool lzp = layout.neighbor_is_internal(0,0,idb,axis_z,+1);
	  passed = passed && (lxp == (idbx != db3[0]-1));
 	  passed = passed && (lxm == (idbx != 0));
 	  passed = passed && (lyp == (idby != db3[1]-1));
 	  passed = passed && (lym == (idby != 0));
 	  passed = passed && (lzp == (idbz != db3[2]-1));
 	  passed = passed && (lzm == (idbz != 0));
	}
      }
    }

    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("neighbor_process");

    passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
	  bool lxm = layout.neighbor_is_internal(0,0,idb,axis_x,-1);
	  bool lxp = layout.neighbor_is_internal(0,0,idb,axis_x,+1);
	  bool lym = layout.neighbor_is_internal(0,0,idb,axis_y,-1);
	  bool lyp = layout.neighbor_is_internal(0,0,idb,axis_y,+1);
	  bool lzm = layout.neighbor_is_internal(0,0,idb,axis_z,-1);
	  bool lzp = layout.neighbor_is_internal(0,0,idb,axis_z,+1);
	  int ipbxm = layout.neighbor_process(0,0,idb,axis_x,-1);
	  int ipbxp = layout.neighbor_process(0,0,idb,axis_x,+1);
	  int ipbym = layout.neighbor_process(0,0,idb,axis_y,-1);
	  int ipbyp = layout.neighbor_process(0,0,idb,axis_y,+1);
	  int ipbzm = layout.neighbor_process(0,0,idb,axis_z,-1);
	  int ipbzp = layout.neighbor_process(0,0,idb,axis_z,+1);
	  if (lxm) passed = passed && (ipbxm == 0);
 	  if (lxp) passed = passed && (ipbxp == 0);
 	  if (lym) passed = passed && (ipbym == 0);
 	  if (lyp) passed = passed && (ipbyp == 0);
 	  if (lzm) passed = passed && (ipbzm == 0);
 	  if (lzp) passed = passed && (ipbzp == 0);
	}
      }
    }
    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("neighbor_thread");

    passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
	  int itxm = layout.neighbor_thread(0,0,idb,axis_x,-1);
	  int itxp = layout.neighbor_thread(0,0,idb,axis_x,+1);
	  int itym = layout.neighbor_thread(0,0,idb,axis_y,-1);
	  int ityp = layout.neighbor_thread(0,0,idb,axis_y,+1);
	  int itzm = layout.neighbor_thread(0,0,idb,axis_z,-1);
	  int itzp = layout.neighbor_thread(0,0,idb,axis_z,+1);
	  passed = passed && (itxm == 0);
	  passed = passed && (itxp == 0);
	  passed = passed && (itym == 0);
	  passed = passed && (ityp == 0);
	  passed = passed && (itzm == 0);
	  passed = passed && (itzp == 0);
	}
      }
    }

    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("extent");
    passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
	  double lower_extent[3],upper_extent[3];
	  layout.extent(0,0,idb,lower_extent,upper_extent);
 	  passed = passed && (fabs(lower_extent[axis_x] - 1.0*idbx/db3[0])<TOL);
 	  passed = passed && (fabs(lower_extent[axis_y] - 1.0*idby/db3[1])<TOL);
  	  passed = passed && (fabs(lower_extent[axis_z] - 1.0*idbz/db3[2])<TOL);
	  passed = passed && (fabs(upper_extent[axis_x] - 1.0*(idbx+1)/db3[0])<TOL);
	  passed = passed && (fabs(upper_extent[axis_y] - 1.0*(idby+1)/db3[1])<TOL);
	  passed = passed && (fabs(upper_extent[axis_z] - 1.0*(idbz+1)/db3[2])<TOL);
	}
      }
    }
    // NOTE: results sensitive to roundoff
    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("array_indices");

    // array size

    int nx = 158;
    int ny = 720;
    int nz = 48;

    passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
	  int index_lower[3],index_upper[3];
	  layout.array_indices(0,0,idb,nx,ny,nz,index_lower,index_upper);
	  passed = passed && (index_lower[axis_x] == nx * idbx/db3[0]);
 	  passed = passed && (index_lower[axis_y] == ny * idby/db3[1]);
 	  passed = passed && (index_lower[axis_z] == nz * idbz/db3[2]);
 	  passed = passed && (index_upper[axis_x] == nx * (idbx+1)/db3[0]);
 	  passed = passed && (index_upper[axis_y] == ny * (idby+1)/db3[1]);
 	  passed = passed && (index_upper[axis_z] == nz * (idbz+1)/db3[2]);
	}
      }
    }

    unit_assert (passed);

  }  

  //----------------------------------------------------------------------
  // parallel blocked layout: (processes,threads,data blocks) = (NP,1,NB)
  //----------------------------------------------------------------------

  {

    int pb3[3] = {5,3,7};  // processor blocks
    int db3[3] = {4,1,9};  // data blocks

    int npb = pb3[0]*pb3[1]*pb3[2];
    int ndb = db3[0]*db3[1]*db3[2];

    unit_func("Layout");
    Layout layout;
    unit_assert (true);

    layout.set_periodic(axis_x,false);
    layout.set_periodic(axis_y,false);
    layout.set_periodic(axis_z,false);

    layout.set_processes(pb3[0],pb3[1],pb3[2]);
    layout.set_blocks   (db3[0],db3[1],db3[2]);

    unit_func("processes");

    int p0,p1,p2;
    unit_assert (layout.processes(&p0,&p1,&p2) == npb);
    unit_assert (p0*p1*p2==npb);

    unit_func("threads");

    int t0,t1,t2;
    unit_assert (layout.threads(&t0,&t1,&t2)  == 1);
    unit_assert (t0*t1*t2==1);

    unit_func("blocks_per_process");

    int bp0,bp1,bp2;
    unit_assert (layout.blocks_per_process(&bp0,&bp1,&bp2)  == ndb);

    unit_func("blocks_per_thread");

    int bt0,bt1,bt2;
    unit_assert (layout.blocks_per_thread(&bt0,&bt1,&bt2)  == ndb);

    unit_func("is_periodic");
    unit_assert (layout.is_periodic(axis_x) == false);
    unit_assert (layout.is_periodic(axis_y) == false);
    unit_assert (layout.is_periodic(axis_z) == false);

    //----------------------------------------------------------------------

    unit_func("neighbor_is_internal");

    // periodic, so all neighbors should be internal

    bool passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
    for (int ipbz=0; ipbz<pb3[2]; ipbz++) {
      for (int ipby=0; ipby<pb3[1]; ipby++) {
	for (int ipbx=0; ipbx<pb3[0]; ipbx++) {
	  int ipb = ipbx + pb3[0]*(ipby + pb3[1]*ipbz);
	  bool lxm = layout.neighbor_is_internal(ipb,0,idb,axis_x,-1);
	  bool lxp = layout.neighbor_is_internal(ipb,0,idb,axis_x,+1);
	  bool lym = layout.neighbor_is_internal(ipb,0,idb,axis_y,-1);
	  bool lyp = layout.neighbor_is_internal(ipb,0,idb,axis_y,+1);
	  bool lzm = layout.neighbor_is_internal(ipb,0,idb,axis_z,-1);
	  bool lzp = layout.neighbor_is_internal(ipb,0,idb,axis_z,+1);
	  passed = passed && (lxp == ! (idbx == db3[0]-1 && ipbx == pb3[0]-1));
	  passed = passed && (lxm == ! (idbx == 0 && ipbx == 0));
	  passed = passed && (lyp == ! (idby == db3[1]-1 && ipby == pb3[1]-1));
	  passed = passed && (lym == ! (idby == 0 && ipby == 0));
	  passed = passed && (lzp == ! (idbz == db3[2]-1 && ipbz == pb3[2]-1));
	  passed = passed && (lzm == ! (idbz == 0 && ipbz == 0));
	}
      }
    }
	}
      }
    }

    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("neighbor_process");

    passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
    for (int ipbz=0; ipbz<pb3[2]; ipbz++) {
      for (int ipby=0; ipby<pb3[1]; ipby++) {
	for (int ipbx=0; ipbx<pb3[0]; ipbx++) {
	  int ipb = ipbx + pb3[0]*(ipby + pb3[1]*ipbz);
	  bool lxm = layout.neighbor_is_internal(ipb,0,idb,axis_x,-1);
	  bool lxp = layout.neighbor_is_internal(ipb,0,idb,axis_x,+1);
	  bool lym = layout.neighbor_is_internal(ipb,0,idb,axis_y,-1);
	  bool lyp = layout.neighbor_is_internal(ipb,0,idb,axis_y,+1);
	  bool lzm = layout.neighbor_is_internal(ipb,0,idb,axis_z,-1);
	  bool lzp = layout.neighbor_is_internal(ipb,0,idb,axis_z,+1);
	  int ipbxm = layout.neighbor_process(ipb,0,idb,axis_x,-1);
	  int ipbxp = layout.neighbor_process(ipb,0,idb,axis_x,+1);
	  int ipbym = layout.neighbor_process(ipb,0,idb,axis_y,-1);
	  int ipbyp = layout.neighbor_process(ipb,0,idb,axis_y,+1);
	  int ipbzm = layout.neighbor_process(ipb,0,idb,axis_z,-1);
	  int ipbzp = layout.neighbor_process(ipb,0,idb,axis_z,+1);
	  if (lxm) passed = passed && 
		     (((idbx == 0) && (ipb+ipbxm == ipbx - 1)) ||
		      ((idbx != 0) && (ipb+ipbxm == ipbx)));
	  if (lxp) passed = passed && 
		     (((idbx == db3[0]-1) && (ipb+ipbxp == ipbx + 1)) ||
		       ((idbx != db3[0]-1) && (ipb+ipbxp == ipbx)));
	  if (lym) passed = passed && 
		     (((idby == 0) && (ipb+ipbym == ipby - 1)) ||
		       ((idby != 0) && (ipb+ipbym == ipby)));
	  if (lyp) passed = passed && 
		     (((idby == db3[1]-1) && (ipb+ipbyp == ipby + 1)) ||
		       ((idby != db3[1]-1) && (ipb+ipbyp == ipby)));
	  if (lzm) passed = passed && 
		     (((idbz == 0) && (ipb+ipbzm == ipbz - 1)) ||
		       ((idbz != 0) && (ipb+ipbzm == ipbz)));
	  if (lzp) passed = passed && 
		     (((idbz == db3[2]-1) && (ipb+ipbzp == ipbz + 1)) ||
		       ((idbz != db3[2]-1) && (ipb+ipbzp == ipbz)));
	}
      }
    }
	}
      }
    }
    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("neighbor_thread");

    passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
    for (int ipbz=0; ipbz<pb3[2]; ipbz++) {
      for (int ipby=0; ipby<pb3[1]; ipby++) {
	for (int ipbx=0; ipbx<pb3[0]; ipbx++) {
	  int ipb = ipbx + pb3[0]*(ipby + pb3[1]*ipbz);
	  int itxm = layout.neighbor_thread(ipb,0,idb,axis_x,-1);
	  int itxp = layout.neighbor_thread(ipb,0,idb,axis_x,+1);
	  int itym = layout.neighbor_thread(ipb,0,idb,axis_y,-1);
	  int ityp = layout.neighbor_thread(ipb,0,idb,axis_y,+1);
	  int itzm = layout.neighbor_thread(ipb,0,idb,axis_z,-1);
	  int itzp = layout.neighbor_thread(ipb,0,idb,axis_z,+1);
	  passed = passed && (itxm == 0);
	  passed = passed && (itxp == 0);
	  passed = passed && (itym == 0);
	  passed = passed && (ityp == 0);
	  passed = passed && (itzm == 0);
	  passed = passed && (itzp == 0);
	}
      }
    }
	}
      }
    }

    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("extent");
    passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
    for (int ipbz=0; ipbz<pb3[2]; ipbz++) {
      for (int ipby=0; ipby<pb3[1]; ipby++) {
	for (int ipbx=0; ipbx<pb3[0]; ipbx++) {
	  int ipb = ipbx + pb3[0]*(ipby + pb3[1]*ipbz);
	  double lower_extent[3],upper_extent[3];
	  layout.extent(ipb,0,idb,lower_extent,upper_extent);
  	  passed = passed && (fabs(lower_extent[axis_x]-
 				   1.0*(idbx + db3[0]*ipbx)/(pb3[0]*db3[0]))<TOL);
  	  passed = passed && (fabs(lower_extent[axis_y]-
 				   1.0*(idby + db3[1]*ipby)/(pb3[1]*db3[1]))<TOL);
  	  passed = passed && (fabs(lower_extent[axis_z]-
 				   1.0*(idbz + db3[2]*ipbz)/(pb3[2]*db3[2]))<TOL);

  	  passed = passed && (fabs(upper_extent[axis_x]-
 				   1.0*((idbx+1) + db3[0]*ipbx)/(pb3[0]*db3[0]))<TOL);
  	  passed = passed && (fabs(upper_extent[axis_y]-
 				   1.0*((idby+1) + db3[1]*ipby)/(pb3[1]*db3[1]))<TOL);
  	  passed = passed && (fabs(upper_extent[axis_z]-
 				   1.0*((idbz+1) + db3[2]*ipbz)/(pb3[2]*db3[2]))<TOL);
	}
      }
    }
	}
      }
    }
    // NOTE: results sensitive to roundoff
    unit_assert(passed);

    //----------------------------------------------------------------------

    unit_func("array_indices");

    // array size

    int nx = 158;
    int ny = 720;
    int nz = 48;

    passed = true;

    for (int idbz=0; idbz<db3[2]; idbz++) {
      for (int idby=0; idby<db3[1]; idby++) {
	for (int idbx=0; idbx<db3[0]; idbx++) {
	  int idb = idbx + db3[0]*(idby + db3[1]*idbz);
    for (int ipbz=0; ipbz<pb3[2]; ipbz++) {
      for (int ipby=0; ipby<pb3[1]; ipby++) {
	for (int ipbx=0; ipbx<pb3[0]; ipbx++) {
	  int ipb = ipbx + pb3[0]*(ipby + pb3[1]*ipbz);
	  int index_lower[3],index_upper[3];
	  layout.array_indices(ipb,0,idb,nx,ny,nz,index_lower,index_upper);
  	  passed = passed && (index_lower[axis_x] == 
			      nx*(idbx + db3[0]*ipbx)/(pb3[0]*db3[0]));
 	  passed = passed && (index_lower[axis_y] == 
 			      ny*(idby + db3[1]*ipby)/(pb3[1]*db3[1]));
   	  passed = passed && (index_lower[axis_z] == 
 			      nz*(idbz + db3[2]*ipbz)/(pb3[2]*db3[2]));

 	  passed = passed && (index_upper[axis_x] == 
 			      nx*((idbx+1) + db3[0]*ipbx)/(pb3[0]*db3[0]));
   	  passed = passed && (index_upper[axis_y] == 
 			      ny*((idby+1) + db3[1]*ipby)/(pb3[1]*db3[1]));
   	  passed = passed && (index_upper[axis_z] == 
 			      nz*((idbz+1) + db3[2]*ipbz)/(pb3[2]*db3[2]));
	}
      }
    }
	}
      }
    }

    unit_assert (passed);

  }  

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Layout.def.h)
