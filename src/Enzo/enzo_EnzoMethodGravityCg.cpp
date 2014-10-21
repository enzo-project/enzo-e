// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityCg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoMethodGravityCg class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodGravityCg::EnzoMethodGravityCg (int iter_max, double res_tol) 
  : Method(), iter_max_(iter_max), res_tol_(res_tol)
{
}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | iter_max_;
  p | res_tol_;

}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::compute ( CommBlock * comm_block) throw()
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);
  Field field = enzo_block->block()->field();

  // density field

  const int id = field.field_id("density");
  int gd3[3];
  field.ghosts(id,&gd3[0],&gd3[1],&gd3[2]);
  int nd3[3];
  field.size(&nd3[0],&nd3[1],&nd3[2]);
  int md3[3] = {nd3[0] > 1 ? nd3[0]+2*gd3[0] : 1,
		nd3[1] > 1 ? nd3[1]+2*gd3[1] : 1,
		nd3[2] > 1 ? nd3[2]+2*gd3[2] : 1};
  void * density = field.values(id);

  // gravitational potential field
  const int ip = field.field_id("potential");
  int gp3[3];
  field.ghosts(ip,&gp3[0],&gp3[1],&gp3[2]);
  int np3[3];
  field.size(&np3[0],&np3[1],&np3[2]);
  int mp3[3] = {np3[0] > 1 ? np3[0]+2*gp3[0] : 1,
		np3[1] > 1 ? np3[1]+2*gp3[1] : 1,
		np3[2] > 1 ? np3[2]+2*gp3[2] : 1};
  void * potential = field.values(ip);

  // precision
  const int p = field.precision(id);

  if      (p == precision_single)
    compute_  ((float*)  density,   md3,nd3,
	       (float*)  potential, mp3,np3);
  else if (p == precision_double)
    compute_ ((double*)  density,   md3,nd3,
	      (double*)  potential, mp3,np3);
  else if (p == precision_quadruple)
    compute_ ((long double*)  density,   md3,nd3,
	      (long double*)  potential, mp3,np3);
  else 
    ERROR1("EnzoMethodGravityCg()", "precision %d not recognized", p);
}

//======================================================================

template <class T>
void EnzoMethodGravityCg::compute_ 
(T * density,   int md3[3], int nd3[3],
 T * potential, int mp3[3], int np3[3]) const throw()
// function [x] = conjgrad(A,b,x)
//     r=b-A*x;
//     p=r;
//     rsold=r'*r;
//  
//     for i=1:1e6
//         Ap=A*p;
//         alpha=rsold/(p'*Ap);
//         x=x+alpha*p;
//         r=r-alpha*Ap;
//         rsnew=r'*r;
//         if sqrt(rsnew)<1e-10
//               break;
//         end
//         p=r+rsnew/rsold*p;
//         rsold=rsnew;
//     end
// end
// nabla ^ 2 (potential) = 4 pi G density
// 
//
// B = 4 * pi * G * density
// X = X0
// REFRESH(X)
// R = A*X
// R = b - R
// P = R
// rs0 = DOT(R,R)
// while (not converged) {
//    REFRESH(P)
//    AP = A*P
//    pap = DOT(P,AP)
//    alpha = rs0/pap
//    X = X * alpha*P
//    R = R - alpha*AP
//    rs = DOT(R,R);
//    if (rs < res_tol) break;
//    P = R + rs/rs0*P
//    rs0 = rs
//  }
//  potential = x
{

  printf ("potential size %d %d %d  dim %d %d %d\n",
	  np3[0],np3[1],np3[2],
	  mp3[0],mp3[1],mp3[2]);

  printf ("density size %d %d %d  dim %d %d %d\n",
	  nd3[0],nd3[1],nd3[2],
	  md3[0],md3[1],md3[2]);

}
