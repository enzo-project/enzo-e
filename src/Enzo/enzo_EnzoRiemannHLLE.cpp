#include "cello.hpp"
#include "enzo.hpp"

// Some Notes
// Stone+ (08) has a typo in equation 52. It should be q_R-q_L instead of
// (q_i-q_{i-q}).


//----------------------------------------------------------------------

void EnzoRiemannHLLE::solve (Block *block, std::vector<int> &priml_ids,
			     std::vector<int> &primr_ids,
			     std::vector<int> &flux_ids, int dim,
			     EnzoEquationOfState *eos)
{
  // Do stuff!
  const int length = primr_ids->size();
  // The following may not be legal
  enzo_float wl[length], wr[length], Ul[length], Ur[length];
  enzo_float Fl[length], Fr[length];
  enzo_float bp, bm;

  // Declare the block and get field object
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // Load in the fields
  std::vector<enzo_float*> wl_arrays;
  std::vector<enzo_float*> wr_arrays;
  std::vector<enzo_float*> flux_arrays;
  for (int i=0;i<length;i++){
    wl_arrays.push_back((enzo_float *) field.values(priml_ids[i]));
    wr_arrays.push_back((enzo_float *) field.values(primr_ids[i]));
    flux_arrays.push_back((enzo_float *) field.values(flux_ids[i]));
  }

  // get integration limits 
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];
  
  // these values are cell-centered, they need to be updated to face-centered
  if (dim == 0) {
    mx++;
  } else if (dim == 1) {
    my++;
  } else {
    mz++;
  }

  // For PLM we only care about fluxes for the third cell in
  // For Nearest-Neighbor, we care about the second cell in
  for (int iz=1; iz<mx-1; iz++) {
    for (int iy=1; iy<my-1; iy++) {
      for (int ix=1; ix<mz-1; ix++) {
	// compute the index
	int i = ix + mx*(iy + my*iz);

	// first, handle the left interfaces

	// get the primitive fields
	for (int field_ind=0; field_ind<length; field_ind++){
	  wl[field_ind] = wl_arrays[field_ind][i];
	}

	// compute the conserved quantity
	eos->conservative_from_primitive(wl,Ul);

	// compute the interface flux
	interface_flux_ (wl, Ul, Fl, dim);

	// Next, handle the left interfaces
	for (int field_ind=0; field_ind<length; field_ind++){
	  wr[field_ind] = wr_arrays[field_ind][i];
	}

	// compute the conserved quantity
	eos->conservative_from_primitive(wr,Ur);

	// compute the interface flux
	interface_flux_ (wr, Ur, Fr, dim);

	// now compute the wavespeeds
	wave_speeds_ (wl, wr, eos, dim, &bp, &bm);

	// Now compute the Riemann Flux

	for (int field_ind=0; field_ind<length; field_ind++){
	  if (field_ind == 5+dim){
	    flux_arrays[field_ind][i] = 0;
	  }
	  // fill in the following line properly
	  flux_arrays[field_ind][i] = (((bp*Fl[field_ind] - bm*Fr[field_ind])
				        / (bp - bm)) +
				       ((Ul[field_ind] - Ur[field_ind])*bp*bm
					/(bp - bm)));

	}
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoRiemannHLLE::wave_speeds_ (enzo_float *wl, enzo_float *wr,
				    EnzoEquationOfState *eos, int dim,
				    enzo_float *bp, enzo_float *bm)
{
  // do stuff

}



void EnzoRiemannHLLE::interface_flux_ (enzo_float *prim, enzo_float* cons,
				       enzo_float *fluxes,
				       EnzoEquationOfState *eos, int dim)
{
  // This assumes that MHD is included
  // This may be better handled by the EquationOfState
  enzo_float rho, vx, vy, vz, p, Bx, By, Bz, Etot, B2, Bv, vi, Bi;
  rho = prim[0];
  vx = prim[1];
  vy = prim[2];
  vz = prim[3];
  p  = prim[4];
  Bx = prim[5];
  By = prim[6];
  Bz = prim[7];

  Etot = cons[4];

  B2 = Bx*Bx + By*By + Bz*Bz;
  Bv = Bx*vx + By*vy + Bz*vz;

  // Compute Fluxes
  fluxes[0] = cons[1+dim];

  
  if (dim == 0){
    vi = vx;
    Bi = Bx;
    // Fluxes for Mx, My, Mz
    fluxes[1] = cons[1+dim]*vx - Bx*Bi+ p + B2/2;
    fluxes[2] = cons[1+dim]*vy - By*Bi;
    fluxes[3] = cons[1+dim]*vz - Bz*Bi;
    fluxes[5] = 0;
    fluxes[6] = By*vi - Bi*vy;
    fluxes[7] = Bz*vi - Bi*vz;
  } else if (dim == 1){
    vi = vy;
    Bi = By;
    // Fluxes for Mx, My, Mz
    fluxes[1] = cons[1+dim]*vx - Bx*Bi;
    fluxes[2] = cons[1+dim]*vy - By*Bi + p + B2/2;
    fluxes[3] = cons[1+dim]*vz - Bz*Bi;
    // Fluxes for Bx,By,Bz
    fluxes[5] = Bx*vi - Bi*vx;
    fluxes[6] = 0;
    fluxes[7] = Bz*vi - Bi*vz;
  } else {
    vi = vz;
    Bi = Bz;
    // Fluxes for Mx, My, Mz
    fluxes[1] = cons[1+dim]*vx - Bx*Bi;
    fluxes[2] = cons[1+dim]*vy - By*Bi + p + B2/2;
    fluxes[3] = cons[1+dim]*vz - Bz*Bi;
    // Fluxes for Bx,By,Bz
    fluxes[5] = Bx*vi - Bi*vx;
    fluxes[6] = By*vi - Bi*vy;
    fluxes[7] = 0;
  }

  fluxes[4] = (cons[4] + p + B2/2)*vi - Bv*Bi;
  
}
