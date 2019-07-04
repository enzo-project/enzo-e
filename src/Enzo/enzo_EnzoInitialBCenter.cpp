/// See LICENSE_CELLO file for license and copyright information
/// @file     enzo_EnzoInitialBCenter.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed Jul 3 2019
/// @brief    [\ref Problem] VLCT Bfield initialization method implementation

#include "enzo.hpp"

//======================================================================

EnzoInitialBCenter::EnzoInitialBCenter
(Parameters * parameters,int cycle, double time) throw ()
  : Initial(cycle,time)
{
  parameters->group_set(0,"Initial");
  parameters->group_set(1,"vlct_bfield");


  // Check if values are specified for any component of the vector potential
  //    - If only a subset of values are specified, than the unspecified
  //      components are assumed to be constant.
  //    - If none of the components are specified, then bfields are assumed to
  //      be pre-calculated and only the cell-centered values are computed
  std::string names[3] = {"Ax","Ay","Az"};
  for (int i = 0; i < 3; i++){
    if (parameters->type(names[i]) == parameter_unknown){
      values_[i] = nullptr;
    } else {
      if ( (parameters->type(names[i]) == parameter_list) &&
	   (parameters->list_length(names[i]) == 0)){
	values_[i] = nullptr;
      }
      values_[i] = new Value(parameters, names[i]);
    }
  }
}

//----------------------------------------------------------------------

void EnzoInitialBCenter::initialize_bfield_center( Block * block )
{
  // Simply sets the values of the cell-centered B-fields based on the
  // previously initialized face-centered B-fields
  EnzoConstrainedTransport ct;
  Grouping bfieldc_group, bfieldi_group;
  bfieldc_group.add("bfield_x","bfield");
  bfieldc_group.add("bfield_y","bfield");
  bfieldc_group.add("bfield_z","bfield");

  bfieldi_group.add("bfieldi_x", "bfield");
  bfieldi_group.add("bfieldi_y", "bfield");
  bfieldi_group.add("bfieldi_z", "bfield");
  ct.compute_center_bfield(block, 0, bfieldc_group, bfieldi_group);
  ct.compute_center_bfield(block, 1, bfieldc_group, bfieldi_group);
  ct.compute_center_bfield(block, 2, bfieldc_group, bfieldi_group);
}

//----------------------------------------------------------------------

void bfieldi_helper_(EFlt3DArray &bfield,
                     CelloArray<double,3> &Aj,
                     CelloArray<double,3> &Ak,
		     int dim, double dj, double dk)
{

  // Aj is centered on edges of i and k dimension but cell-centered along j
  // Ak is centered on edges of i and j dimension but cell-centered along k
  // Aj_right(k,j,i) = Aj(k+1/2,j,i-1/2)    Aj_left(k,j,i) = Aj(k-1/2,j,i-1/2)
  // Ak_right(k,j,i) = Ak(k,j+1/2,i-1/2)    Ak_left(k,j,i) = Ak(k,j-1/2,i-1/2)

  // get dimensions of interface field
  int fc_mx = bfield.shape(2);
  int fc_my = bfield.shape(1);
  int fc_mz = bfield.shape(0);

  CelloArray<double,3> Ak_left, Ak_right,Aj_right, Aj_left;

  if (dim == 0){
    // Aj_right = Ay(iz+1/2,iy,ix-1/2), Ak_right = Az(iz, iy+1/2,ix-1/2)
    Aj_right = Aj.subarray(CSlice(1, fc_mz+1), CSlice(0, fc_my),
                           CSlice(0, fc_mx));
    Ak_right = Ak.subarray(CSlice(0, fc_mz), CSlice(1, fc_my+1),
                           CSlice(0, fc_mx));
  } else if (dim == 1){
    // Aj_right = Az(iz,iy-1/2,ix+1/2), Ak_right = Ax(iz+1/2,iy-1/2,ix)
    Aj_right = Aj.subarray(CSlice(0, fc_mz), CSlice(0, fc_my),
                           CSlice(1, fc_mx+1));
    Ak_right = Ak.subarray(CSlice(1, fc_mz+1), CSlice(0, fc_my),
                           CSlice(0, fc_mx));
  } else {
    // Aj_right = Ax(iz-1/2,iy+1/2,ix), Ak_right = Ay(iz-1/2,iy,ix+1/2)
    Aj_right = Aj.subarray(CSlice(0, fc_mz), CSlice(1, fc_my+1),
                           CSlice(0, fc_mx));
    Ak_right = Ak.subarray(CSlice(0, fc_mz), CSlice(0, fc_my),
                           CSlice(1, fc_mx+1));
  }

  Aj_left = Aj.subarray(CSlice(0, fc_mz), CSlice(0, fc_my), CSlice(0, fc_mx));
  Ak_left = Ak.subarray(CSlice(0, fc_mz), CSlice(0, fc_my), CSlice(0, fc_mx));

  for (int iz=0;iz<fc_mz;iz++){
    for (int iy=0; iy<fc_my; iy++){
      for (int ix=0; ix<fc_mx; ix++){
        // Bi(k,j,i-1/2) =
        //    ( Ak(    k, j+1/2, i-1/2) - Ak(    k, j-1/2, i-1/2) )/dj -
        //    ( Aj(k+1/2,     j, i-1/2) - Aj(k+1/2,     j, i-1/2) )/dk
        bfield(iz,iy,ix) =
	  (enzo_float)(((Ak_right(iz,iy,ix) - Ak_left(iz,iy,ix))/dj -
			(Aj_right(iz,iy,ix) - Aj_left(iz,iy,ix))/dk));
      }
    }
  }
}

//----------------------------------------------------------------------

// Components of the magnetic field from the vector potential                   
// Bx(k, j, i-1/2) =
//    ( Az(    k, j+1/2, i-1/2) - Az(    k, j-1/2, i-1/2) )/dy -
//    ( Ay(k+1/2,     j, i-1/2) - Ay(k+1/2,     j, i-1/2) )/dz
// By(k, j-1/2, i) =
//    ( Ax(k+1/2, j-1/2,     i) - Ax(k+1/2, j-1/2,     i) )/dz -
//    ( Az(    k, j-1/2, i+1/2) - Az(    k, j-1/2, i-1/2) )/dx
// Bz(k-1/2, j, i) =
//    ( Ay(k-1/2,     j, i+1/2) - Ay(k-1/2,     j, i-1/2) )/dx -
//    ( Ax(k-1/2, j+1/2,     i) - Ax(k-1/2, j-1/2,     i) )/dy
void EnzoInitialBCenter::initialize_bfield_interface( Block * block,
						      CelloArray<double,3> &Ax,
						      CelloArray<double,3> &Ay,
						      CelloArray<double,3> &Az)
{

  EFlt3DArray bfieldi_x, bfieldi_y, bfieldi_z;
  EnzoFieldArrayFactory array_factory(block);
  bfieldi_x = array_factory.from_name("bfieldi_x");
  bfieldi_y = array_factory.from_name("bfieldi_y");
  bfieldi_z = array_factory.from_name("bfieldi_z");

  double dx,dy,dz;
  block->data()->field_cell_width(&dx,&dy,&dz);

  bfieldi_helper_(bfieldi_x, Ay, Az, 0, dy,dz);
  bfieldi_helper_(bfieldi_y, Az, Ax, 1, dz,dx);
  bfieldi_helper_(bfieldi_z, Ax, Ay, 2, dx,dy);
}

//----------------------------------------------------------------------

void EnzoInitialBCenter::enforce_block( Block * block,
					const Hierarchy * hierarchy ) throw()
{

  // if specified, optionally initialize the face-centered magnetic field from
  // the vector potential
  if (values_[0] != nullptr || values_[1] != nullptr || values_[2] != nullptr){

    Data* data = block->data();
    Field field = data->field();

    double t = block->time();
    int nx, ny, nz; // number of cells per axis in the active zone
    field.size(&nx, &ny, &nz);
    int gx, gy, gz;
    field.ghost_depth(0, &gx, &gy, &gz);

    int mx, my, mz;
    mx = nx+2*gx;             my = ny+2*gy;             mz = nz+2*gz;

    // allocate corner-centered arrays for the magnetic vector potentials
    // Ax, Ay, and Az are always cell-centered along the x, y, and z dimensions
    double* data_ptrs[3];
    data_ptrs[0] = new double[(mz+1)*(my+1)*(mx)]();
    data_ptrs[1] = new double[(mz+1)*(my)*(mx+1)]();
    data_ptrs[2] = new double[(mz)*(my+1)*(mx+1)]();
    
    double *xc, *yc, *zc, *xf, *yf, *zf;
    xc = new double [mx];     yc = new double [my];     zc = new double [mz];
    xf = new double [mx+1];   yf = new double [my+1];   zf = new double [mz+1];
    data->field_cells (xc, yc, zc, gx, gy, gz);
    data->field_cell_faces (xf, yf, zf, gx, gy, gz);


    // if specified, evaluate the vector potentials
    // (if not, assume they are zero - the values were all initialized to zero)
    if (values_[0] != nullptr){
      values_[0]->evaluate(data_ptrs[0], t,
			   mx  , mx  , xc,
			   my+1, my+1, yf,
			   mz+1, mz+1, zf);
    }
    if (values_[1] != nullptr){
      values_[1]->evaluate(data_ptrs[1], t,
                           mx+1, mx+1, xf,
                           my  , my  , yc,
                           mz+1, mz+1, zf);
    }
    if (values_[0] != nullptr){
      values_[2]->evaluate(data_ptrs[2], t,
			   mx+1, mx+1, xf,
			   my+1, my+1, yf,
			   mz  , mz  , zc);
    }

    // initialize arrays wrapping the data
    CelloArray<double,3> Ax(data_ptrs[0],mz+1,my+1,mx);
    CelloArray<double,3> Ay(data_ptrs[1],mz+1,my,mx+1);
    CelloArray<double,3> Az(data_ptrs[2],mz,my+1,mx+1);

    EnzoInitialBCenter::initialize_bfield_interface(block,Ax,Ay,Az);

    delete[] data_ptrs[0];    delete[] data_ptrs[1];    delete[] data_ptrs[2];
    delete[] xc;              delete[] yc;              delete[] zc;
    delete[] xf;              delete[] yf;              delete[] zf;
  }

  // initialize cell-centered values from initialized face-centered values
  EnzoInitialBCenter::initialize_bfield_center(block);
}

  
