#ifndef MESH_BLOCK_HPP
#define MESH_BLOCK_HPP

class Block {

public:

  Block(int ix, int iy, int iz,
	int nx, int ny, int nz,
	double xm, double ym, double zm,
	double hx, double hy, double hz,
	int num_field_blocks) throw();

  Block() {printf ("%s:%d Oops\n",__FILE__,__LINE__);};

  Block (CkMigrateMessage *m) {printf ("%s:%d Oops\n",__FILE__,__LINE__);};

  virtual ~Block();

protected:

  int index_[3];
  double lower_[3];
  double upper_[3];

};

#endif /* MESH_BLOCK_HPP */

