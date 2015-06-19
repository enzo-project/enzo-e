// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputImage.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Io] Declaration for the OutputImage class

#ifndef IO_OUTPUT_IMAGE_HPP
#define IO_OUTPUT_IMAGE_HPP

class Factory;
class FieldDescr;


enum mesh_color_type {
  mesh_color_unknown,
  mesh_color_level,
  mesh_color_process,
  mesh_color_neighbor,
  mesh_color_age
};

class OutputImage : public Output {

  /// @class    OutputImage
  /// @ingroup  Io
  /// @brief [\ref Io] class for writing images

public: // functions

  /// Empty constructor for Charm++ pup()
  OutputImage() throw() {}

  /// Create an uninitialized OutputImage object
  OutputImage(const FieldDescr * field_descr,
	      int index,
	      const Factory * factory,
	      int process_count,
	      int nx0, int ny0, int nz0,
	      int nxb, int nyb, int nzb,
	      int max_level,
	      std::string image_type,
	      int image_size_x, int image_size_y,
	      std::string image_reduce_type,
	      std::string image_mesh_color,
	      int         image_block_size,
	      int face_rank,
	      bool image_log,
	      bool ghost,
	      bool specify_bounds,
	      double min, double max) throw();

  /// OutputImage destructor: free allocated image data
  virtual ~OutputImage() throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(OutputImage);

  /// Charm++ PUP::able migration constructor
  OutputImage (CkMigrateMessage *m) : Output (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  // Set the image colormap
  void set_colormap
  (int n, double * map_r, double * map_g, double * map_b)
  throw();

  // Set the axis for projecting
  void set_axis (axis_type axis) throw()
  { axis_ = axis; };

public: // virtual functions

  /// Prepare for accumulating block data
  virtual void init () throw();

  /// Open (or create) a file for IO
  virtual void open () throw();

  /// Close file for IO
  virtual void close () throw();

  /// Cleanup after output
  virtual void finalize () throw();

  /// Write block-related field data
  virtual void write_block
  ( const Block * block,
    const FieldDescr * field_descr) throw();

  /// Write fields
  virtual void write_field_data
  ( const FieldData * field_data, 
    const FieldDescr * field_descr,
    int field_index) throw();

  /// Prepare local array with data to be sent to remote chare for processing
  virtual void prepare_remote (int * n, char ** buffer) throw();

  /// Accumulate and write data sent from a remote processes
  virtual void update_remote  ( int n, char * buffer) throw();

  /// Free local array if allocated; NOP if not
  virtual void cleanup_remote (int * n, char ** buffer) throw();

private: // functions

  /// value associated with the given mesh level
  double mesh_color_(int level, int age) const;

  bool type_is_mesh () const
  { return (image_type_ == "mesh" || image_type_ == "data+mesh"); }

  bool type_is_data () const
  { return (image_type_ == "data" || image_type_ == "data+mesh"); }

  /// Create the png file object
  void png_create_ (std::string filename) throw();

  /// Delete the png object
  void png_close_() throw();

  /// Create the image data object
  void image_create_ () throw();

  /// Generate PNG image, using given min and max for colormap
  void image_write_ (double min=0.0, double max=0.0) throw();

  /// Close the image data
  void image_close_ () throw();

   /// Generate a PNG image of array data
  void reduce_point_ ( double * data, 
		       double value, double alpha=1.0) throw();

  void extents_img_ (const Block * block,
		     int *ixm, int *ixp,
		     int *iym, int *iyp,
		     int *izm, int *izp ) const;

  void reduce_line_(double * data, int ixm, int ixp, int iym, int iyp, 
		    double value, double alpha=1.0);
  void reduce_line_x_(double * data, int ixm, int ixp, int iy, 
		      double value, double alpha=1.0);
  void reduce_line_y_(double * data, int ix, int iym, int iyp, 
		      double value, double alpha=1.0);
  void reduce_box_(double * data, int ixm, int ixp, int iym, int iyp, 
		   double value, reduce_type reduce, double alpha=1.0);
  void reduce_cube_(double * data, int ixm, int ixp, int iym, int iyp, 
		    double value, double alpha=1.0);

  double data_(int i) const ;

private: // attributes

  /// Color map
  std::vector<double> map_r_;
  std::vector<double> map_g_;
  std::vector<double> map_b_;

  /// Current image for data
  double * image_data_;

  /// Current image for mesh
  double * image_mesh_;

  /// Reduction operation
  reduce_type op_reduce_;

  /// Color
  int mesh_color_type_;

  /// Axis along which to reduce
  axis_type axis_;

  /// Whether to use given or computed min/max for colormap
  bool specify_bounds_;

  /// Minimum and maximum values if specified
  double min_;
  double max_;

  /// Current image size (depending on axis_)
  int nxi_,nyi_,nzi_;

  /// Current pngwriter
  pngwriter * png_;

  /// Image type: data or mesh
  std::string image_type_;

  /// Minimal rank of faces to include face level indicators 
  int face_rank_;

  /// Whether to plot the log of the field
  int image_log_;

  /// Whether to include ghost zones
  bool ghost_;

  /// Maximum mesh level
  int max_level_;
};

#endif /* IO_OUTPUT_IMAGE_HPP */
