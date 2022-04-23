// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputImage.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Io] Declaration for the OutputImage class

#ifndef IO_OUTPUT_IMAGE_HPP
#define IO_OUTPUT_IMAGE_HPP

class Factory;

enum mesh_color_type {
  mesh_color_unknown,
  mesh_color_level,
  mesh_color_process,
  mesh_color_neighbor,
  mesh_color_age,
  mesh_color_order
};

class OutputImage : public Output {

  /// @class    OutputImage
  /// @ingroup  Io
  /// @brief [\ref Io] class for writing images

public: // functions

  /// Empty constructor for Charm++ pup()
  OutputImage() throw() {}

  /// Create an uninitialized OutputImage object
  OutputImage(int index,
	      const Factory * factory,
	      int process_count,
	      const int root_size[3],
	      const int root_blocks[3],
	      int min_level, int max_level,
              int leaf_only,
	      std::string image_type,
	      int         image_size[2],
	      std::string image_reduce_type,
	      std::string image_mesh_color,
	      std::string image_mesh_order,
	      std::string image_color_particle_attribute,
	      double      image_lower[],
	      double      image_upper[],
	      int face_rank,
	      int axis,
	      bool image_log,
	      bool image_abs,
	      bool include_ghost,
              bool use_min_max,
	      double min_value, double max_value) throw();

  /// OutputImage destructor: free allocated image data
  virtual ~OutputImage() throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(OutputImage);

  /// Charm++ PUP::able migration constructor
  OutputImage (CkMigrateMessage *m)
    : Output (m),
      image_data_(NULL),
      image_mesh_(NULL),
      op_reduce_(reduce_unknown),
      mesh_color_type_(mesh_color_unknown),
      mesh_color_order_(),
      color_particle_attribute_(""),
      axis_(axis_all),
      min_value_(std::numeric_limits<double>::max()),
      max_value_(-std::numeric_limits<double>::max()),
      png_(NULL),
      image_type_(""),
      face_rank_(0),
      image_log_(false),
      image_abs_(false),
      include_ghost_(false),
      min_level_(0),
      max_level_(0),
      leaf_only_(false)
  {
    colormap_[0].clear();
    colormap_[1].clear();
    colormap_[2].clear();
    for (int axis=0; axis<3; axis++) {
      image_lower_[axis] = -std::numeric_limits<double>::max();
      image_upper_[axis] =  std::numeric_limits<double>::max();
    }
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  // Set the image colormap
  void set_colormap (std::vector<float> colormap[3]);

public: // virtual functions

  /// Prepare for accumulating block data
  virtual void init () throw();

  /// Open (or create) a file for IO
  virtual void open () throw();

  /// Close file for IO
  virtual void close () throw();

  /// Cleanup after output
  virtual void finalize () throw();

  /// Write block-related field and particle data
  virtual void write_block ( const Block * block ) throw();

  /// Write fields
  virtual void write_field_data
  ( const FieldData * field_data, 
    int index_field) throw();

  /// Write particles
  virtual void write_particle_data
  ( const ParticleData * particle_data, 
    int index_particle) throw();

  /// Prepare local array with data to be sent to remote chare for processing
  virtual void prepare_remote (int * n, char ** buffer) throw();

  /// Accumulate and write data sent from a remote processes
  virtual void update_remote  ( int n, char * buffer) throw();

  /// Free local array if allocated; NOP if not
  virtual void cleanup_remote (int * n, char ** buffer) throw();

private: // functions

  /// value associated with the given mesh level
  double mesh_color_(const Block *, int level ) const;

  bool type_is_mesh_ () const
  { return (image_type_ == "mesh" || image_type_ == "data+mesh"); }

  bool type_is_data_ () const
  { return (image_type_ == "data" || image_type_ == "data+mesh"); }

  bool is_active_ (const Block * block) const;

  /// Create the png file object
  void png_create_ (std::string filename) throw();

  /// Delete the png object
  void png_close_() throw();

  /// Create the image data object
  void image_create_ () throw();

  /// Generate PNG image, using given min and max for colormap
  void image_write_ () throw();

  /// Close the image data
  void image_close_ () throw();

   /// Generate a PNG image of array data
  void reduce_point_
  ( double * data,  int ix, int iy, double value, double alpha=1.0) throw();

  void reduce_line_(double * data, int ixm, int ixp, int iym, int iyp, 
		    double value, double alpha=1.0);
  void reduce_line_x_(double * data, int ixm, int ixp, int iy, 
		      double value, double alpha=1.0);
  void reduce_line_y_(double * data, int ix, int iym, int iyp, 
		      double value, double alpha=1.0);
  void reduce_box_(double * data, int ixm, int ixp, int iym, int iyp, 
		   double value, reduce_type reduce, double alpha=1.0);
  void reduce_box_filled_(double * data, int ixm, int ixp, int iym, int iyp, 
		    double value, double alpha=1.0);

  double data_(int i) const ;

private: // attributes

  /// Color map
  std::vector<float> colormap_[3];

  /// Current image for data
  double * image_data_;

  /// Current image for mesh
  double * image_mesh_;

  /// Reduction operation
  reduce_type op_reduce_;

  /// Color type for mesh images
  int mesh_color_type_;

  /// Block int Scalar for ordering when mesh_color_type_ = mesh_color_order
  std::string mesh_color_order_;

  /// Particle attribute defining color (default -1: constant)
  std::string color_particle_attribute_;

  /// Axis along which to reduce
  axis_type axis_;

  /// Minimum and maximum values if specified
  bool use_min_max_;
  double min_value_;
  double max_value_;

  /// Current image size in pixels
  int image_size_[2];
  
  /// Current pngwriter
  pngwriter * png_;

  /// Image type: data or mesh
  std::string image_type_;

  /// Minimal rank of faces to include face level indicators 
  int face_rank_;

  /// Whether to plot the log of the field
  int image_log_;

  /// Whether to plot the absolute value of the field
  int image_abs_;

  /// Whether to include ghost zones
  bool include_ghost_;

  /// Maximum mesh level of Block to output
  int min_level_;

  /// Maximum mesh level of Block to output
  int max_level_;

  /// Whether to restrict Blocks to only leaf nodes
  bool leaf_only_;

  /// Lower and upper bounds on image (can be used for slices)
  double image_lower_[3];
  double image_upper_[3];
  
};

#endif /* IO_OUTPUT_IMAGE_HPP */
