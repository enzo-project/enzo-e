#ifndef INDEX3_HPP
#define INDEX3_HPP

class Index3 {

public:

  Index3();

  Index3(int iax, int iay, int iaz);

  bool operator == (const Index3 & index) const;

  bool operator != (const Index3 & index) const;

  inline int operator [] (std::size_t i) const
  { return v_[i]; }

  inline int & operator [] (std::size_t i)
  { return v_[i]; }

  /// Set the Index3 according to raw bit values
  inline void set_values (int ix, int iy, int iz)
  {
    v_[0] = ix;
    v_[1] = iy;
    v_[2] = iz;
  }

  /// Return the packed bit index for the given axis
  inline void values (int v3[3]) const
  { v3[0] = v_[0];
    v3[1] = v_[1];
    v3[2] = v_[2];
  }
  /// Comparison operator required for Charm++ pup()
  friend bool operator < (const Index3 & x, const Index3 & y) {
    if (x.v_[2] < y.v_[2]) return true;
    if (x.v_[2] > y.v_[2]) return false;
    // else x.v_[2] == y.v_[2]
    if (x.v_[1] < y.v_[1]) return true;
    if (x.v_[1] > y.v_[1]) return false;
    // else x.v_[1] == y.v_[1]
    return  (x.v_[0] < y.v_[0]);
  }

  ///--------------------
  /// PACKING / UNPACKING
  ///--------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  char * load_data (char * buffer);

private: // methods

private: // attributes

   int v_[3];
};

#ifndef TEST
  PUPbytes(Index3)

  class CkArrayIndexIndex3 : public CkArrayIndex {
    Index3 * index_;
  public:
    CkArrayIndexIndex3(const Index3 &in)
    {
      index_ = new (index) Index3(in);
      nInts=sizeof(Index3)/sizeof(int);
    }
  };
#endif

#endif /* INDEX3_HPP */
