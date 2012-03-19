  /// Return the ith metadata item associated with the object
  virtual void meta_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Return the ith data item associated with the object
  virtual void data_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * nxd=0, int * nyd=0, int * nzd=0,
   int * nx=0,  int * ny=0,  int * nz=0) throw();

