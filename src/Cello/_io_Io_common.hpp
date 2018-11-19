  /// Return the ith metadata item associated with the object
  virtual void meta_value 
  (int index, 
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Return the ith data item associated with the object
  virtual void field_array 
  (int index, 
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0,
   int * nx=0,  int * ny=0,  int * nz=0) throw();

   virtual void particle_array 
   (int it, int ib, int ia,
    void ** buffer, std::string * name, int * type,
    int * n, int * k) throw();
