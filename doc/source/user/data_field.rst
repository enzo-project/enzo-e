.. include:: ../roles.incl

*************
Enzo-E Fields
*************

**[ This page is under development ]**

Cello allows field data to be created and operated on in Blocks of an
adaptive mesh hierarchy.  Enzo-E defines specific fields, as well as groups of
related fields.  This page documents what fields are accessed by different Methods,
recommended usage of Fields when writing Methods, and a reference of the Field
API.

-----------------------
Using Fields in Methods
-----------------------

---------
Field API
---------

:Method:    :p:`Field::Field(FieldDescr *, FieldData *)`
:Summary:   :s:`Create a Field object given a Field descriptor and a Field Data object`
:Return:    :t:`none`

:e:`This constructor creates a new Field object given a FieldDescr (field descriptor) and
FieldData (field data) object.`

----
   
:Method:  :p:`Field::Field(FieldDescr \*, FieldData \*)`
:Method:  Field(const Field & field) 
:Summary:   :s:`Copy constructor`
:Return:   :t:`none`
		       
----

:Method:  :p:`Field::operator= (const Field & field)`
:Summary:   :s:`Assignment operator`
:Return:   :t:`Field &`

----

:Method:  :p:`Field::~Field()`
:Summary:   :s:`Destructor`
:Return:   :t:`none`

----

:Method:  :p:`Field::pup (PUP::er &p)`
:Summary:   :s:`CHARM++ Pack / Unpack function`
:Return:   :t:`void`
  
----

:Method:  :p:`Field::field_descr()`
:Summary:   :s:`Return the field descriptor for this field`
:Return:   :t:`FieldDescr \*`

----

:Method:  :p:`Field::field_data()`
:Summary:   :s:`Return the field data for this field`
:Return:   :t:`FieldData \*`

----

:Method:   :p:`Field::set_field_descr(FieldDescr \* field_descr)`
:Summary:  :s:`Set the field descriptor for the Field object`
:Return:   :t:`void`

----

:Method:  :p:`Field::set_field_data(FieldData \* field_data)`
:Summary:  :s:`Set the field data object for the Field object`
:Return:   :t:`void`

  
Field Descriptor
----------------

:Method:  :p:`Field::set_alignment(int alignment)`
:Summary:   :s:`Set alignment`
:Return:   :t:`void`

----

:Method:   :p:`Field::set_padding(int padding)`
:Summary:   :s:`Set padding`
:Return:   :t:`void`

----

:Method:    :p:`Field::set_centering(int id, int cx, int cy=0, int cz=0)`
:Summary:   :s:`Set centering for a field`
:Return:   :t:`void`

----

:Method:   :p:`Field::set_ghost_depth(int id, int gx, int gy=0, int gz=0)`
:Summary:   :s:`Set ghost_depth for a field`
:Return:   :t:`void`
    
----

:Method:   :p:`Field::set_precision(int id, int precision)`
:Summary:   :s:`Set precision for a field`
:Return:   :t:`void`
    
----

:Method:   :p:`Field::insert_permanent(const std::string & name)`
:Summary:   :s:`Insert a new field`
:Return:   :t:`int`

----

:Method:   :p:`Field::insert_temporary(const std::string & name = "")`
:Summary:   :s:`Insert a new field`
:Return:   :t:`int`

----

:Method:   :p:`Field::field_count() const`
:Summary:   :s:`Return the number of fields`
:Return:   :t:`int`

----

:Method:  :p:`Field::field_name(int id) const`
:Summary:   :s:`Return name of the ith field`
:Return:   :t:`std::string`

----

:Method:  :p:`Field::is_field(const std::string & name) const`
:Summary:   :s:`Return whether the field has been inserted`
:Return:   :t:`bool`

----

:Method:  :p:`Field::field_id(const std::string & name) const`
:Summary:   :s:`Return the integer handle for the named field`
:Return:   :t:`int`

Properties
----------

:Method:  :p:`Field::groups()`
:Summary:  :s:`Return the grouping object for Fields`
:Return:   :t:`Grouping \*`

----

:Method:  :p:`alignment() const`
:Summary:   :s:`alignment in bytes of fields in memory`
:Return:   :t:`int`

----

:Method:   :p:`padding() const`
:Summary:   :s:`padding in bytes between fields in memory`
:Return:   :t:`int`

----

:Method:   :p:`centering(int id, int \* cx, int \* cy = 0, int \* cz = 0) const`
:Summary:   :s:`centering of given field`
:Return:   :t:`void`
    
----

:Method:   :p:`is_centered(int id) const`
:Summary:   :s:`return whether the field variable is centered in the cell`
:Return:   :t:`bool`

----

:Method:   :p:`ghost_depth(int id, int \* gx, int \* gy = 0, int \* gz = 0) const`
:Summary:   :s:`depth of ghost zones of given field`
:Return:   :t:`void`
    
----

:Method:   :p:`precision(int id) const`
:Summary:   :s:`Return precision of given field`
:Return:   :t:`int`

----

:Method:   :p:`bytes_per_element(int id) const`
:Summary:   :s:`Number of bytes per element required by the given field`
:Return:   :t:`int`

----

:Method:   :p:`is_permanent (int id_field) const`
:Summary:   :s:`Whether the field is permanent or temporary`
:Return:   :t:`bool`

----

:Method:   :p:`num_permanent() const`
:Summary:   :s:`Return the number of permanent fields`
:Return:   :t:`int`


History
-------

:Method:   :p:`set_history (int num_history)`
:Summary:   :s:`Set the history depth for storing old field values`
:Return:   :t:`void`

----

:Method:   :p:`num_history () const`
:Summary:   :s:`Return the number of history generations to store`
:Return:   :t:`int`

----
  
:Method:   :p:`save_history (double time)`
:Summary:   :s:`Copy "current" fields to history`
:Return:   :t:`void`

:e:`Copy "current" fields to history = 1 fields (saving time), and
push back older generations up to num_history()`
  
----

:Method:   :p:`history_time (int ih) const`
:Summary:   :s:`Return time for given history`
:Return:   :t:`double`

Units
-----

:Method:   :p:`units_scale_cgs (int id, double amount)`
:Summary:   :s:`scale the field to cgs units given the unit scaling factor`
:Return:   :t:`void`

:e:`if it's already in cgs, then leave as-is
except if it's in cgs but the scaling factor has changed (e.g. due to
expansion) then adjust for the new scaling factor`
    
----

:Method:   :p:`units_scale_code (int id, double amount)`
:Summary:   :s:`convert the field to "code units" given the unit scaling factor`
:Return:   :t:`void`

:e:`if it's already in code units, leave it as-is warning if scaling factor has changed.`
	  
----

:Method:   :p:`units_scaling (const FieldDescr \*, int id)`
:Summary:   :s:`Return the current scaling factor of the given Field`
:Return:   :t:`double`

:e:`1.0 if in code units, or the scaling factor if in cgs`

FieldData
---------

:Method:   :p:`size(int \* nx, int \* ny = 0, int \* nz = 0) const`
:Summary:   :s:`Return size of fields on the data, assuming centered`
:Return:   :t:`void`

----

:Method:   :p:`dimensions(int id_field,int \* mx, int \* my = 0, int \* mz = 0) const`
:Summary:   :s:`Return dimensions of fields on the data, assuming centered`
:Return:   :t:`void`

----

:Method:   :p:`values (int id_field, int index_history=0)`
:Method:   :p:`values (std::string name, int index_history=0)`
:Summary:   :s:`Return full array of values for the corresponding field`
:Return:   :t:`char \*`

:e:`Return array for the corresponding field, which may or may not contain ghosts depending on if they're allocated`

----

:Method:   :p:`unknowns (int id_field, int index_history=0)`
:Method:  :p:`unknowns (std::string name, int index_history=0)`
:Summary:   :s:`Return array for the corresponding field`
:Return:   :t:`char \*`

:e:`Return array for the corresponding field, which does not contain ghosts whether they're allocated or not`
	     
----

:Method:   :p:`permanent ()  const`
:Summary:  :s:`Return the array of all permanent fields`
:Return:   :t:`const char \*`

----

:Method:  :p:`cell_width(double xm,   double xp,   double \* hx,
 		  double ym=0, double yp=0, double \* hy=0,
		  double zm=0, double zp=0, double \* hz=0) const`
:Summary:   :s:`Return width of cells along each dimension`
:Return:   :t:`void`

----

:Method:   :p:`clear (float value = 0.0, int id_first = -1, int id_last  = -1)`
:Summary:   :s:`Clear specified array(s) to specified value`
:Return:   :t:`void`
 
----

:Method:   :p:`permanent_allocated() const`
:Summary:   :s:`Return whether array is allocated or not`
:Return:   :t:`bool`

----

:Method:   :p:`permanent_size() const`
:Summary:   :s:`Return whether array is allocated or not`
:Return:   :t:`size_t`

----

:Method:   :p:`allocate_permanent(bool ghosts_allocated = false)`
:Summary:   :s:`Allocate storage for the field data`
:Return:   :t:`void`

----

:Method:   :p:`allocate_temporary(int id)`
:Summary:   :s:`Allocate storage for the temporary fields`
:Return:   :t:`void`

----

:Method:   :p:`deallocate_temporary(int id)`
:Summary:   :s:`Deallocate storage for the temporary fields`
:Return:   :t:`void`

----

:Method:   :p:`reallocate_permanent(bool ghosts_allocated = false)`
:Summary:   :s:`Reallocate storage for the field data`
:Return:   :t:`void`

:e:`Reallocate storage for the field data, e.g. when changing from ghosts to non-ghosts [ costly for large blocks ]`

----

:Method:   :p:`deallocate_permanent()`
:Summary:   :s:`Deallocate storage for the field data`
:Return:   :t:`void`

----

:Method:   :p:`ghosts_allocated() const`
:Summary:   :s:`Return whether ghost cells are allocated or not.`
:Return:   :t:`bool`

----

:Method:   :p:`field_size (int id, int \*nx=0, int \*ny=0, int \*nz=0) const`
:Summary:   :s:`Return the number of elements (nx,ny,nz) along each axis`
:Return:   :t:`int`

:e:`Return the number of elements (nx,ny,nz) along each axis, and total number of bytes n`

Debugging
---------  

:Method:   :p:`print (const char \* message, bool use_file = false) const`
:Summary:   :s:`Print basic field characteristics for debugging`
:Return:   :t:`void`

-------------
Enzo-E Fields
-------------

* acceleration_x
* acceleration_y
* acceleration_z
* B
* bfieldx
* bfieldx_rx
* bfieldx_ry
* bfieldx_rz
* bfieldy
* bfieldy_rx
* bfieldy_ry
* bfieldy_rz
* bfieldz
* bfieldz_rx
* bfieldz_ry
* bfieldz_rz
* cooling_time
* density
* density_total
* dens_rx
* dens_ry
* dens_rz
* DI_density
* DII_density
* driving_x
* driving_y
* driving_z
* e_density
* gamma
* H2I_density
* H2II_density
* HDI_density
* HeI_density
* HeII_density
* HeIII_density
* HI_density
* HII_density
* HM_density
* internal_energy
* metal_density
* potential
* pressure
* species_De
* species_DI
* species_DII
* species_H2I
* species_H2II
* species_HDI
* species_HeI
* species_HeII
* species_HeIII
* species_HI
* species_HII
* species_HM
* temperature
* total_energy
* velocity_x
* velocity_y
* velocity_z
* velox
* velox_rx
* velox_ry
* velox_rz
* veloy
* veloy_rx
* veloy_ry
* veloy_rz
* veloz
* veloz_rx
* veloz_ry
* veloz_rz
* X
