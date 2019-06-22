****************
Using CelloArray
****************

==========
Background
==========

Historically, to access elements at a given location, ``(ix, iy, iz)`` of
an array in Enzo, the developer would need to explicitly calculate the
index of the pointer using knowledge of the underlying shape of the
array represented by the pointer.

To simplify and enhance the readability of code in Enzo-E, we have
implemented ``CelloArray``, a Multi-dimensional Array class template
that wraps the data and encapsulate array operations. This class
template draw’s loose inspiration from Athena++’s ``AthenaArray`` and
numpy’s ``ndarray``.

See the first two cases listed in :ref:`array-examples`
for comparisons of snipets written using ``CelloArray`` and traditional
pointer operations. These examples reflect operations performed in
Enzo-E.


Throughout the Enzo portion of the codebase, we extensively use the type
``EFlt3DArray`` which acts is an alias for ``CelloArray<enzo_float,3>``.

============
Design Goals
============

The design of the class was primarily driven by the following specifications:

  * Emphasize fast access to array elements by passing the index along each
    dimension to the ``operator()`` method.

    * This method can be inlined within for-loops and for 3D arrays it has
      the same complexity as that of ``AthenaArray``.

    * Simple benchmarks show that the current implementation achieves
      performance comparable to c-style array access

  * The CelloArray needs tp be able to allocate and manage its own memory AND
    wrap existing pointers (namely the pointers allocated by the Cello’s Field
    framework)

    * This allows code using the ``CelloArray`` to coexist alongside code which
      use pointers in a more conventional way.

  * ``CelloArray`` needs be able to represent a view of a (mostly contiguous)
    subarray of a pre-existing instance of ``CelloArray``. This facilitates the
    encapsulation of a directional mesh operation in a single generalized
    function (e.g. writing a single flux function for all directions rather
    than separate functions to compute flux along the x, y, and z directions).


==============
User Interface
==============

The class template is formally defined as ``CelloArray<T,D>`` where
``T`` is the contained type (frequently ``enzo_float``) and ``D`` is
the number of dimensionsions of the array. The interface is inspired
by numpy’s ``ndarray`` .

It is more straightforward to describe the user interface for
``CelloArray<T,D>`` by describing different operations with examples
rather than providing a detailed API.

Array Creation
--------------
Simplest initialization:

  * Use the constructor ``CelloArray(Args... args)`` to construct an
    array of 0s of shape ``(arg0, arg1, ... arg{D-1})``.  The resulting
    array owns the underlying memory and deallocation is entirely
    taken care of

  * Examples:

    * Construct an array of shape ``(2,3,4)`` that holds doubles:

    .. code-block:: c++

       CelloArray<double,3> arr(2,3,4);


    * Construct an array of shape ``(5,)`` that holds ints:
       
    .. code-block:: c++

       CelloArray<int,1> arr(5); 
       CelloArray<int,1> arr2 = CelloArray<int,1>(5); // yields same result

Wrap a pre-existing pointer:

  * Use the constructor ``CelloArray(T* array, Args... args);`` to wrap the
    pointer ``array`` which represents an array with shape
    ``(arg0, arg1, ... arg{D-1})``.

  * Example: Construct an array representing ``[[0,1,2],[3,4,5]]``:

  .. code-block:: c++

     int data[] = {0,1,2,3,4,5};
     CelloArray<int,2> arr(data,2,3);

We can also forward declare an array and assign values to it later.

.. code-block:: c++

   int data[] = {0,1,2,3,4,5};
   CelloArray<int,1> arr; 
   arr = CelloArray<int,2>(data,2,3);


Dimension Size
--------------

To get the length along a dimension (or axis), call
``arr.shape(unsigned int dim)``, where ``dim`` is the number of the
dimension. Dimensions numbers start at ``0`` and are ordered with
increasing indexing speed (``dim=D-1`` is the dimension with fastest
indexing).

Element Access
--------------

To access an element pass indices to the ``operator()(Args... args)``
method. As many indices should be specified as there are dimensions in
the array (the number of args **must** match the number of dimensions.

The ``operator()(Args... args)`` method returns a reference or copy
(depending on the circumstance) of the element.

**Example:** print element ``(0,2)`` of the array ``[[0,1,2],[3,4,5]]``:

.. code-block:: c++

   int data[] = {0,1,2,3,4,5};
   CelloArray<int,2> arr(data,2,3);
   printf("%d\n", arr(0,2)); // prints "2"
   // printf("%d\n", arr(2));       This would fail to compile
   // printf("%d\n", arr(0,0,2));   This would fail to compile


Simple Assignment - Shallow/Deep Copies
---------------------------------------

Shallow copies are produced via ordinary assignment.
.. code-block:: c++

   int data[] = {0,1,2,3,4,5};
   CelloArray<int,2> a(data,2,3);
   CelloArray<int,2> b = a; // b is now a shallow copy of arr
   CelloArray<int,2> c(2,2); // c represents [[0,0],[0,0]]
   CelloArray<int,2> d = c; // d is now a shallow copy of c
   c = a; // c is now a shallow copy of a

When ``c`` is assinged the contents of ``a``, ``c`` becomes a shallow
copy of ``a``. However the contents of ``d`` are unaffected.  It still
represents the array ``[[0,0],[0,0]]``.

To perform a deepcopy, assign the the results of the ``deepcopy`` method.

.. code-block:: c++

   int data[] = {0,1,2,3,4,5};
   CelloArray<int,2> a(data,2,3);
   CelloArray<int,2> b = a; // b is now a shallow copy of arr
   CelloArray<int,2> c(2,2); // c represents [[0,0],[0,0]]
   CelloArray<int,2> d = c; // d is now a shallow copy of c
   c = a; // c is now a shallow copy of a
   
Modifications to the contents of ``e`` will not be reflected in ``a``
or ``data`` (and vice-versa)


Creating Subarrays
------------------
Calling ``arr.subarray(Args... args)`` returns a (mostly contiguous) view
of a subarray specified by ``args``, where ``args`` represent the slices
along each dimension. Each ``arg`` should be an instance of ``CSlice`` and
the number of ``args`` **must** match the number of dimensions of the array.
Calling ``arr.subarray()`` without any arguments returns a shallow copy

``CSlice`` is a class that represents the start and stop points
along a given dimension. The constructor standard is simply:
``CSlice(int start, int stop)``.

Subarray Examples
~~~~~~~~~~~~~~~~~

We present an extanded example below. We start by defining a subarray,
``sub`` of an array ``arr`` (which wraps an existing pointer of data
and represents the array ``[[0,1,2],[3,4,5]]``).

.. code-block:: c++

   int data[] = {0,1,2,3,4,5};
   CelloArray<int,2> arr(data,2,3);
   CelloArray<int,2> sub = arr.subarray(CSlice(0,2),CSlice(1,3));
   printf("%d\n", sub(1,0)) // prints "4";

At this point ``sub`` represents the subarray ``[[1,2],[4,5]]``
of the full array held by ``arr``. ``sub`` is truly a "view" of
``arr``. Modifications to the elements of ``sub`` and
modifications to elements in ``arr`` (if it lies in the subarray),
are reflected in both locations.

.. code-block:: c++

   arr(1,3) *= -3;
   sub(0,0) = -100;

After executing the above block of code, ``arr`` now represents
``[[0,-100,2],[3,4,-15]]`` and ``sub`` represents the subarray
``[[-100,2],[4,-15]]``.

``CelloArray`` also provides support for taking subarrays of
subarrays (or taking subarrays of shallow copies). If we define
a subarray of ``sub`` the result will represent a view of the
same underlying data

.. code-block:: c++

   CelloArray<int,2> sub_of_sub = sub.subarray(CSlice(0,2),CSlice(0,1));
   sub_of_sub(1,0) +=8;

After the above operations, ``arr`` now reflects the full array
``[[0,-100,2],[3,12,-15]]``, while ``sub`` and ``sub_of_sub``
represent the subarrays ``[[-100,2],[12,-15]]`` and ``[[-100],[12]]``.
Continuing to make shallow copies or subarrays of ``sub_of_sub`` and
its derivatives will still yield views of the original array.

If we assign ``arr`` the value of an unrelated array, the data
tracked by all subarrays and subcopies are unaffected.

.. code-block:: c++

   CelloArray<int,2> sub2 = arr.subarray(CSlice(1,2),CSlice(0,3));
   arr = CelloArray<int, 2>(3,3); // setting arr equal to another array
   sub(1,0) /= -2;

After execution of the preceeding block of code, ``sub`` represents
``[[-100,2],[-6,-15]]`` of the full array,
``sub_of_sub`` represents ``[[-100],[-6]]``, and ``sub2`` represents
``[[3,-6,-15]]`` (at this point the ``data`` pointer holds
``[0, -100, 2, 3, -6, -15]``).

The fact that ``arr`` originally wrapped ``data`` has no bearing on
the outcomes described above for each instance of ``CelloArray``.
We illustrate this below with an analogous abreviated example, where
the analog to ``arr``, called ``array``, originally owns its data.

.. code-block:: c++

   CelloArray<int,2> array(2,3);
   array(0,0) = 0;    array(0,1) = 1;    array(0,2) = 2;
   array(1,0) = 3;    array(1,1) = 4;    array(1,2) = 5;
   CelloArray<int,2> subarray = array.subarray(CSlice(0,2), CSlice(1,3));
   array(1,3) *= -3;
   subarray(0,0) = -100;
   CelloArray<int,2> subarray_of_subarray = subarray.subarray(CSlice(0,2),
                                                              CSlice(0,1));
   subarray_of_subarray(1,0) += 8;

After executing the preceeding block of code, ``array`` reflects
``[[0,-100,2],[3,12,-15]]``, while ``subarray`` and
``subarray_of_subarray`` represent the subarrays
``[[-100,2],[12,-15]]`` and ``[[-100],[12]]``. If this was all the
code we executed, the memory of ``array`` would be freed after its
destructor and the destructors of all of subarrays or shallowcopies
are called.

If we reassign ``array`` to a different array, just like before, the values
of its subarrays and shallow copies will be unaffected.

.. code-block:: c++

   CelloArray<int,2> subarray2 = array.subarray(CSlice(1,2),CSlice(0,3));
   array = CelloArray<int, 2>(3,3);
   subarray(1,0) /= -2;

Now, ``subarray`` represents ``[[-100,2],[-6,-15]]`` from the full
array, ``subarray_of_subarray`` represents ``[[-100],[-6]]``, and
``subarray2`` represents ``[[3,-6,-15]]``. We note that no memory
has been deallocated. The memory will only be deallocated after
``subarray``, ``subarray_of_subarray``, and ``subarray2`` have
all had their deconstructor called and/or been assigned unrelated
arrays, assuming no additional subarrays or shallowcopies of any of
the 3 variables are made in the meantime (in that case the memory
would still not be deallocated until any additional
subarrays/shallowcopies that view the original data are destroyed).

Additional CSlice features
~~~~~~~~~~~~~~~~~~~~~~~~~~
``CSlice`` provides two additional features to simplify code when
the generating subarrays of a ``CelloArray`` instance. These are

  1. The constructor supports negative indexing. For example
     ``CSlice(1,-1)`` represents a slice starting at the second
     element and stopping at (does not include) the last element
     along a dimension. Additionally, ``CSlice(-3,-1)`` represents
     starting from the third-to-last and stopping at the last
     element along a given dimension.
  2. The constructor accepts the ``NULL`` and ``nullptr`` as the
     ``stop`` argument and understands it to mean that the last element
     along the axis. For example, ``CSlice(1, NULL)`` and
     ``CSlice(1,nullptr)`` both represent slices from the second
     element through the last element of the dimension.
     ``CSlice(-3,NULL)`` and ``CSlice(-3,nullptr)`` both represent
     slices extending from the third-to-last element through the last
     element of a dimension. Additionally,  if ``NULL`` or ``nullptr``
     are passed as the ``start`` argument, they are understood to mean
     that the slice starts at the first element
     (``CSlice(0,NULL)``, ``CSlice(0,nullptr)``, ``CSlice(NULL,NULL)``, &
     ``CSlice(nullptr,nullptr)`` are all equivalent). 

Finally, we note that ``CSlice`` provides a default constructor to
simplify the construction of arrays of slices. However, to help avoid
bugs, we require that any default-constructed ``CSlice`` must be
assigned a non-default constructed value (or an error will be raised).


Elementwise Assignment
----------------------

We also provide elementwise assignment (copying elements between arrays).
To invoke this, ``arr.subarray(Arg... args)`` must appear on the LHS
(left-hand side) of an ``=``. The expression on the RHS (right-hand) side
can be either:

  * A scalar of the type contained by ``arr``. In this case, all elements in
    the resulting subarray are set equal to the scalar.

  * Another array or subarray that contains the same type of elements as
    ``arr`` (The dimensions & shape of the RHS array must match the LHS
    subarray). In this case, the elements in the LHS subarray are each
    set equal to the corresponding elements of the array on the RHS.

An example is illustrated below:

.. code-block:: c++

   int data[] = {0,1,2,3,4,5,6,7,8,9,10,11};
   CelloArray<int,2> arr(data,3,4);
   // arr reflects: [[0,1,2,3],[4,5,6,7],[8,9,10,11]]
   CelloArray<int,2> arr2(2,2); // arr2 is initially [[0,0],[0,0]]
   arr2.subarray(CSlice(0,2),
                 CSlice(0,2)) = 7; //arr2 is now [[7,7],[7,7]]
   // The previous 2 lines could be re-written as: arr2.subarray() = 7;
   arr.subarray(CSlice(1,3), CSlice(0,2)) = arr2;
   // arr now reflects: [[0,1,2,3],[7,7,6,7],[7,7,10,11]]
   arr.subarray(CSlice(0,3), 
                CSlice(1,3)) = arr.subarray(CSlice(0,3), CSlice(2,4));
   // arr now reflects: [[0,2,3,3],[7,6,7,7],[7,10,11,11]]
   CelloArray<int,2> arr3 = arr.subarray(CSlice(1,3), CSlice(2,4));
   // arr3 is a view of the subarray: [[7,7],[11,11]] of arr
   //arr3 = 17;   // This will not compile
   arr3 = arr2; // arr3 is now a shallow copy of arr2 & arr is unchanged


===========
Convenience
===========

In the Enzo layer of the codebase, we provide several short-cuts for
performing frequent actions related to the ``CelloArray`` to reduce
boilerplate code.

  * We define and make extensive use of the type ``EFlt3DArray`` which
    is an alias for ``CelloArray<enzo_float,3>``.

  * We define the class ``EnzoFieldArrayFactory`` which drastically
    reduces the boilerplate code associated with the initialization of
    instances of ``CelloArray`` that wrap Cello fields.

  * We define the class ``EnzoPermutedCoordinates`` convenience class
    which helps reduce boilerplate code associated with writing
    functions using instances of ``CelloArray`` that are generalized
    with respect to dimension.

Two additional, features that can be enabled at compile-time to assist
with debugging by defining macros before the inclusion of the ``CelloArray``
header file.

  * Defining the ``CHECK_BOUNDS`` macro, will cause checks of the validity of
    indices every time an element is accessed and will raise an error when it
    detects that an element that lies outside of the array bounds.

  * Defining the ``CHECK_FINITE_ELEMENTS`` macro will cause a check during
    retrieval of array elements that they are not ``NaN`` or ``inf``

.. _array-examples:

========
Examples
========

Below, we show some factored out, simplified examples, ways in which how
``CelloArray`` might simplify code:

Copying Elements
----------------

This example illustrates how ``CelloArray`` simplifies the code
required to copy elements between arrays. (We illustrate how one might
write Nearest Neighbor reconstruction along the x-direction).

This code assumes a mesh with shape ``(mz, my, mx)``. Suppose we have:

  * An ``(mz,my,mx)`` array of cell-centered primitives ``w``

  * An ``(mz,my,mx-1)`` array of left reconstructed values, ``wl``

  * An ``(mz,my,mx-1)`` array of right reconstructed values, ``wr`` 

First is an the ``CelloArray`` version:
    
.. code-block:: c++

   typedef double enzo_float;
   typedef CelloArray<enzo_float,3> EFlt3DArray;

   void reconstruct_NN_x(EFlt3DArray &w, EFlt3DArray &wl, 
                         EFlt3DArray &wr){
       wl.subarray() = w.subarray(CSlice(0,w.shape(0)),
                                  CSlice(0,w.shape(1)),
                                  CSlice(0,-1));
       wr.subarray() = w.subarray(CSlice(0,w.shape(0)),
                                  CSlice(0,w.shape(1)),
                                  CSlice(1,w.shape(2)));
   }

The analogous code using conventional pointer operations is:

.. code-block:: c++

   typedef double enzo_float;

   void reconstruct_NN_x(enzo_float *w, enzo_float *wl, enzo_float *wr,
                      int mx, int my, int mz){
     int offset = 1;
     for (int iz=0; iz<mz-1; iz++) {
       for (int iy=0; iy<my-1; iy++) {
         for(int ix=0; ix<mx-1; ix++) {
           int i = (iz*my + iy)*mx + ix;
           int i_xf = (iz*my + iy)*(mx-1) + ix; 
           wl[i_xf] = w[i];
           Wr[i_xf] = w[i + offset];
         }
       }
     }
   }

Adding Flux Divergence
----------------------

We show a factored out, slightly simplified version of the code used
to add the flux divergence in an unsplit manner. This example is one
of the more notable cases where the ``CelloArray`` leads to more
transparent code.

This code assumes a mesh with shape (mz, my, mx). Suppose we have:

  * An ``(mz,my,mx)`` array of cell-centered conserved quantities ``u``

  * An ``(mz,my,mx-1)`` array of x-face centered fluxes in the x-direction,
    ``xflux``

  * An ``(mz,my-1,mx)`` array of y-face centered fluxes in the y-direction,
    ``yflux``

  * An ``(mz-1,my,mx)`` array of y-face centered fluxes in the z-direction,
    ``zflux``

  * The timestep is ``dt``, and the size of cells along the x, y, and z
    directions are ``dx``, ``dy``, ``dz``

  * We set place the updated values in ``out`` (which may be a
    reference to the same array as ``u`` or to a different array)

.. code-block:: c++

   typedef double enzo_float;
   typedef CelloArray<enzo_float,3> EFlt3DArray;

   void  update_cons(EFlt3DArray &u, EFlt3DArray &out,
                     EFlt3DArray &xflux, EFlt3DArray &yflux,
                     EFlt3DArray &zflux, enzo_float dt, enzo_float dx,
                     enzo_float dy, enzo_float dz){
     enzo_float dtdx = dt/dx;
     enzo_float dtdy = dt/dy;
     enzo_float dtdz = dt/dz;

     for (int iz=1; iz<u.shape(0)-1; iz++) {
       for (int iy=1; iy<u.shape(1)-1; iy++) {
         for (int ix=1; ix<u.shape(2)-1; ix++) {
           out(iz,iy,ix) = (u(iz,iy,ix) -
                            dtdx*(xflux(iz,iy,ix) - xflux(iz,iy,ix-1)) -
                            dtdy*(yflux(iz,iy,ix) - yflux(iz,iy-1,ix)) -
                            dtdz*(zflux(iz,iy,ix) - zflux(iz-1,iy,ix)));
         }
       }
     }
   }

The analogous function using conventional pointer operations is provided below:

.. code-block:: c++

   typedef double enzo_float;
   typedef CelloArray<enzo_float,3> EFlt3DArray;

   void update_cons(enzo_float *u, enzo_float *out, 
                    enzo_float *xflux, enzo_float *yflux,
                    enzo_float *zflux, enzo_float dt, 
                    enzo_float dx, enzo_float dy, enzo_float dz,
                    int mx, int my, int mz){
     enzo_float dtdx = dt/dx;
     enzo_float dtdy = dt/dy;
     enzo_float dtdz = dt/dz;

     int x_offset = 1;
     int y_offset = mx;
     int z_offset = my*mx;

     for (int iz=1; iz<mz-1; iz++) {
       for (int iy=1; iy<my-1; iy++) {
         for (int ix=1; ix<mx-1; ix++) {
           int i = (iz*my + iy)*mx + ix;
           int i_zf = i;
           int i_yf = (iz*(my-1) + iy) * mx + ix;
           int i_xf = (iz*my + iy) * (mx-1) + ix;

           out[i] = (u[i] 
                     - dtdx * (xflux[i_xf] - xflux[i_xf - x_offset])
                     - dtdy * (yflux[i_yf] - yflux[i_yf - y_offset])
                     - dtdz * (zflux[i_zf] - zflux[i_zf - z_offset]));
         }
       }
     }
   }



Direction Generalized Functions
-------------------------------

This example illustrates how subarrays allows functions using
``CelloArray`` to be written so that they are generalized with respect
to Cartesian direction. Due to the simplicity of the example, code
with conventional pointer operations is comparable to the code using
arrays (however arrays make more complex examples more understandable)

In the van Leer + Constrained Transport scheme, we need to update
update the cell-centered B-field component along a given direction by
averaging the same components of the B-field stored at cell
interfaces. We track Bx at the x-faces, By at the y-faces and Bz at
the z-faces.

This code assumes a mesh with shape ``(mz, my, mx)``. Suppose we have:

  * An array of cell-centered B-field values (along a given component ) ``bc``

  * An array of interface B-field values (for the same component) ``bi``.
    This array includes values of cell faces on the exterior of the mesh (e.g.
    for values centered along the x-axis the shape would be ``(mz,my,mx+1)``).

  * The direction of the component of the B-field is passed in with ``dim``.
    The values 0,1 & 2 map to x, y, and z

.. code-block:: c++

   typedef double enzo_float;
   typedef CelloArray<enzo_float,3> EFlt3DArray;

   void calc_center_bfield(EFlt3DArray &bc, EFlt3DArray &bi, int dim){
     EFlt3DArray bi_l = bi;

     // The following is a repeating pattern that gets factored out into 
     // a helper function
     EFlt3DArrau bi_r;
     if (dim == 0) {
       bi_r = bi.subarray(CSlice(0,NULL), CSlice(0,NULL), CSlice(1,NULL));
     } else if (dim == 1) {
       bi_r = bi.subarray(CSlice(0,NULL), CSlice(1,NULL), CSlice(0,NULL));
     } else {
       bi_r = bi.subarray(CSlice(0,NULL), CSlice(1,NULL), CSlice(0,NULL));
     }

     for (int iz=0; iz<bc.shape(0); iz++) {
       for (int iy=0; iy<bc.shape(1); iy++) {
         for(int ix=0; ix<bc.shape(2); ix++) {
           bc(iz,iy,ix) = 0.5 * (bi_l(iz,iy,ix) + bi_r(iz,iy,ix));
         }
       }
     }
   }
