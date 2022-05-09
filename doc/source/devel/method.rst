***************
Writing Methods
***************

*[ This page is under development ]*

All method classes are subclasses of the ``Method`` class. The virtual
methods that concrete classes override are shown below:

.. doxygenfunction:: Method::compute

.. doxygenfunction:: Method::name

.. doxygenfunction:: Method::timestep

.. doxygenfunction:: Method::compute_resume

.. warning::

   This page is very incomplete. Among other things, we have not discussed
   ``pup`` routines.
