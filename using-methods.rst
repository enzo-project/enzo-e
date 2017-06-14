**************
Enzo-P methods
**************

*[ This page is under development ]*
  
This section decribes all Methods available in Enzo-P: what they do,
what fields / particles they access or update, solvers used, and how
different Methods are intended to be used with each other.

Current Enzo-P methods available include those listed below.

``"ppm"``
   PPM hydrodynamics solver

``"ppml"``
   PPML ideal MHD solver

``"pm_deposit"``
   Particle-mesh ("PM") method component to deposit of field and
   particle mass into a "total density" field
   
``"pm_update"``
   Particle-mesh ("PM") method component to update particle positions given acceleration fields
   
``"heat"``
   A sample Method for implementing forward-euler to solve the heat equation.   
   
``"grackle"``
   Calls methods provided by the external Grackle 3.0 chemistry and cooling library.
   
``"turbulence"``
   Turbulence driving.

``"gravity"``
   Particle-mesh ("PM") method component to compute gravitational potential given a total density field, and calculate associated acceleration fields.
   
``"trace"``
   Moves tracer particles given the velocity field.    
   

