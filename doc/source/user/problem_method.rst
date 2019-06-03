.. _using-methods:

**************
Enzo-E Methods
**************

*[ This page is under development ]*
  
This section decribes all Methods available in Enzo-E: what they do,
what method-specific parameters they have, what fields / particles
they access or update, what solvers if any used, and how different
Methods are intended to be used with each other.

Current Enzo-E methods available include those listed below.

``"ppm"``: hydrodynamics
========================

This implements the modified piecewise parabolic method (PPM) in Enzo.

parameters
----------

.. list-table:: Method ``ppm`` parameters
   :widths: 10 5 1 30
   :header-rows: 1
   
   * - Parameter
     - Type
     - Default
     - Description
   * - ``"density_floor"``
     - `float`
     - `1.0e-6`
     - `Lower limit on density`
   * - ``"diffusion"``
     - `logical`
     - `false`
     - `PPM diffusion parameter`
   * - ``"dual_energy"``
     - `logical`
     - `false`
     - `Whether to use dual-energy formalism`
   * - ``"dual_energy_eta_1"``
     - `float`
     - `0.001`
     - `Dual energy parameter eta 1`
   * - ``"dual_energy_eta_2"``
     - `float`
     - `0.1`
     - `Dual energy parameter eta 2`
   * - ``"flattening"``
     - `integer`
     - `3`
     - `PPM flattening parameter`
   * - ``"minimum_pressure_support_parameter"``
     - `integer`
     - `100`
     - `Enzo's MinimumPressureSupportParameter`
   * - ``"number_density_floor"``
     - `float`
     - `1.0e-6`
     - `Lower limit on number density`
   * - ``"pressure_floor"``
     - `float`
     - `1.0e-6`
     - `Lower limit on pressure`
   * - ``"pressure_free"``
     - `logical`
     - `false`
     - `Pressure-free flag`
   * - ``"steepening"``
     - `logical`
     - `false`
     - `PPM steepening parameter`
   * - ``"temperature_floor"``
     - `float`
     - `1.0e-6`
     - `Lower limit on temperature`
   * - ``"use_minimum_pressure_support"``
     - `logical`
     - `false`
     - `Minimum pressure support`
   * - ``"mol_weight"``
     - `float`
     - `0.6`
     - `Mean molecular mass`

fields
------

.. list-table:: Method ``ppm`` fields
   :widths: 5 5 1 30
   :header-rows: 1

   * - Field
     - Type
     - Read/Write
     - Description   
   * - ``"density"``
     - ``enzo_float``
     - [rw]
     -    
   * - ``"velocity_x"``
     - ``enzo_float``
     - [rw]
     -
   * - ``"velocity_y"``
     - ``enzo_float``
     - [rw]
     - if rank ≥ 2
   * - ``"velocity_z"``
     - ``enzo_float``
     -  [rw]
     - if rank ≥ 3
   * - ``"total_energy"``
     - ``enzo_float``
     - [rw]
     -
   * - ``"acceleration_x"``
     - ``enzo_float``
     - [r]
     - if gravity
   * - ``"acceleration_y"``
     - ``enzo_float``
     - [r]
     - if gravity and rank ≥ 2
   * - ``"acceleration_z"``
     - ``enzo_float``
     - [r]
     - if gravity and rank ≥ 3
   * - ``"pressure"``
     - ``enzo_float``
     - [w]
     - computed from ``total_energy``

``"ppml"``: MHD
===============

PPML ideal MHD solver

``"pm_deposit"``: particle-mesh
===============================

Particle-mesh ("PM") method component to deposit of field and particle
mass into a "total density" field

parameters
----------

.. list-table:: Method ``ppm`` parameters
   :widths: 10 5 1 30
   :header-rows: 1
   
   * - Parameter
     - Type
     - Default
     - Description
   * - ``"alpha"``
     - `float`
     - `0.5`
     - `Deposit mass at time t + alpha * dt`

fields
------

particles
---------
   
``"pm_update"``: particle-mesh
==============================

Particle-mesh ("PM") method component to update particle positions
given acceleration fields
   
``"heat"``: heat equation
=========================

A sample Method for implementing forward-euler to solve the heat equation.   
   
``"grackle"``: chemistry/cooling
================================

Calls methods provided by the external Grackle 3.0 chemistry and
cooling library.
   
``"comoving_expansion"``: comoving expansion
============================================

Adds the comoving expansion terms to the physical variables.
   
``"turbulence"``: driving
=========================

Turbulence driving.

``"gravity"``: particle-mesh
============================

Particle-mesh ("PM") method component to compute gravitational
potential given a total density field, and calculate associated
acceleration fields.
   
``"trace"``: tracer particles
=============================

Moves tracer particles given the velocity field.    
   

