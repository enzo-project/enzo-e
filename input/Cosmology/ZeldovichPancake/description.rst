.. _zeldovich-pancake-test:

------------------
Zel'dovich Pancake
------------------

The Zel'dovich Pancake `(Zel'dovich 1970) <https://ui.adsabs.harvard.edu/abs/1970A%2526A.....5...84Z>`__ test problem is a standard 1D test of cosmology.

* In his original work, Zel'dovich applied linear perturbation theory to solve the evolution of a pressure-free fluid (i.e. dark matter) in an isotropic universe, when a sinusoidal perturbation.
  He arrived at an analytic solution that is exact during the linear phase, up until the point where the caustic form (this is when the pancake forms).
  You can think of the formation of the caustic as a 1D analog for collapse into filaments and clusters.

* We primarily consider the description of the solution presented in `Anninos et al. (1994) <https://adsabs.harvard.edu/abs/1994ApJ...436...11A>`__.
  The density and velocity profiles of a pressure-free fluid at some arbitrary redshift :math:`z_0`, that occurs before :math:`z_c` (the redshift of collapse or caustic formation) are given by

  .. math::
     :nowrap:

     \begin{eqnarray}
       \rho(x_l) &=& \rho_{\rm bkg} \left[1 - \frac{1+z_c}{1+z_0}\cos(k x_l)\right]^{-1}, \\
       v(x_l) &=& -H_0 \frac{1+z_c}{\sqrt{1+z_0}}\, \frac{\sin(k x_l)}{k}.
     \end{eqnarray}

  In these equations :math:`H_0` is the Hubble constant's value at :math:`z=0` and :math:`k=2\pi/\lambda` is the wave-number of the perturbation.
  These equations are functions of :math:`x_l`, the Lagrangian position of the fluid elements or the positions of the fluid elements at the redshift when the perturbation is imposed.
  When you consider the position on the grid at :math:`z_0`, you are considering a Eulerian coordinate, :math:`x_e`.
  The conversion between these two quantities is:

  .. math::

     x_e = x_l - \frac{1 + z_c}{1+z_0} \frac{sin(k x_l)}{k}.

In the test, we start a simulation at some redshift :math:`z_{\rm 0}`, and initialize an ordinary ideal gas such that it obeys the above density and velocity profiles (of the pressure-free fluid).
The internal energy of this gas is initialized so that entropy is constant everywhere.

This problem is useful because it checks the minimal set of physics required for a cosmological simulation (self-gravity, hydrodynamics, expansion).
It also serves as a test of the dual energy formalism.

The actual test we run is described as follows:

* We consider a flat universe with :math:`\Omega_b=1` with :math:`h=0.5` in a narrow box that extends over :math:`64\, {\rm Mpc}\, h^{-1}`.
  We initialize the simulation at :math:`z_0 = 20`, with :math:`z_c = 1`, :math:`\rho_0 = \rho_c` (the critical density), and a background temperature of :math:`100 {\rm K}`.

* While this test technically is just 1D, the VL integrator can only run in 3D.
  For that reason, we consider a long narrow box that is thin along the other axes (the side of each cell is constant).
  We effectively run multiple identical versions of this problem in parallel.

**Note on timing:** The ppm integrator will run this test in fewer cycles because of differences in the minimum courant factor.


.. COMMENT BLOCK

   We may wish to introduce other variants of this problem. For example we could modify the conditions to introduce Bfields.

     * This is similar to how `Li et al. (2008) <https://adsabs.harvard.edu/abs/2008ApJS..174....1L>`__ and `Collins et al. (2010) <https://adsabs.harvard.edu/abs/2010ApJS..186..308C>`__ extend the conditions from  `Ryu et al. (1993) <https://adsabs.harvard.edu/abs/1993ApJ...414....1R>`__

     * I believe the initializer from Enzo shows us how to do this

   We may also want to run a version of the problem with AMR, and a version with drift velocity.

   May also want to try a dual-pancake version... (for a more




