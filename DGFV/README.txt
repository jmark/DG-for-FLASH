Implementation of the discontinuous Galerkin solver in FLASH documented in

  Markert, Johannes, Stefanie Walch, and Gregor Gassner.
  "A Discontinuous Galerkin Solver in the FLASH Multi-Physics Framework."
  arXiv preprint arXiv:2112.11318 (2021). url: https://arxiv.org/abs/2112.11318

DISCLAIMER: This module is still considered experimental and may give unexpected results!

For questions or feedback, please write an email to

  johannes.markert@jmark.de

!! --------------------- !!

This is the Discontinuous-Galerkin-Finite-Volume (DGFV) blending scheme
for (ideal) magneto-hydrodynamics.

Per design this module implements a high-order *unsplit* DG solver which may
lead to more accurate results at the price of more timesteps.

The module is supposed to be fully compatible to the interface provided by the
*split* solver unit without any changes to the FLASH framework or your setup
configuration. There is one exception!

When you are running (real) astrophysical simulations with extreme initial
conditions (in a numerical sense) such as densities at scales of say 1e-24
[g/cm^3] or spatial dimensions of 1e19 [pc] than the DGFV solver gives
erroneous results. To alleviate this problem you have to define characteristic
scales in your 'flash.par'. This allows the solver to internally transform the
physical states to numerically sensible scales. For example:

  dgfv_characteristic_time        = 1e14
  dgfv_characteristic_density     = 1e-24
  dgfv_characteristic_dimension   = 1e19

!! --------------------- !!

The DGFV module allows to run in three modes. In DG-only, FV-only and the
blending scheme. For example, if you want to run in FV-only mode then define
the following in your 'flash.par':

  dgfv_fv_active = .true.
  dgfv_dg_active = .false.
  dgfv_runge_kutta = "RALSTON"

If you want to run in DG-only method then define
the following in your 'flash.par':

  dgfv_fv_active = .false.
  dgfv_dg_active = .true.
  dgfv_runge_kutta = "LOW_STORAGE_45"

For more options please consult the 'Config' file of this module.

!! --------------------- !!

For installation unpack 'DGFV.zip' and copy the folder 'DGFV' to:

    ${FLASH}/source/physics/Hydro/HydroMain/split

Use the DGFV solver by requiring the unit

    physics/Hydro/HydroMain/split/DGFV

in your setup.
