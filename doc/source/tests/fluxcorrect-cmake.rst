-----------------
FluxCorrect tests
-----------------

Tests for FluxCorrect method. This currently runs a 3D simulation with
a pure-hydro inclined entropy/contact wave and checks how well the
various fields are conserved.

inclined_contact_smr_ppm
========================
Integrates the entropy/contact wave using the PPM method. This effectively uses static mesh refinement.

inclined_contact_smr_vl
=======================
Integrates the entropy/contact wave using the VL+CT method in hydro mode. This effectively uses static mesh refinement.
