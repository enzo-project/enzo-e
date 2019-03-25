
#include "enzo.hpp"

//======================================================================

EnzoInitialBCenter::EnzoInitialBCenter
(int cycle, double time) throw ()
  : Initial(cycle,time)
{ }

//----------------------------------------------------------------------

void EnzoInitialBCenter::initialize_bfield_center( Block * block )
{
  // Simply sets the values of the cell-centered B-fields based on the
    // previously initialized face-centered B-fields
    EnzoConstrainedTransport ct;
    Grouping bfieldc_group, bfieldi_group;
    bfieldc_group.add("bfield_x","bfield");
    bfieldc_group.add("bfield_y","bfield");
    bfieldc_group.add("bfield_z","bfield");

    bfieldi_group.add("bfieldi_x", "bfield");
    bfieldi_group.add("bfieldi_y", "bfield");
    bfieldi_group.add("bfieldi_z", "bfield");
    ct.compute_center_bfield(block, 0, bfieldc_group, bfieldi_group);
    ct.compute_center_bfield(block, 1, bfieldc_group, bfieldi_group);
    ct.compute_center_bfield(block, 2, bfieldc_group, bfieldi_group);
}

//----------------------------------------------------------------------

void EnzoInitialBCenter::enforce_block( Block * block,
					const Hierarchy * hierarchy ) throw()
{
  EnzoInitialBCenter::initialize_bfield_center(block);
}

  
