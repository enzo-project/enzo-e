#!/bin/python

# runs a generic Grackle test where we compare summary statistics about the
# fields at the final snapshot

import argparse
import os.path
import sys

from run_endtime_grackle_test import run_grackle_test, parent_parser

_LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
_TOOLS_DIR = os.path.join(_LOCAL_DIR, "../../tools")
if os.path.isdir(_TOOLS_DIR):
    sys.path.insert(0, _TOOLS_DIR)
    from field_summary import compare_against_reference

else:
    raise RuntimeError(
        f"expected testing utilities to be defined in {_TOOLS_DIR}, but that "
        "directory does not exist"
    )




parser = argparse.ArgumentParser(
    description = ("Runs a general Grackle test that compares summary "
                   "statistics at a completion time with some reference "
                   "values."),
    parents = [parent_parser]
)

parser.add_argument(
    "--prec", required = True, type = str, choices = ["single", "double"],
    help = "Specifies the precision of Enzo-E."
)

if __name__ == '__main__':
    args = parser.parse_args()

    test_passes = run_grackle_test(
        launcher = args.launch_cmd,
        nominal_config_path = os.path.join(_LOCAL_DIR,
                                           'method_grackle_general.in'),
        generate_config_path = args.generate_config_path,
        grackle_data_file = args.grackle_data_file,
        dump_path = args.output_dump
    )
    # note that there aren't actually any built-in tests in this test problem.
    # Thus, if test_passes is False, that means that Enzo-E crashed
    if not test_passes:
        raise RuntimeError("Enzo-E crashed")


    # now check field values against the reference values
    if args.prec == 'double':
        _ref_tab = os.path.join(_LOCAL_DIR, 'ref_general_grackle-double.csv')
        atol = 0
        # these should be flexible for different compiler versions
        rtol = {"min" : 5e-15, "max" : 5e-6, "mean" : 5e-8,
                "standard_deviation" : 5e-8}
    else:
        _ref_tab = os.path.join(_LOCAL_DIR, 'ref_general_grackle-single.csv')
        atol = 0
        # the following may need to be relaxed for different compiler versions
        rtol = dict((k, 1e-7) for k in ["min","max","mean",
                                        "standard_deviation"])

    test_passes = compare_against_reference(
        './GeneralGrackle-500.00/GeneralGrackle-500.00.block_list',
        ref_summary_file_path = _ref_tab,
        atol = atol, rtol = rtol, report_path = None
    )

    if test_passes:
        sys.exit(0)
    else:
        sys.exit(1)
