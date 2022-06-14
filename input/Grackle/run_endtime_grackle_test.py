import argparse
import os.path
import sys

_LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
_TOOLS_DIR = os.path.join(_LOCAL_DIR, "../../tools")
if os.path.isdir(_TOOLS_DIR):
    sys.path.insert(0, _TOOLS_DIR)
    from gen_grackle_testing_file import generate_grackle_input_file
    from run_cpp_test import run_test_and_check_success

else:
    raise RuntimeError(
        f"expected testing utilities to be defined in {_TOOLS_DIR}, but that "
        "directory does not exist"
    )

parent_parser = argparse.ArgumentParser(add_help=False)

parent_parser.add_argument(
    "--grackle-data-file", required = True, type = str,
    help = ("Specifies the path to the grackle data file that is to be used in "
            "the simultaion.")
)
parent_parser.add_argument(
    "--generate-config-path", required = True, type = str,
    help = ("Specifies the path to the configuration file that is generated by "
            "this program. The generated file includes the contents of "
            "--nominal-config-path and overwrites path to the grackle data "
            "file based on --grackle-data-path")
)
parent_parser.add_argument(
    '--launch_cmd', required = True, type = str,
    help = "Specifies the commands used to launch the Enzo-E simulation"
)
parent_parser.add_argument(
    "--output-dump", action = "store", default = None,
    help = ("Specifies the path where a copy of the standard output stream "
            "from the execution of the program should optionally be dumped "
            "(the data is still written to the standard output stream). The "
            "contents of this file may be used to determine the outcome of "
            "the tests.")
)

_description = '''\
Runs a test Grackle-related test that succeeds or fails based on the completion
time of the test. The success/failure of the test is reflected by the return 
code of this program (an exit code of 0 indicates the test was entirely 
successful).
'''

_epilog = '''\
In more detail, this function expects most of the test problem's parameters to 
be specified by the file at the location given by --nominal_config_path.

The program will generate a new configuration file that uses all of the 
parameters from the --nominal-config-path file but overwrites the parameter 
used to specify the grackle data file with the value specified by
--grackle-data-file. The generated config file is written to the path given by
--generate-config-path.

Finally, the program executes Enzo-E with this generated configuration
file and reports whether the tests have passed.
'''

parser = argparse.ArgumentParser(description = _description, epilog = _epilog,
                                 parents = [parent_parser])
parser.add_argument(
    "--nominal-config-path", required = True, type = str,
    help = ("Specifies the path to the configuration file that specifies most "
            "parameters for the test problem.")
)



def run_grackle_test(launcher, nominal_config_path, generate_config_path,
                     grackle_data_file, dump_path = None):
    """
    Runs an enzo-e simulation test problem involving Grackle.

    In detail, this function:
      - expects most of the test problem's parameters to be specified by the 
        file at the location given by `nominal_config_path`.
      - generates a new configuration file at the path given by
        `generate_config_path`. This file includes all of the parameters from
        the `nominal_config_path`, but overwrites the parameter used to specify
        the grackle data file with the value specified by `grackle_data_dir`.
      - the function then executes enzo-e with this generated configuration
        file and reports whether all tests built into the simulation (e.g. an
        expected completion time) have passed, if there are any.

    Parameters
    ----------
    launcher: str
        Specifies the command used to launch enzo-e.
    nominal_config_path: str
        Specifies the path to the config file that specifies the bulk of the 
    generate_config_path: str
        Specifies the path where the temporary input file should be written.
    grackle_data_file: str
        Specifies the path to the grackle data file that is to be used by the
        test problem.
    dump_path: str, optional
        Path to a file where the output of the simulation should be written.
        If this is None (the default), the output is written to a temporary
        file.

    Returns
    -------
    tests_pass: bool
        Specifies whether the simulation ran successfully and whether all tests
        that are built into the simulation have passed (if there are any).
    """

    print("generating config file at {} that uses the grackle data file at {}"\
          .format(generate_config_path, grackle_data_file))
    generate_grackle_input_file(
        include_path = nominal_config_path,
        data_path = grackle_data_file,
        use_abs_paths = True,
        output_fname = generate_config_path
    )

    print("Executing Enzo-E")
    test_passes = run_test_and_check_success(
        command = launcher, args_for_command = [generate_config_path],
        dump_path = "cooling-test.in"
    )
    return test_passes

if __name__ == '__main__':
    args = parser.parse_args()

    test_passes = run_grackle_test(
        launcher = args.launch_cmd,
        nominal_config_path = args.nominal_config_path,
        generate_config_path = args.generate_config_path,
        grackle_data_file = args.grackle_data_file,
        dump_path = args.output_dump
    )
    if test_passes:
        sys.exit(0)
    else:
        sys.exit(1)