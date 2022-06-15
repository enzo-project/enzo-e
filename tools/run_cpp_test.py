import argparse
import subprocess
import sys
import tempfile

_description = '''\
Runs a test that is directly encoded in a C++ binary (the binary is specified
by COMMAND and optional arguments are specified by ARGS), that makes use of the
Cello testing machinery. After the test completes, this program determines
whether it was succesful and returns an appropriate exit code (an exit code of
0 indicates the test was entirely successful).
'''


parser = argparse.ArgumentParser(description = _description)
parser.add_argument(
    "--output-dump", action = "store", default = None,
    help = ("Specifies the path where a copy of the standard output stream "
            "from the execution of the program should optionally be dumped "
            "(the data is still written to the standard output stream). The "
            "contents of this file may be used to determine the outcome of "
            "the tests.")
)
parser.add_argument(
    "command", metavar = "COMMAND", action = "store",
    help = "the C++ binary that is to be executed"
)
parser.add_argument(
    "args_for_command", metavar = "ARGS", action = "store",
    nargs = argparse.REMAINDER, default = [],
    help = "the arguments to the C++ binary that are to be executed"
)

def execute_command(command, args, output_dump):
    """
    Executes the command. This will write a copy of stdout to the file 
    specified by output_dump (while continuing to write to stdout).

    Returns
    -------
    success: bool
        True indicates the program completed succesfully
    """
    arg_l = [command] + args + ['|', 'tee', output_dump]
    exit_code = subprocess.call(' '.join(arg_l), shell = True)
    return exit_code == 0

def log_suggests_test_failure(log_path):
    """
    Returns whether the log has any contents suggestive of a test failure.

    Notes
    -----
    This logic is ported from the build.sh script from earlier versions of
    Enzo-E (for concreteness, we considered the file from commit 
    3f7f33f4c6254f718dc3deead5c8119a50f26318)
    """

    def _num_occurences(pattern, skip_binary_file_search = False):
        # employs grep to count the number of occurences of a pattern
        if skip_binary_file_search:
            flags = ""
        else:
            flags = "-I"
        command = 'grep {} "{}" {} | wc -l'.format(flags, pattern, log_path)
        return int(subprocess.check_output(command, shell = True))

    # check the log for any lines explicitly reporting failures (lines 184-189)
    #  - I suspect that passing the -I flag to grep in the original shell
    #    script was a typo (which tells grep to assume that a binary file
    #    doesn't contain a match) 
    n_fail_lines = _num_occurences("^ FAIL", skip_binary_file_search = True)
    #n_incomplete_lines = _num_occurences("^ incomplete",
    #                                     skip_binary_file_search = True)
    #n_pass_lines = _num_occurences("^ pass", skip_binary_file_search = True)
    if n_fail_lines > 0:
        return True

    # check for signs that a test may have crashed (lines 199-202)
    num_crashes = (_num_occurences("UNIT TEST BEGIN") -
                   _num_occurences("UNIT TEST END"))
    if num_crashes != 0:
        return True

    # check other indications of a failure (lines 242-243)
    # - the following logic was originally used to determine failures in a
    #   nicely formated text report but did not affect the exit code.
    # - the original logic suggests that there was a problem with the test if
    #   _num_occurences("BEGIN") == 0 (but again didn't reflect that in the
    #   exit code)
    has_BEGIN = _num_occurences("BEGIN") > 0
    has_END_CELLO = _num_occurences("END CELLO") > 0
    if has_BEGIN and not has_END_CELLO:
        return True # the test failed
    else:
        return False

def run_test_and_check_success(command, args_for_command = [],
                               dump_path = None):
    if dump_path is None:
        delete_output = True
        dump_path = tempfile.mktemp() # path for a temporary file
    else:
        delete_output = False

    command_success = execute_command(command = command,
                                      args = args_for_command,
                                      output_dump = dump_path)
    success = command_success and not log_suggests_test_failure(dump_path)

    # cleanup the temporary file
    if delete_output:
        os.remove(dump_path)
    return success

if __name__ == '__main__':
    args = parser.parse_args()

    test_passes = run_test_and_check_success(
        command = args.command, args_for_command = args.args_for_command,
        dump_path = args.output_dump
    )

    if test_passes:
        sys.exit(0)
    else:
        sys.exit(1)
