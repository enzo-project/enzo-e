#!/usr/bin/python
# this currently works with python 2 or 3

import argparse
import io
import os
import os.path
import shutil
import subprocess
import sys
import tempfile

from cello_parse import load, dump, isstr
from test_report import create_test_report

def _subprocess_exit_success(args):
    """
    Launches a subprocess and indicates whether the return code is 0. (This is
    a pretty hacky - we could probably do better)
    """
    try:
        # subprocess.check_output will hide stdout from the user
        subprocess.check_output(args)
        return True
    except subprocess.CalledProcessError: # exit code is nonzero
        return False

def _command_line_utility_exists(name):
    """ Returns whether a given command line utility exists """
    return _subprocess_exit_success(['which',name])

class EnzoE:
    """
    Wraps the Enzo-E binary. This can be used to execute Enzo-E from an 
    arbitrary directory.

    Parameters
    ----------
    enzoe_binary: string
        Path to the enzo-e binary.
    charmrun_binary: string, optional
        Path to the charmrun binary. If this is not specified, the simulation 
        can't be run in parallel.
    parallel_args: str or list, optional
        List of parallel arguments that are passed to charmrun.
    """
    def __init__(self, enzoe_binary, charmrun_binary = None,
                 parallel_args = None,):
        self._program_args = []

        if charmrun_binary is not None:
            if os.path.isfile(charmrun_binary):
                self._program_args.append(os.path.abspath(charmrun_binary))
            elif _command_line_utility_exists(charmrun_binary):
                self._program_args.append(charmrun_binary)
            else:
                raise ValueError('{} does not refer to an existing file'\
                                 .format(charmrun_binary))

        if (parallel_args is not None) and (len(parallel_args) > 0):
            if charmrun_binary is None:
                raise ValueError("parallel_args can't be specified when no "
                                 "value is specified for charmrun_binary")
            elif isinstance(parallel_args,list):
                self._program_args += parallel_args
            elif isstr(parallel_args):
                self._program_args += parallel_args.strip().split()
            else:
                raise TypeError('parallel_args is expected to be a string or '
                                'a list of strings.')

        if os.path.isfile(enzoe_binary):
            self._program_args.append(os.path.abspath(enzoe_binary))
        elif _command_line_utility_exists(enzoe_binary):
            self._program_args.append(enzoe_binary)
        else:
            raise ValueError('{} does not refer to an existing file'\
                             .format(enzoe_binary))

    def restart_from_ckpt(self, ckpt_name,
                          stdout = None, stderr = None):
        assert os.path.isdir(ckpt_name)
        args = self._program_args + ['+restart', ckpt_name]
        print('Executing: {}'.format(' '.join(args)))
        subprocess.check_call(args, stdout = stdout, stderr = stderr)

    def run(self, param_file, stdout = None, stderr = None):
        assert os.path.isfile(param_file)
        args = self._program_args + [param_file]
        print('Executing: {}'.format(' '.join(args)))
        subprocess.check_call(args, stdout = stdout, stderr = stderr)

def _gen_config_file(template_input_path, config_path, ckpt_cycle, stop_cycle,
                     include_directive_root = None):
    """
    Generates the configuration file

    Notes
    -----
    At the time of writing, there appears to be some bugs in Cello with output 
    directories. At the moment, we can't specify nested output directories.
    For example, we run into problems if we specify something like 
    'test_run/data_%d' or an absolute path to a directory.

    Presently, we can only specify a single level of output directories that 
    are produced within the directory that Cello/Enzo-E is executed in.
    """
    assert stop_cycle > ckpt_cycle
    assert os.path.isfile(template_input_path)

    if include_directive_root is None:
        include_directive_root = ''

    # read in the template configuration file
    with io.open(template_input_path, encoding='utf-8', mode = 'r') as f:
        config = load(f, filename = template_input_path,
                      includedir = include_directive_root)

    # Delete the Stopping and Output sections if they already exist
    for section in ['Stopping', 'Output']:
        if section in config:
            del config[section]

    # Identify the field list
    field_list = config['Field']['list']

    # Now add Stopping and Output Parameters to the configuration
    config['Stopping'] = {'cycle': stop_cycle}
    config['Output'] = {
        'list' : ['checkpoint', 'cycled'],

        'checkpoint' : {
            'type' : 'checkpoint',
            'dir' : ["ckpt-%d", "cycle"],
            'schedule' : {'var' : 'cycle', 'list' : [ckpt_cycle]},
        },

        'cycled' : {
            'type' : 'data',
            'field_list' : field_list,
            'schedule' : {'var' : 'cycle', 'list' : [stop_cycle]},
            'dir' : ['data_%d', 'cycle'],
            'name' : ['data-%d.h5', 'proc']
        }
    }

    # Finally, write the configuration to disk
    with open(config_path,'w') as f:
        dump(config, f)

    expected_checkpoint = 'ckpt-%d' % ckpt_cycle
    expected_output = 'data_%d' % stop_cycle
    return expected_checkpoint, expected_output

_SIM_STDOUT_FILE = 'stdout.log'
_SIM_STDERR_FILE = 'stderr.log'

def _run_restart_sims(enzoe_wrapper, config_fname, output_dirs,
                      expected_ckpt_dir, expected_output_dir):
    """
    Run the checkpoint restart simulations.

    This also makes sure all of the expected files/directories are produced and
    moves/copies files/directories out of the current directory and places them
    in the output directories
    """
    for i,output_dir in enumerate(output_dirs):
        os.mkdir(output_dir)
        # run the actual simulations
        stdout_f = open(_SIM_STDOUT_FILE,'w')
        stderr_f = open(_SIM_STDERR_FILE,'w')
        if i == 0:
            print('Launching first sim')
            enzoe_wrapper.run(config_fname,
                              stdout = stdout_f, stderr = stderr_f)
        else:
            print('Launching restart')
            enzoe_wrapper.restart_from_ckpt(
                expected_ckpt_dir, stdout = stdout_f, stderr = stderr_f
            )
        stdout_f.close()
        stderr_f.close()

        # check that the expected files exist and then move/copy them to the
        # dedicated output directory
        for fname in ['parameters.libconfig', 'parameters.out',
                      _SIM_STDOUT_FILE, _SIM_STDERR_FILE]:
            if os.path.isfile(fname):
                shutil.move(fname, os.path.join(output_dir, fname))

        assert os.path.isdir(expected_ckpt_dir)
        assert os.path.isdir(expected_output_dir)
        shutil.move(expected_output_dir,
                    os.path.join(output_dir, expected_output_dir))

def _compare_files(original_fname, restart_fname):
    assert os.path.abspath(original_fname) != os.path.abspath(restart_fname)
    return _subprocess_exit_success(['diff', original_fname, restart_fname])

def _compare_block_outputs(block_name, original_fname, restart_fname):
    assert os.path.abspath(original_fname) != os.path.abspath(restart_fname)
    args = ['h5diff', '-c', original_fname, restart_fname, '/' + block_name]
    try:
        result = subprocess.check_output(args)
        return len(result) == 0
    except subprocess.CalledProcessError: # exit code is nonzero
        return False

def _read_block_list_file(fname):
    with open(fname,'r') as f:
        return dict([line.split() for line in f])

def identical_data_outputs(pre_restart_path, post_restart_path, output_dir):
    """
    Returns whether the contents of {pre_restart_path}/{output_dir} and
    {post_restart_path}/{output_dir} are identical

    Parameters
    ----------
    pre_restart_path : str
        Path to directory that stores sets of outputs written by the original
        simulation run.
    post_restart_path : str
        Path to directory that stores sets of outputs written by the 
        simulation that was restarted from the checkpoint.
    output_dir : str
        Name of the subdirectory within pre_restart_path and post_restart_path
        that contains the files that are to be compared.

    Returns
    -------
    success: bool
        Indicates whether the data outputs are identical.

    Note
    ----
    This function currently employs the `h5diff` and `diff` command line
    utilities. If either of these can't be found, then a RuntimeError is
    raised.

    We currently don't compare the root-level attributes of each hdf5 output
    files. While this information should be compared, it's not obvious how to
    use h5diff to compare just that information.

    We may want to return object types to represent different classes of errors
    (i.e. to indicate that command line utilities are missing vs an actual
    output error)
    """
    print("Comparing Simulation Outputs")
    if not _command_line_utility_exists('h5diff'):
        raise RuntimeError("This test couldn't be completed because the "
                           "h5diff command line utility can't be found")
    elif not _command_line_utility_exists('diff'):
        raise RuntimeError("This test couldn't be completed because the "
                           "diff command line utility can't be found")

    original_path = os.path.join(pre_restart_path, output_dir)
    restart_path = os.path.join(post_restart_path, output_dir)

    # check for the existence of text files and make sure the libconfig and
    # parameters files are identical
    for suffix in ['libconfig', 'parameters', 'block_list','file_list']:
        basename = '.'.join([output_dir,suffix])
        original_fname = os.path.join(original_path, basename)
        restart_fname = os.path.join(restart_path, basename)
        for fname in [original_fname,restart_fname]:
            if not os.path.isfile(fname):
                print('{} is missing'.format(fname))
                return False

        if suffix in ['libconfig', 'parameters']:
            if not _compare_files(original_fname, restart_fname):
                print('the contents of {} are different'.format(basename))
                return False

    # we skip the comparison between file_list files because that information
    # is redundant with the file_list information (and if a different number of
    # processes are used after restart, they will be different)

    # compare the contents of the block_lists. We can't just use diff because
    # after restarting:
    #   - the ordering of blocks may have been shuffled
    #   - the blocks could have been shuffled between processes and saved to
    #     different processes
    original_fname = os.path.join(original_path, output_dir + '.block_list')
    restart_fname  = os.path.join(restart_path, output_dir + '.block_list')

    original_block_storage = _read_block_list_file(original_fname)
    restart_block_storage = _read_block_list_file(restart_fname)
    if set(original_block_storage) != set(restart_block_storage):
        print("Different sets of blocks were saved to disk before after "
              "restarting")
        return False

    # Finally, let's confirm that the data saved for each block is identical
    for block_name in original_block_storage.keys():
        original_h5_fname = os.path.join(original_path,
                                         original_block_storage[block_name])
        restart_h5_fname = os.path.join(restart_path,
                                        restart_block_storage[block_name])
        print(_compare_block_outputs(block_name, original_h5_fname,
                                      restart_h5_fname))
        for fname in [original_h5_fname, restart_h5_fname]:
            if not os.path.isfile(fname):
                print('{} is missing'.format(fname))
                return False
        if not _compare_block_outputs(block_name, original_h5_fname,
                                      restart_h5_fname):
            print("The data saved for {} changes after the restart"\
                  .format(block_name))
            return False
    return True

_PRERESTART_DIR_TEMPLATE = '{}_PreRestart'
_POSTRESTART_DIR_TEMPLATE = '{}_PostRestart'
_INPUT_FILE_TEMPLATE = '{}_Input.in'

def ckpt_restart_test(test_name, template_input_path, enzoe_wrapper,
                      ckpt_cycle = 2, stop_cycle = 4, clobber = False,
                      include_directive_root = None, delete_data = None):
    """
    Executes the checkpoint restart test

    Parameters
    ----------
    test_name: string
        The name of the test. This determines the output directories.
    template_input_filepath: string
        The path to the enzo-e input file. This file must not include
        scalar/logical expressions or append lists to parameters. If "Stopping"
        "Output" sections are present, they will be overwritten. Not all of
        these conditions can be diagnosed.
    enzoe_wrapper: instance of `EnzoE` class
        Object used to execute simulations
    ckpt_cycle: int, optional
        The cycle where the checkpoint is written to disk. Default is 2.
    stop_cycle: int, optional
        The cycle where the simulation terminates. Default is 4.
    clobber: bool, optional
        Indicates whether to delete existing output directories that are
        presumably left over from a previous run.
    include_directive_root: string, optional
        Root directory used to resolve include directives found in the input
        config file; all relative paths in the directives are treated as though
        they originate from this directory. By default, this is set to the
        current working directory.
    delete_data: bool, optional
        Indicates whether to delete the data afterwards. The default behavior
        is to just delete data when the test fails. When True, the data is
        always deleted. False means that the data should never be deleted.

    Returns
    -------
    success: bool
        Indicates whether the test has passed.
    """

    # make sure that the following files actually exist
    prerestart_dir = _PRERESTART_DIR_TEMPLATE.format(test_name)
    postrestart_dir = _POSTRESTART_DIR_TEMPLATE.format(test_name)
    for output_dir in (prerestart_dir, postrestart_dir):
        if os.path.isdir(output_dir):
            if clobber:
                shutil.rmtree(output_dir)
            else:
                raise RuntimeError(
                    ("The {} directory already exists. Can't overwrite unless "
                     "the clobber argument is specified.").format(output_dir)
                )

    # generate the input file (and get the expected output directories)
    config_fname = _INPUT_FILE_TEMPLATE.format(test_name)
    expected_ckpt_dir, expected_output_dir = _gen_config_file(
        template_input_path, config_fname, ckpt_cycle, stop_cycle,
        include_directive_root = include_directive_root
    )

    # generate temporary directory and change to that directory so that we
    # can execute enzo-e from there
    orig_dir = os.getcwd()
    tmp_dir = tempfile.mkdtemp(prefix = 'ckpt-restart-scratch-',
                               dir = orig_dir)

    os.chdir(tmp_dir)
    try:
        _run_restart_sims(enzoe_wrapper,
                          config_fname = os.path.join('..', config_fname),
                          output_dirs = [os.path.join('..', prerestart_dir),
                                         os.path.join('..', postrestart_dir)],
                          expected_ckpt_dir = expected_ckpt_dir,
                          expected_output_dir = expected_output_dir)
    except:
        os.chdir(orig_dir)
        raise
    os.chdir(orig_dir)

    # delete the temporary directory
    shutil.rmtree(tmp_dir)

    # determine whether the results from the original run and after restarting
    # are identical
    success = identical_data_outputs(prerestart_dir, postrestart_dir,
                                     expected_output_dir)

    # cleanup after ourselves
    if delete_data or (success and delete_data is None):
        for output_dir in (prerestart_dir, postrestart_dir):
            if os.path.isdir(output_dir):
                shutil.rmtree(output_dir)
        os.remove(config_fname)

    return success


def run_ckpt_restart_test(test_name, template_input_path, enzoe_wrapper,
                          test_file, ckpt_cycle = 2, stop_cycle = 4,
                          clobber = False, include_directive_root = None,
                          delete_data = None, command = None):
    """
    Executes the checkpoint restart test and writes the result to an output
    file that is understood by the unit testing framework
    """

    with create_test_report(test_file, clobber = clobber) as test_report:
        if command is not None:
            test_report.write("The command line arguments are as follows:\n")
            test_report.write(command)
            test_report.write('\n', flush = True)
        success = ckpt_restart_test(
            test_name = test_name, template_input_path = template_input_path,
            enzoe_wrapper = enzoe_wrapper, ckpt_cycle = ckpt_cycle,
            stop_cycle = stop_cycle, clobber = clobber,
            include_directive_root = include_directive_root,
            delete_data = delete_data
        )
        if success:
            test_report.passing("The ckpt-restart test has passed.")
            print('Test Passed')
        else:
            test_report.fail("The ckpt-restart test has not passed.")



# Setup the parser for the CLI

_default_enzoe_binary = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),'../bin/enzo-e'
)

_description = """\
Perform a checkpoint-restart test.

The test performs the following procedure:
  1. The specified Cello/Enzo-E input file is copied and modified such that the
     simulation writes a checkpoint output to disk at RESTART_CYCLE, runs until
     STOP_CYCLE and writes a data output including all fields at STOP_CYCLE.
  2. The simulation is run using the modified input file and all outputs are
     moved to a separate directory (a copy of the checkpoint is left in the
     original directory).
  3. The simulation is then restarted and the new outpus are moved to a
     separate directory.
  4. Finally, the program ensures that the output files written at STOP_CYCLE
     contain identical data before and after the restart.

Note: This program executes Cello/Enzo-E in a temporary directory, so that
multiple simultaneous executions of this test should not interfere with each
other.
"""

_TEST_REPORT_TEMPLATE = 'test_{}.unit'

parser = argparse.ArgumentParser(
    description = _description,
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument('test_name', help = 'specifies the name of the test')
parser.add_argument('-i','--input', required = True,
                    help = ('Required option used to specify the template '
                            'input file. This can\'t have any scalar/logical '
                            'expressions or append lists to parameters.'))
parser.add_argument('-o','--output', default = None,
                    help = ('Name of the file where to write the test report. '
                            + 'Default is ' +
                            _TEST_REPORT_TEMPLATE.format('{TEST_NAME}')))
parser.add_argument('-r','--restart_cycle', type = int, default = 2,
                    help = ('The cycle from which the simulation restarts. '
                            'Default is 2.'))
parser.add_argument('-s','--stop_cycle', type = int, default = 4,
                    help = ('The cycle at which data outputs (from before and '
                            'before and after the restart) are compared. '
                            'Default is 4.'))
parser.add_argument('--enzo', default = _default_enzoe_binary,
                    help = ('the path to the enzo-e binary. The default value '
                            'is: "{}"'.format(_default_enzoe_binary)))
parser.add_argument('--charm', default = None,
                    help = ('the path to charmrun. If not specified, the '
                            'simulation will be run serially.'))
parser.add_argument('--parallel-args', nargs = '*', default = [],
                    help = 'arguments to pass to charmrun')
parser.add_argument(
    '--no-clobber', action = 'store_false', dest = 'clobber',
    help = ('Stop the executable if any files/directories persist from '
            'previous executions of this test (including the test report). By '
            'default, they are always removed before the test runs.')
)
parser.add_argument(
    '--include-directive-root', default = None,
    help = ('Root directory used to resolve include directives found in the '
            'input config file; all relative paths in the directives are '
            'treated as though they originate from this directory. By '
            'default, this is the current working directory.')
)

# specify whether to cleanup from the test
cleanup_group = parser.add_mutually_exclusive_group()
cleanup_msg = ()
cleanup_group.add_argument(
    '--cleanup', action = 'store_const', const = True, default = None,
    help = ('By default, the output simulation files are only cleaned up if '
            'the test passes. This option indicates that the outputs should '
            'ALWAYS be cleaned up, while the --no-cleanup option indicates '
            'that the outputs should NEVER be cleaned up.')
)
cleanup_group.add_argument('--no-cleanup', action = 'store_const',
                           const = False, default = None, dest = 'cleanup',
                           help = argparse.SUPPRESS)

if __name__ == '__main__':
    args = parser.parse_args()

    enzoe_wrapper = EnzoE(args.enzo, charmrun_binary = args.charm,
                          parallel_args = args.parallel_args)
    ckpt_cycle = args.restart_cycle
    assert ckpt_cycle > 0
    stop_cycle = args.stop_cycle
    assert stop_cycle > 0

    if args.output is None:
        test_file = _TEST_REPORT_TEMPLATE.format(args.test_name)
        print('No output specified. Writing test report to ' + test_file)
    else:
        test_file = args.output

    run_ckpt_restart_test(
        args.test_name, template_input_path = args.input,
        enzoe_wrapper = enzoe_wrapper, test_file = test_file,
        ckpt_cycle = ckpt_cycle, stop_cycle = stop_cycle,
        clobber = args.clobber,
        include_directive_root = args.include_directive_root,
        delete_data = args.cleanup, command = ' '.join(sys.argv)
    )
