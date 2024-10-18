import contextlib
from dataclasses import dataclass
import os
import os.path
import shutil
import signal
import subprocess
import tempfile
import time
from typing import Optional

def _effective_parameter_fname(run_dir, parameter_fname, extra_options = []):
    """
    This assists with tweaking the parameters used by an existing file.

    This is inspired by athena++. In the future we should allow enzo-e to
    directly mutate the parameter values by letting us pass extra_options as
    command line arguments. Instead, for the moment we just create a new 
    parameter file with the updated values.

    Parameters
    ----------
    run_dir: str
        Path to the directory where enzo-e will be executed
    parameter_fname: str
        Path to the nominal parameter file
    extra_options: list[str]
        Each element is a key-value pair

    Returns
    -------
    out: str
        Returns the name of the newly created parameter file should be used to
        run the simulation. If a new file is created, it's called 
        f'{run_dir}/parameters.in'
    """
    if len(extra_options) == 0:
        return os.path.abspath(parameter_fname)

    new_parameter_fname = os.path.join(os.path.abspath(run_dir),
                                       'parameters.in')
    with open(new_parameter_fname,'w') as f:
        f.write('# This file was automatically generated\n\n')
        f.write(f'include "{os.path.abspath(parameter_fname)}"\n')
        for elem in extra_options:
            if '=' not in elem:
                raise RuntimeError(f"'{elem}' can't be used as an option. It "
                                   "doesn't have a '=' character")
            param_name, value = elem.split('=', maxsplit=1)

            separator = '='
            if param_name[-1] == '+':
                separator = '+='
                param_name = param_name[:-1]

            name_parts = param_name.split(':')
            assert all(len(part) > 0 for part in name_parts)

            f.write(
                '{ '.join(name_parts) +
                f' {separator} {value};' +
                (' }' * (len(name_parts) - 1))
                + '\n'
            )
    return new_parameter_fname

def _run_simulation(command, max_runtime = float('inf'), sim_name = None,
                    cwd = None, buffer_outputs_on_disk = False):
    """
    Executes a command that runs a simulation

    Parameters
    ----------
    command: str
        The command to execute
    max_runtime: float, Default: float('inf')
    sim_name: str, Default: None
    cwd: str, Default: None
        Forwarded to the `subprocess.Popen` constructor. This overrides the
        current working directory
    buffer_outputs_on_disk: bool, Default: False
        When True, the simulation's stdout and stderr are piped to
        {prefix}/log.out and {prefix}/log.err, respectively, where prefix
        is replaced with cwd (or '.' if cwd is None). If the simulation
        fails, these are forwarded to stdout and stderr, respectively.

    Returns
    -------
    elapsed: float
        The elapsed times in seconds
    """

    if sim_name is None:
        sim_name = 'Simulation'

    with contextlib.ExitStack() as stack:
        if buffer_outputs_on_disk:
            _log_prefix = '%s/log.' % ('.' if cwd is None else cwd)
            f_err, f_out = [
                stack.enter_context(open(_log_prefix + suffix, 'w+'))
                for suffix in ('err', 'out')]
        else:
            f_err, f_out = None, None

        proc = subprocess.Popen(
            command, shell = True, close_fds = True,
            preexec_fn = os.setsid,
            cwd = cwd, stderr = f_err, stdout = f_out)

        stime = time.time()
        while proc.poll() is None:
            if ((time.time() - stime) > max_runtime):
                os.killpg(proc.pid, signal.SIGUSR1)
                raise RuntimeError(
                    f"{sim_name} exceeded max runtime of {max_runtime} "
                    "seconds.")
            time.sleep(1)
        elapsed = time.time() - stime

        if proc.returncode != 0:
            if buffer_outputs_on_disk:
                print("Dumpting stdout:")
                for line in f_out.seek(0):
                    print(line, end = '')
                print("Dumpting stderr:")
                for line in f_err.seek(0):
                    print(line, end = '')

            raise RuntimeError(
                f"{sim_name} exited with nonzero return code "
                f"{proc.returncode}.")
    return elapsed

def _format_command_with_launcher(charmrun_path, enzoe_command, ncpus = 1):
    if charmrun_path is None:
        if ncpus > 1:
            raise ValueError(
                "Can't use more than 1 CPU when charmrun is not provided"
            )
        return enzoe_command
    else:
        if os.path.basename(charmrun_path) == 'charmrun':
            return f"{charmrun_path} ++local +p{ncpus} {enzoe_command}"
        else:
            raise ValueError("Not currently equipped to handle a launcher "
                             "other than charmrun")

@dataclass(frozen = True)
class EnzoEDriver:
    """
    Wraps the Enzo-E binary. This can be used to execute Enzo-E from an 
    arbitrary directory.

    Parameters
    ----------
    enzoe_path: string
        Path to the enzo-e binary.
    charmrun_path: string, optional
        Path to the charmrun binary. If this is not specified, the simulation 
        can't be run in parallel.
    """

    enzoe_path: str
    # charmrun_path should probably be a little more generic. Historically, on
    # clusters when charm++ is built atop MPI, we need to use mpi launcher
    charmrun_path: Optional[str] = None

    def __post_init__(self):

        # check invariants!
        def _check_path(var_name, val):
            if not os.path.isabs(val):
                raise ValueError(f"{var_name} was assigned {val}, which isn't "
                                 "an absolute path")
            elif not os.path.isfile(val):
                raise ValueError(f"the path assigned to {var_name}, {val}, "
                                 "doesn't point to a file")

        _check_path('enzoe_path', self.enzoe_path)
        if self.charmrun_path is not None:
            _check_path('charmrun_path', self.charmrun_path)

    def run(self, parameter_fname, max_runtime = float('inf'), ncpus = 1,
            sim_name = None, cwd = None, extra_options = [],
            buffer_outputs_on_disk = False):
        """
        Run the Enzo-E simulation

        Parameters
        ----------
        parameter_fname: str
        max_runtime: float, Default: float('inf')
        ncpus: int, Default: 1
        sim_name: str, Default: None
        cwd: str, Default: None
            Forwarded to the `subprocess.Popen` constructor. This overrides the
            current working directory
        extra_options: list of strs, Default: []
            List of extra options to use in the simulations (these can
            overwrite values in the parameter file
        buffer_outputs_on_disk: bool, Default: False
            When True, the simulation's stdout and stderr are piped to
            {prefix}/log.out and {prefix}/log.err, respectively, where prefix
            is replaced with cwd (or '.' if cwd is None). If the simulation
            fails, these are forwarded to stdout and stderr, respectively.
        Returns
        -------
        elapsed: float
            The elapsed times in seconds
        """

        if len(extra_options) > 0:
            parameter_fname = _effective_parameter_fname(
                run_dir = cwd if cwd is not None else './',
                parameter_fname = parameter_fname,
                extra_options = extra_options)

        command = _format_command_with_launcher(
            charmrun_path = self.charmrun_path,
            enzoe_command = f"{self.enzoe_path} {parameter_fname}",
            ncpus = 1)

        return _run_simulation(command, max_runtime = max_runtime,
                               sim_name = sim_name, cwd = cwd,
                               buffer_outputs_on_disk = buffer_outputs_on_disk)

    def run_charmrun_restart(self, ckpt_path, ncpus = 1, **kwargs):
        """
        Run the Enzo-E simulation, using a charm++-style restart from the
        specified file

        Parameters
        ----------
        cpkt_path: str
            Path to the checkpoint to restart the simulation from
        ncpus: int, Default: 1
            The number of cpus to use
        **kwargs
            kwargs from _launch_enzoe (other than command)

        Returns
        -------
        elapsed: float
            The elapsed times in seconds
        """
        ckpt_path = os.path.abspath(ckpt_path)

        command = _format_command_with_launcher(
            charmrun_path = self.charmrun_path,
            enzoe_command = f"{self.enzoe_path} +restart {ckpt_path}",
            ncpus = 1)

        return _run_simulation(command, **kwargs)

    def query_uses_double_precision(self):
        # we need to pass ++quiet to prevent charm++'s diagnostic messages from
        # clogging things up
        command = _format_command_with_launcher(
            charmrun_path = self.charmrun_path,
            enzoe_command = f"{self.enzoe_path} ++quiet -precision",
            ncpus = 1)

        rslt = subprocess.run(command, shell = True,
                              capture_output=True).stdout.rstrip()
        if rslt == b'double':
            return True
        elif rslt == b'single':
            return False
        else:
            raise RuntimeError(f"Something went horribly wrong: {rslt}")

    def query_has_grackle(self):
        command = _format_command_with_launcher(
            charmrun_path = self.charmrun_path,
            enzoe_command = f"{self.enzoe_path} -grackle-version",
            ncpus = 1)

        return subprocess.run(command, shell = True).returncode == 0

    def aggregate_params_in_file(self, parameter_fname, cwd = None):
        """
        Execute Enzo-E as a dryrun. In the current directory, this writes 
        parameters.out, an aggregated parameter file, and parameters.libconfig,
        which encodes the same info but in the libconfig format.

        Notes
        -----
        It was initially our intention to write a function that would write the
        libconfig at an arbitrary location, but it turns out that this is 
        somewhat more involved than anticipated if we want to automatically
        avoid overwriting existing files. 
        """

        charmrun_path = None
        if self.charmrun_path is not None:
            charmrun_path = os.path.abspath(self.charmrun_path)
        parameter_fname = os.path.abspath(parameter_fname)
        enzoe_path = os.path.abspath(self.enzoe_path)

        command = _format_command_with_launcher(
            charmrun_path = charmrun_path,
            enzoe_command = f"{enzoe_path} -dryrun {parameter_fname}",
            ncpus = 1)

        rslt = subprocess.run(command, shell = True, cwd = cwd,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.STDOUT)
        if rslt.returncode != 0:
            print("Generation of aggregated parameter file exited with "
                  f"returncode: {rslt.returncode}.\n"
                  f"We tried to execute:{command}.\n"
                  "Combined STDOUT and STDERR follows down below:")
            print(rslt.stdout.decode("utf-8"))
            rslt.check_returncode() # raise a CalledProcessError

def create_symlinks(dst_dir, target_iterable):
    """
    Create a symlink to each path each element of targets_iter within dst_dir

    Parameter
    ---------
    dst_dir: str
        The path to the directory where the symlinks will be made
    target_iterable: iterable of str
        Contains strings that each hold the path to a file or directory for
        which a symlink will be produced.

    Notes
    -----
    This is commonly used when setting up an Enzo-E simulation
    """
    abs_dst_dir = os.path.abspath(dst_dir)
    for src in map(lambda e: os.path.abspath(e), target_iterable):
        dst = os.path.join(abs_dst_dir, os.path.basename(src))
        assert not os.path.exists(dst)
        os.symlink(src, dst)


def grackle_symlink_targets(grackle_input_data_dir):
    """
    generator that yields paths to files in grackle_input_data_dir for which 
    symlinks should be created.

    Notes
    -----
    This is commonly used when setting up an Enzo-E simulation
    """
    with os.scandir(grackle_input_data_dir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            yield entry.path


def create_enzoe_driver_from_args():
    pass
