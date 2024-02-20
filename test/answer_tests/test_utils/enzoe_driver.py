import contextlib
from dataclasses import dataclass
import os
import os.path
import signal
import subprocess
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

        # build the enzo-e part of the command
        enzoe_command = f"{self.enzoe_path} {parameter_fname}"

        if self.charmrun_path is None:
            if ncpus > 1:
                raise ValueError(
                    "Can't use more than 1 CPU when charmrun is not provided"
                )
            command = enzoe_command
        else:
            if os.path.basename(self.charmrun_path) == 'charmrun':
                command = (
                    f"{self.charmrun_path} ++local +p{ncpus} {enzoe_command}")
            else:
                raise ValueError("Not currently equipped to handle a launcher "
                                 "other than charmrun")

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

    def query_uses_double_precision(self):
        # we need to pass ++quiet to prevent charm++'s diagnostic messages from
        # clogging things up
        if self.charmrun_path is None:
            command = f"{self.enzoe_path} ++quiet -precision"
        else:
            command = (f"{self.charmrun_path} ++local +p1 ++quiet " +
                       f"{self.enzoe_path} -precision")

        rslt = subprocess.run(command, shell = True,
                              capture_output=True).stdout.rstrip()
        if rslt == b'double':
            return True
        elif rslt == b'single':
            return False
        else:
            raise RuntimeError(f"Something went horribly wrong: {rslt}")

    def query_has_grackle(self):
        if self.charmrun_path is None:
            command = f"{self.enzoe_path} -grackle-version"
        else:
            command = (f"{self.charmrun_path} ++local +p1 {self.enzoe_path} " +
                       "-grackle-version")
        return subprocess.run(command, shell = True).returncode == 0


def create_symlinks(dst_dir, src_l):
    """
    Create a symlink for each item in src_l

    Notes
    -----
    This is commonly used when setting up an Enzo-E simulation
    """
    abs_dst_dir = os.path.abspath(dst_dir)
    for src in map(lambda e: os.path.abspath(src), src_l):
        dst = os.path.join(abs_dst_dir, os.path.basename(src))
        assert not os.path.exists(dst)
        os.symlink(src, dst)

def create_enzoe_driver_from_args():
    pass


