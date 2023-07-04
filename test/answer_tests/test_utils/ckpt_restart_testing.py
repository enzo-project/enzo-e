import os
import os.path
import shutil
import subprocess
import tempfile

import libconf
import numpy as np
import h5py

from .enzoe_driver import create_symlinks

def ckpt_block_file_map(filelist_path):
    """
    Builds the h5object-map for hdf5 files that serve as an Enzo-E checkpoint
    """

    out = {}
    h5_fnames = []

    # open the filelist file - it should be called check.file_list
    with open(filelist_path, 'r') as f_filelist:
        dirname = os.path.dirname(os.path.abspath(filelist_path))

        # the first line should specify the number of h5 files
        n_files = int(f_filelist.readline())

        # the next lines specify the prefix shared by hdf5/block_list files
        for fname_prefix in map(lambda line: line.rstrip(), f_filelist):
            full_prefix = os.path.join(dirname, fname_prefix)

            h5_fnames.append(full_prefix + '.h5')
            assert os.path.isfile(h5_fnames[-1])

            # determine which files are held by current h5 file by reading
            # the blocklist file (alternatively, we could just check h5 file)
            with open(full_prefix + '.block_list', 'r') as f_blocklist:
                for line in f_blocklist:
                    block, block_level = line.rstrip().split(' ')
                    assert block not in out # SANITY CHECK!
                    out[block] = h5_fnames[-1]
    assert n_files == len(h5_fnames) # sanity check!
    return out

def legacy_output_block_file_map(blocklist_path):
    """
    Builds the h5object-map for hdf5 files were written through the legacy
    output approach.

    The legacy output system created a unified .block_list that contained info
    about the locations of all blocks (the .file_list file is redundant!)
    """
    with open(blocklist_path, 'r') as f:
        dirname = os.path.dirname(os.path.abspath(blocklist_path))

        def block_h5fname_pair(line): # line fmt: "<block_name> <h5basename>\n"
            block, h5_basename = line.rstrip().split(' ')
            return (block, os.path.join(dirname, h5_basename))
        return dict(block_h5fname_pair(line) for line in f)

def _compare_block_outputs(block_name, f1, f2):
    # we choose to use h5diff instead of writing custom python code, because
    # I'm more confident in the accuracy of h5diff

    assert os.path.abspath(f1) != os.path.abspath(f2)
    args = ['h5diff', '-c', f1, f2, '/' + block_name]
    try:
        result = subprocess.check_output(args)
        return len(result) == 0
    except subprocess.CalledProcessError: # exit code is nonzero
        return False

def _dict_differences(ref_dict, other_dict):

    def _is_equal(a, b):
        a_isarr, b_isarr = isinstance(a, np.ndarray), isinstance(b, np.ndarray)
        if (a_isarr and b_isarr and np.array_equal(a, b, equal_nan = True)):
            return True
        return (not (a_isarr or b_isarr)) and (a == b)

    missing_keys, unequal_keys = [], []
    for k,v in ref_dict.items():
        if k not in other_dict:
            missing_keys.append(k)
        elif not _is_equal(v, other_dict[k]):
            unequal_keys.append(k)

    if (len(missing_keys) > 0) or (len(ref_dict) != len(other_dict)):
        return missing_keys + [k for k in other_dict if k not in ref_dict]
    return unequal_keys

def _report_h5_attr_diff(f, f2, dset, attr_l):
    def _to_str(v):
        if isinstance(v, np.ndarray):
            return np.array2string(v, floatmode = 'unique')
        return repr(v)

    fname1, fname2 = f.filename, f2.filename

    if ((attr_l[0] in f[dset].attrs) and (attr_l[0] in f2[dset].attrs)):

        msg = [f"attributes of <{dset}> have discrepant values in " +
               f"{fname1} & {fname2}"]                        
        for a in attr_l:
            msg.append(f'   "{a}": {_to_str(f[dset].attrs[a])}  ' +
                       f'{_to_str(f2[dset].attrs[a])}')
        raise AssertionError('\n'.join(msg))
    else:
        raise AssertionError('The following attribute(s) appear in the '
                             f'<{dset}> of only {fname1} or {fname2} (they '
                             f'should be in both): [{", ".join(attr_l)}]')

def _unique(itr):
    tmp = set()
    for e in itr:
        if e not in tmp:
            tmp.add(e)
            yield e

def assert_equal_distributed_h5(actual_h5obj_map, desired_h5obj_map):
    """
    Checks that the hdf5 objects (and attributes) distributed across 2 distinct
    sets of hdf5 files are all identical.

    This requires that the root-level attributes are identical across ALL
    objects and that there are NO root level datasets.

    Internally, this makes calls to the h5diff command-line utility.
    """

    actual_fnames = list(_unique(actual_h5obj_map.values()))

    # check consistency of attributes in at the root level in all hdf5 files
    with h5py.File(actual_fnames[0], 'r') as f:

        # check consistency of attributes at the root level against other hdf5
        # files in the checkpoint group (this is a sanity check!)
        for fname in actual_fnames[1:]:
            with h5py.File(fname, 'r') as f2:
                problem_attrs = _dict_differences(f.attrs, f2.attrs)
                if len(problem_attrs) > 0:
                    _report_h5_attr_diff(f, f2, '/', problem_attrs)

        # check consistency of attributes at the root level against hdf5 files
        # written after the restart
        for fname in _unique(desired_h5obj_map.values()):
            with h5py.File(fname, 'r') as f2:
                problem_attrs = _dict_differences(f.attrs, f2.attrs)
                if len(problem_attrs) > 0:
                    _report_h5_attr_diff(f, f2, '/', problem_attrs)

    # In an ideal world, we would also confirm that:
    # - there are no objects at the root level other than the ones specified in
    #   the object maps

    # now, let's check consistency of everything else (we need to do this
    # comparison on an h5object-by-h5object basis since h5objects can be
    # grouped together differently in the "actual" and "desired" cases
    assert len(actual_h5obj_map) == len(desired_h5obj_map)

    for h5obj_name, actual_fname in actual_h5obj_map.items():
        desired_fname = desired_h5obj_map[h5obj_name]
        if not _compare_block_outputs(h5obj_name, actual_fname, desired_fname):
            raise AssertionError(
                f"There are differences in the <{h5obj_name}> data saved in "
                f"the {actual_fname} and {desired_fname}"
            )

def enforce_assumptions_and_query(enzoe_driver, parameter_fname, scratch_dir,
                                  disallowed_output_subsections = [],
                                  use_charm_restart = False):
    """
    This function parses the parameter files and checks some assumptions and
    queries some basic information.

    Note: to actually query this information, we need to actually a dryrun of
    Enzo-E in a directory where all of the symlinks are properly set up. This
    is intended to be done in a scratch directory (to avoid overwriting
    something important in the future if the default behavior of Enzo-E 
    accidentally changes).

    While tools/cello_parse.py could be used for this, that tools doesn't
    handle a few edge cases (this approach should be a little more robust).

    Returns
    -------
    method_list:
        list of methods in the parameter file. This is returned with the
        intention of letting us modify the order of Methods, if necessary.
    field_list:
        list of fields
    """

    # convert to libconf-format, parse the libconf file, and delete temp files
    enzoe_driver.aggregate_params_in_file(parameter_fname, cwd = scratch_dir)
    with open(os.path.join(scratch_dir, 'parameters.libconfig')) as f:
        config = libconf.load(f)

    # Check our 3 assumptions
    if any(k != 'cycle' for k in config.get('Stopping', {}).keys()):
        raise ValueError(
            "There must not be any Stopping criterion, OR if it exists, it "
            "by cycled-based (so that it can be overwritten)."
        )

    method_list = config.get('Method', {}).get('list', [])
    if use_charm_restart and ("check" in method_list):
        raise ValueError(
            '"check" can\'t be in Method:list for charm-based restarts')
    elif (not use_charm_restart) and (("order_morton" in method_list) or
                                      ("check" in method_list)):
        if ( (len(method_list) < 3) or (method_list[0] != "order_morton") or
             (method_list[1] != "check") ):
            raise ValueError(
                'If either "order_morton" and "check" is present in '
                'Method:list, we require them to be The first and second '
                'elements respectively'
            )

        for name in ["order_morton", "check"]:
            schedule = config["Method"][name].get("schedule", {})
            if (any(k not in ['var', 'list'] for k in schedule.keys()) or
                (schedule.get('var', 'cycle') != 'cycle')):
                raise ValueError(
                    f'The "{name}" method can only be configured with a '
                    'schedule that uses a list of cycles. Alternatively, it '
                    'shouldn\'t be configured with any schedule')
    elif (not use_charm_restart):
        # we can hopefully relax this requirement in the future!
        raise ValueError(
            'currently "order_morton" and "check" must both be in Method:list '
            'to run new-style checkpoint-restart tests')

    output_section = config.get('Output', {})
    if output_section.get('list',[]):
        raise ValueError("Ouput:list must be unset or be empty")
    for subsection in disallowed_output_subsections:
        if len(output_section.get(subsection, [])) > 0:
            raise ValueError("No parameters are allowed to be set within "
                             f"the Output:{subsection} parameter group")

    field_list = config.get('Field', {}).get('list', [])
    if len(field_list) == 0 and use_charm_restart:
        raise ValueError("Field:list should contain at least 1 element when "
                         "testing charm-based checkpoints")
    return method_list, field_list

def _list_to_paramstr(l):
    """Convert a list to a string that can be parsed in parameter file"""
    def _tostr(arg):
        if isinstance(arg,str):
            return f'"{arg}"'
        else:
            return f'{arg}'
    return "[" + ", ".join(_tostr(e) for e in l) + "]"

def fetch_extra_ckpt_run_opts(checkpoint_dir_fmt, checkpoint_outputs,
                              stop_cycle = 4, field_list = [],
                              legacy_output_dir_fmt = None,
                              use_charm_restart = False):

    shared_extra = ()
    if (legacy_output_dir_fmt is not None) and (len(field_list) > 0):
        shared_extra = (
            'Output:h5dump:type="data"',
            f'Output:h5dump:field_list={_list_to_paramstr(field_list)}',
            f'Output:h5dump:schedule:list=[{stop_cycle - 1}]',
            'Output:h5dump:schedule:var="cycle"',
            f'Output:h5dump:dir=["{legacy_output_dir_fmt}", "cycle"]',
            'Output:h5dump:name=["data-%d.h5", "cycle"]',
            'Output:list=["h5dump"]'
        )

    if use_charm_restart:
        assert (len(shared_extra) > 0) and (len(checkpoint_outputs) == 1)
        return shared_extra + (
            'Output:list=["checkpoint","h5dump"]',
            'Output:checkpoint:type="checkpoint"',
            f'Output:checkpoint:schedule:list=[{checkpoint_outputs[0]}]',
            'Output:checkpoint:schedule:var="cycle"',
            f'Output:checkpoint:dir=["{checkpoint_dir_fmt}", "cycle"]',
            f'Stopping:cycle={stop_cycle}',
        )
    else:
        _sched_list_str = _list_to_paramstr(checkpoint_outputs)
        return shared_extra + (
            f'Method:order_morton:schedule:list={_sched_list_str}',
            'Method:order_morton:schedule:var="cycle"',
            f'Method:check:schedule:list={_sched_list_str}',
            'Method:check:schedule:var="cycle"',
            f'Method:check:dir=["{checkpoint_dir_fmt}", "cycle"]',
            f'Stopping:cycle={stop_cycle}'
        )

def run_ckpt_restart_test(nominal_input, working_dir, enzoe_driver,
                          use_charm_restart = False, nproc = 1,
                          ckpt_cycle = 2, stop_cycle = 4, symlink_srcs = [],
                          sim_name_prefix = None,
                          legacy_output_dir_fmt = None,
                          buffer_outputs_on_disk = False):
    """
    Runs the checkpoint-restart test.

    This test consists of 3 steps:
        1. Execute the "ckpt-run": execute enzo-e in `{working_dir}/ckpt_run`
           to generate checkpoint files (and possibly other outputs)
        2. Execute the "restart-run": execute enzo-e in
           `{working_dir}/restart_run` using one of the checkpoint file
        3. Compare the outputs
           - When legacy_output_dir_fmt is not `None`, data dumps written with
             the legacy output machinery (during both the the checkpoint and
             restart runs) will also be compared
           - When use_charm_restart is `False`, the contents of the checkpoint
             files are also compared

    Parameters
    ----------
    nominal_input: str
        Path to the nominal configuration file
    working_dir: str
        Path to the directory where we will execute the test
    enzoe_driver: `EnzoE`
        `EnzoE` instance used to drive the test
    use_charm_restart: bool, default: False
        When true, use the checkpoint-restart functionality built into charm++
    nproc: int or tuple of ints
        The number of processors to use to drive the simulations. When this is
        a tuple, the first (second) element specifies how many processors to
        use while executing the ckpt-run (restart-run)
    ckpt_cycle, stop_cycle: int
        Specifies the cycle to perform the restart at and the cycle to end the
        simulation at.
    symlink_srcs: list of str
        Specifies the paths to files/directories that symlinks need to be
        created for in the ckpt_run and restart_run directories
    sim_name_prefix: str, default: None
        Optional simulation name. This is only used when logging and reporting
        errors.
    legacy_output_dir_fmt: str, default None
        When specified, this should be a string used to specify the directory
        name where outputs are written using the legacy-Output machinery. When
        specified the string be printf-style formatting string that expects a
        single integer (they cycle number)
    buffer_outputs_on_disk: bool, default: False
        When True, the simulations' stdout and stderr are piped to
        {working_dir}/{run_dir}/log.out and {working_dir}/{run_dir}/log.err,
        respectively, where {run_dir} is replaced by ckpt_run and restart_run,
        based on which replaced simulations are run. If the tests fail, these
        are forwarded to stdout and stderr. When False, the simulations' stdout
        stderr are directly forwarded to stdout and stderr.
    """
    # Preliminary argument checking!
    try:
        _nproc = int(nproc)
        assert _nproc == nproc
        nproc_ckpt, nproc_restart = _nproc, _nproc
    except TypeError:
        nproc_ckpt, nproc_restart = nproc
    assert (nproc_ckpt > 0) and (nproc_restart > 0)

    assert 0 < ckpt_cycle
    assert ckpt_cycle+2 <= stop_cycle

    # with the new-style restarts, we write 2 checkpoint files since we
    # actually compare the contents of the files
    if use_charm_restart:
        checkpoint_outputs = [ckpt_cycle]
    else:
        checkpoint_outputs = [ckpt_cycle, ckpt_cycle+1]

    if sim_name_prefix is None:
        sim_name_ckpt, sim_name_restart = 'ckpt_run', 'restart_run'
    else:
        sim_name_ckpt = f'{sim_name_prefix}_ckpt_run'
        sim_name_restart = f'{sim_name_prefix}_restart_run'

    # Step 1. setup the directories where we will execute Enzo-E
    # ==========================================================
    def _prep_dir(basedirname):
        dir_path = os.path.join(working_dir, basedirname)
        os.mkdir(dir_path)
        create_symlinks(dir_path, symlink_srcs)
        return dir_path
    ckpt_run_dir = _prep_dir('ckpt_run')
    restart_run_dir = _prep_dir('restart_run')

    # ASIDE: query info about nominal configuration and enforce assumptions
    query_param_dir = _prep_dir('query_dir')
    disallowed_output_subsections = []
    if legacy_output_dir_fmt is not None:
        disallowed_output_subsections.append('h5dump')
    if use_charm_restart:
        disallowed_output_subsections.append('checkpoint')
    in_method_list, in_field_list = enforce_assumptions_and_query(
        enzoe_driver, nominal_input,
        scratch_dir = query_param_dir,
        disallowed_output_subsections = disallowed_output_subsections,
        use_charm_restart = use_charm_restart,
    )

    # Step 2. Execute the checkpoint run!
    # ===================================
    print(f"Executing checkpoint run, in {ckpt_run_dir}")
    _checkpoint_dir_fmt = "Check-%02d" # the argument is cycle
    ckpt_run_extra_opts = fetch_extra_ckpt_run_opts(
        checkpoint_dir_fmt = _checkpoint_dir_fmt,
        checkpoint_outputs = checkpoint_outputs,
        stop_cycle = stop_cycle, field_list = in_field_list,
        legacy_output_dir_fmt = legacy_output_dir_fmt,
        use_charm_restart = use_charm_restart)
    enzoe_driver.run(
        parameter_fname = nominal_input, max_runtime = np.inf,
        ncpus = nproc_ckpt, sim_name = sim_name_ckpt, cwd = ckpt_run_dir,
        extra_options = ckpt_run_extra_opts,
        buffer_outputs_on_disk = buffer_outputs_on_disk)

    # Step 3. Execute the restart run!
    # ================================
    print(f"Executing checkpoint run, in {restart_run_dir}")
    _restart_from = os.path.abspath(os.path.join(
        ckpt_run_dir, _checkpoint_dir_fmt % ckpt_cycle))
    if use_charm_restart:
        enzoe_driver.run_charmrun_restart(
            ckpt_path = _restart_from, ncpus = nproc_restart,
            max_runtime = np.inf, sim_name = sim_name_restart,
            cwd = restart_run_dir,
            buffer_outputs_on_disk = buffer_outputs_on_disk)
    else:
        restart_extra_opts = ('Initial:list=[]', 'Initial:restart=true',
                              f'Initial:restart_dir="{_restart_from}"')
        enzoe_driver.run(
            parameter_fname = nominal_input, max_runtime = np.inf,
            ncpus = nproc_restart, sim_name = sim_name_restart,
            cwd = restart_run_dir,
            extra_options = ckpt_run_extra_opts + restart_extra_opts,
            buffer_outputs_on_disk = buffer_outputs_on_disk)

    # Step 4: compare the outputs
    # ===========================
    if use_charm_restart:
        # there is nothing to compare since we restart from a different
        # directory
        #
        # NOTE: Actually, https://github.com/enzo-project/enzo-e/issues/5
        # and https://github.com/enzo-project/enzo-e/pull/45 suggest that
        # Enzo-E may actually write new parameters.{out|libconf} files
        # following a charm-based restart. But, I guess it may overwrite the
        # original version of the file
        #
        # We may want to address this in the future!
        pass
    else:
        # do NOT compare the root-level parameters.out or parameters.libconf
        # files -> we know for a fact that there will be differences in the
        # Initial section (we could do careful comparisons in the future)

        # Compare the restart files
        print("Comparing outputs from the \"check\" method:")
        for dir_basename in [_checkpoint_dir_fmt % e
                             for e in checkpoint_outputs]:
            print(f"-> {dir_basename}")
            ckpt_h5obj_map = ckpt_block_file_map(
                f'{ckpt_run_dir}/{dir_basename}/check.file_list')
            restart_h5obj_map = ckpt_block_file_map(
                f'{restart_run_dir}/{dir_basename}/check.file_list')
            assert_equal_distributed_h5(ckpt_h5obj_map, restart_h5obj_map)


    if legacy_output_dir_fmt is not None:
        print("Comparing files written by legacy output machinery:")
        for dir_basename in [legacy_output_dir_fmt % (stop_cycle - 1)]:
            print(f"-> {dir_basename}")
            ckpt_h5obj_map = legacy_output_block_file_map(
                f'{ckpt_run_dir}/{dir_basename}/{dir_basename}.block_list')
            restart_h5obj_map = legacy_output_block_file_map(
                f'{restart_run_dir}/{dir_basename}/{dir_basename}.block_list')
            assert_equal_distributed_h5(ckpt_h5obj_map, restart_h5obj_map)
