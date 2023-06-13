import os
import os.path
import shutil
import subprocess

import numpy as np
import h5py

from enzoe_driver import create_symlinks

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
        return dict(block_h5_fname_pair(line) for line in f) 

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

def run_ckpt_restart_test(nominal_input, working_dir, enzoe_driver, nproc = 1,
                          ckpt_cycle = 2, stop_cycle = 4, symlink_srcs = [],
                          sim_name_prefix = None,
                          buffer_outputs_on_disk = False):
    """
    Runs the checkpoint-restart test.

    This test consists of 3 steps:
        1. Execute the "ckpt-run": execute enzo-e in `{working_dir}/ckpt_run`
           to generate checkpoint files (and possibly other outputs)
        2. Execute the "restart-run": execute enzo-e in
           `{working_dir}/restart_run` using one of the checkpoint file
        3. Compare the outputs

    Parameters
    ----------
    nominal_input: str
        Path to the nominal configuration file
    working_dir: str
        Path to the directory where we will execute the test
    enzoe_driver: `EnzoE`
        `EnzoE` instance used to drive the test
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
    checkpoint_outputs = [ckpt_cycle, ckpt_cycle+1]

    if sim_name_prefix is None:
        sim_name_ckpt, sim_name_restart = 'ckpt_run', 'restart_run'
    else:
        sim_name_ckpt = f'{sim_name_prefix}_ckpt_run'
        sim_name_restart = f'{sim_name_prefix}_restart_run'

    # TODO: check internal assumptions about the parameter file!
    # - assume that there is no Stopping condition (or if there is one, it only
    #   specifies a cycle)
    # - assume that "order_morton" is the first Method and "check" is the
    #   second method
    # - assume that "order_morton" & "check" don't have schedules or if they
    #   do, that they are configured with list and var = "cycle"


    # Step 1. setup the directories where we will execute the code
    # ============================================================
    ckpt_run_dir, restart_run_dir = [os.path.join(working_dir,f'{e}_run')
                                     for e in ('ckpt', 'restart')]
    os.mkdir(ckpt_run_dir)
    create_symlinks(ckpt_run_dir, symlink_srcs)
    os.mkdir(restart_run_dir)
    create_symlinks(restart_run_dir, symlink_srcs)

    # Step 2. Run the checkpoint run!
    # ===============================
    print(f"Executing checkpoint run, in {ckpt_run_dir}")
    _checkpoint_dir_fmt = "Check-%02d" # the argument is cycle
    _sched_list_str = f'[{checkpoint_outputs[0]}, {checkpoint_outputs[1]}]'
    shared_extra_opts = (
        f'Method:order_morton:schedule:list={_sched_list_str}',
        'Method:order_morton:schedule:var="cycle"',
        f'Method:check:schedule:list={_sched_list_str}',
        'Method:check:schedule:var="cycle"',
        f'Method:check:dir=["{_checkpoint_dir_fmt}", "cycle"]',
        f'Stopping:cycle={stop_cycle}'
    )
    enzoe_driver.run(
        parameter_fname = nominal_input, max_runtime = np.inf,
        ncpus = nproc_ckpt, sim_name = sim_name_ckpt, cwd = ckpt_run_dir,
        extra_options = shared_extra_opts,
        buffer_outputs_on_disk = buffer_outputs_on_disk)

    # Step 3. Run the restart run!
    # ============================
    print(f"Executing checkpoint run, in {restart_run_dir}")
    _restart_from = os.path.abspath(os.path.join(
        ckpt_run_dir, _checkpoint_dir_fmt % ckpt_cycle))
    restart_extra_opts = ('Initial:list=[]', 'Initial:restart=true',
                          f'Initial:restart_dir="{_restart_from}"')
    enzoe_driver.run(
        parameter_fname = nominal_input, max_runtime = np.inf,
        ncpus = nproc_restart, sim_name = sim_name_restart,
        cwd = restart_run_dir,
        extra_options = shared_extra_opts + restart_extra_opts,
        buffer_outputs_on_disk = buffer_outputs_on_disk)


    # Step 4: compare the outputs
    # ===========================
    # do NOT compare the root level parameters.out or parameters.libconf files
    # - we know for a fact that there will be differences in the Initial
    #   section (we could do careful comparisons in the future)
    # Right now, we're just going to compare the restart files
    # - in the future, it would take minimal work to compare old-style output
    #   files as well as new-style output files (all of the necessary machinery
    #   is in place). We just need to figure out what we're comparing)
    print("Comparing outputs:")
    for dir_basename in [_checkpoint_dir_fmt % e for e in checkpoint_outputs]:
        print(f"-> {dir_basename}")
        ckpt_h5obj_map = ckpt_block_file_map(
            f'{ckpt_run_dir}/{dir_basename}/check.file_list')
        restart_h5obj_map = ckpt_block_file_map(
            f'{restart_run_dir}/{dir_basename}/check.file_list')
        assert_equal_distributed_h5(ckpt_h5obj_map, restart_h5obj_map)
