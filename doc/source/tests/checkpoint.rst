.. _checkpoint-tests:

----------------
Checkpoint Tests
----------------

At the time of writing, there are currently two mechanisms for checkpointing:

1. The new approach that writes checkpoints with the ``"check"`` method, that promises a lot more flexibility.
2. The legacy approach that uses Charm++'s machinery

The intention is to entirely transition to the new approach. But, at this time the new approach has limitations (see `Issue #316 <https://github.com/enzo-project/enzo-e/issues/316>`_). For that reason, both approaches are automatically tested (by the pytest machinery).

The checkpoint-restart tests consist of 4 steps:

1. Perform some setup of the temporary directories where the simulation that creates the checkpoints will be run and a temporary directory where the restarted simulation will be run.
   This setup may involve the creation of modified parameter files and the creation of symlinks (As is common with many test-cases, the symlinks allow the include-directives in parameter files to function properly - plus they facillitate the finding of the appropriate grackle data files).

2. Execute the checkpoint-run (this starts from a parameter file and creates checkpoints).

3. Execute the restart-run (this starts a simulation from the one of the checkpoints)

4. Compare the outputs of the two different runs.
   Currently, we check that outputs are bitwise identical (**aside:** this may be a little limitting for functionality involving reductions like gravity methods)

The testing machinery has two major limitations:

1. it generally requires the usage of incomplete parameter files (the test machinery fills in some missing information as it goes)
2. it makes enforces different assumptions about what should be contained inside of a parameter file when testing the scalable-checkpoint approach and when testing the older Charm++-approach.

The first limitation can be somewhat mitigated through the use of our `tools/ckpt_restart_test.py` script.

* This script implements all of the logic from our checkpoint-restart testing machinery, but can be used outside of the pytest machinery.
  It primarily exists to help debug cases where checkpoint-restart breaks.
  It can also be used as a sanity check for whether the checkpoint-restart functionality works properly when using an arbitrary set of methods.

* For example, to run the test for the new-style machinery, you should invoke the following from the base directory of the repository:

  .. code-block:: bash

      python3 tools/ckpt_restart_test.py \
          --input input/Checkpoint/checkpoint_ppm.in \
          --enzoe <path/to/enzo-e> \
          --charm <path/to/charmrun> \
          --stop_cycle 5 \
          --symlink input \
          --grackle-input-data-dir <path/to/grackle/input/data/dir> \
          --legacy-output \
          --test-dir my_test

  The command is very similar to launch the test of checkpoints that use charm++ machinery:

  .. code-block:: bash

      python3 tools/ckpt_restart_test.py \
          --input input/Checkpoint/legacy/checkpoint_ppm.in \
          --enzoe <path/to/enzo-e> \
          --charm <path/to/charmrun> \
          --stop_cycle 5 \
          --symlink input \
          --grackle-input-data-dir <path/to/grackle/input/data/dir> \
          --legacy-output \
          --charm-restart \
          --test-dir my_test

  There are a few things to note about the flags passed to the command:

    * For the particular simulation considered here, the ``--grackle-input-data-dir`` flag is technically unnecessary.
      It's ONLY necessary when you consider a case that involves Grackle.

    * ``--enzoe`` and ``--charm`` should be passed paths to ``enzo-e`` binary and ``charmrun`` binary, respectively

    * the value passed to ``--test-dir`` is somewhat arbitrary.
      It specifies the directory in which the simulations are executed.
      The script will claim full-ownership over this directory.
      If the directory already exists (and has any contents), the tests won't run unless the ``--clobber`` flag is passed (in which case, the script will clear ALL contents of that directory before doing anything else).

* In each case the full parameter file used to launch the initial run (that generates checkpoint-dumps) can be found in ``my_test/ckpt_run/parameters.in`` (the aggregated parameter file produced by that simulation should be located in ``my_test/ckpt_run/parameters.out``).

* For the new-style restarts, the parameter file used to launch the restart can be found in ``my_test/restart_run/parameters.in`` (or ``my_test/restart_run/parameters.out``). For charm-based restarts, there won't be a parameter file.


List of the test cases in Framework
===================================

The following files are used for testing the legacy Checkpoint functionality

* `input/Checkpoint/legacy/checkpoint_boundary.in`
* `input/Checkpoint/legacy/checkpoint_grackle.in`
* `input/Checkpoint/legacy/checkpoint_ppm.in`
* `input/Checkpoint/legacy/checkpoint_vlct.in`

The following files are used for testing the new-style Checkpoint functionality

* `input/Checkpoint/checkpoint_boundary.in`
* `input/Checkpoint/checkpoint_grackle.in`
* `input/Checkpoint/checkpoint_ppm.in`

.. note::

    Introduce a test of the new-style checkpointing for a simulation using the VL+CT Method.

.. note::

    This list should get updated as more tests get introduced. It may also be nice to add descriptions


Tests outside of the framework
==============================

The files `input/Checkpoint/test_cosmo-check.in` and `input/Checkpoint/test_cosmo-restart.in` show a sample-cosmology simulation that uses the checkpoint-restart functionality.
