.. _new-test:

------------------------
How to create a new test
------------------------

Introduction
============

When adding a new function to the enzo-e project there is a possibility that it may affect the rest of the code, for this reason new code must be able to pass all existing tests as well as tests built for the new code. This means that integration tests are required. These tests are designed to stress the system. When making a new test these steps must be followed:

* Create input parameters
* Create a new test subdirectory
* Integrate the test into test infrastructure
* Running your new test


Create Input Parameters for the New Test.
=========================================


Create a new subdirectory for the new tests in the input directory. Create new ``.in`` file that contains input parameters for testing the new feature. This should be a simple problem such as a 2D implosion problem (See enzo-e/input/Hydro/test_implosion.in). More details on how to set input parameters see :ref:`parameter-file` . As an example of a parameter file see load-balance-4.in in the enzo-e/input/Balance directory. This file will be used later on to set the parameters of the test in the running the test section

Create a New Directory
======================

When creating a new test, you must first create a new subdirectory where that test will live. In this subdirectory you will need to place an SConscript file. For an example SConscript file see test/Balance/SConscript. You can copy this file into your new subdirectory. You will then need to make edits in that file as follows

SConscript files start by importing necessary components ::
  Import('env')

  Import('parallel_run')

  Import('serial_run')

  Import('ip_charm')

  Import('bin_path')

  Import('test_path')

This imports the construction environment, the two ways to run the tests and the paths needed for storing the files as exported from SConstruct (the main SCons file)


Integrating the Test into Test Infrastructure
=============================================

In the defines section three variables are defined ::
  env['CPIN'] = 'touch parameters.out; mv parameters.out ${TARGET}.in'

  env['RMIN'] = 'rm -f parameters.out'

  date_cmd = 'echo $TARGET > test/STATUS; echo "----------------"; date +"%Y-%m-%d %H:%M:%S";'

This sets ``CPIN`` up to create the parameters.out files and ``RMIN`` to remove them. ``date_cmd`` writes the test run and time to test/STATUS


The builders and paths are then defined as follows ::
  run_balance_none = Builder(action = "$RMIN; " + date_cmd + parallel_run + "$SOURCE $ARGS > $TARGET 2>&1; $CPIN; $COPY")

  env.Append(BUILDERS = { 'RunBalanceNone' : run_balance_none } )

  env_mv_none = env.Clone(COPY = 'mkdir -p ' test_path + '/Balance/None; mv `ls *.png *.h5`' + test_path +'/Balance/None')


SConscript files are a type of SCons configuration script, typically used in sub directories.
A Builder is an SConscript object called from an SConscript file. It is used to set the relationship between a source and target object. This user defined Builder is set up to remove old parameters.out files, print the test to test/STATUS, run the test in parallel (This is dependent on whether the test is intended to test parallel functionality or serial functionality), takes in a source and a target, then creates a new parameters.out file in <$TARGET>.in
env is the construction environment. By using append the user defined builder is added to the list of builders within the construction environment.
Clone returns a copy of the construction environment with new specifications added to the copy. In this example to make a new directory ``/Balance/None`` and store any ``.png`` or ``.h5`` files created in the test there.

This will need to be edited so that the builder has an appropriate variable name, it is running the correct way (parallel/serial) and the clone needs to have the new directory be made in the new subdirectory and any expected outputs be moved there.

Running Your New Test
=====================

Use the created builder in conjunction with the clone to create the test run variable. This variable will require the test result name(test_WhatsTested.unit), the version of enzo-e to compile and run and the ARGS set to the directory of the previously created ``.in`` file as the input parameters to be run. ::
  balance_none = env_mv_none.RunBalanceNone(

    'test_balance_none.unit',

     bin_path + '/enzo-e',

     ARGS='input/Balance/load-balance-4.in')

This runs the test in the environment specified above (in Clone). It makes the target ``test_balance_none.unit`` which is created by compiling and running enzo-e as bin/enzo-e using ``input/Balance/load-balance-4.in`` as the input file.
     
Clean the test ::
  Clean(balance_none,
        [Glob('#/' + test_path + '/Balance/balance*png'),
	 '#/input/parameters.out',])
	 
This adds the files to the clean list where on calling ``make clean`` the files specified in the Glob(i.e. matching the pattern) will be removed. 

Then in the test directory SConscript link the newly made SConscript ::
  SConscript(Balance/SConscript)

This allows the test to be seen in the main testing framework.  
  
To run all the tests use ``make test`` in the main enzo directory.
