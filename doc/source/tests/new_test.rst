----------------------
How to create new test
----------------------

Introduction
============

When adding a new function to the enzo-e project there is a possiblity that it may effect the rest of the code, for this reason new code must be able to pass all existing tests as well as tests built for the new code. This means that integration tests are required. These tests are designed to stress the system. 

Design input parameters that test new function.
===============================================

Create new ``.in`` file that contains these input parameters in relevent input subdirectory. As an example of a parameter file see load-balance-4.in in the enzo-e/input/Balance directory. This file will be used later on to set the parameters of the test in the running the test section

Defines for new test
====================

In relevent test subdirectory edit SConscript. Create Builder and append to env and env Clone the path for the tests ::
  run_balance_none = Builder(action = "$RMIN; " + date_cmd + parallel_run + "$SOURCE $ARGS > $TARGET 2>&1; $CPIN; $COPY")
  env.Append(BUILDERS = { 'RunBalanceNone' : run_balance_none } )
  env_mv_none = env.Clone(COPY = 'mkdir -p ' test_path + '/Balance/None; mv `ls *.png *.h5`' + test_path +'/Balance/None')


SConscript files are a type of SCons configuration script, typically used in sub directories.
A Builder is an SConscript object called from an SConscript file. It is used to set the relationship between a source and target object. This user defined Builder is set up to remove old parameters.out files, print the test to test/STATUS, run the test in parallel (This is dependent on whether the test is intended to test parallel functionality or serial functionality), takes in a source and a target, then creates a new parameters.out file in <$TARGET>.in
env is the construction enviroment. By using append the user defined builder is added to the list of builders within the construction enviroment.
Clone returns a copy of the construction eniviroment with new specifications added to the copy. In this example to make a new directory ``/Balance/None`` and store any ``.png`` or ``.h5`` files created in the test there.
  
Running test
============

Use the created builder in conjunction with the clone to create test run varible. Ths varible will require the test result name(test_WhatsTested.unit), the path followed and the ARGS set to the directory of the previously created ``.in`` file. ::
  balance_none = env_mv_none.RunBalanceNone(
     'test_balance_none.unit',
     bin_path + '/enzo-p',
     ARGS='input/Balance/load-balance-4.in')

This runs the test in the enviroment specified above (in Clone). It makes the targets ``test_balance_none.unit`` which is stored in the subdirectory and the bin/enzo-p using ``input/Balance/load-balance-4.in`` as the input file.
     
Clean the test ::
  Clean(balance_none,
        [Glob('#/' + test_path + '/Balance/balance*png'),
	 '#/input/parameters.out',])
	 
This adds the files to the clean list where on calling make clean the files specified in the Glob(i.e matching the pattern) will be removed. 

New Directory
=============

If there is not already an appropriate subdirectory create a new subdirectory and in that create a new SConscript, start this file by importing necessary components ::
  Import('env')
  Import('parallel_run')
  Import('serial_run')
  Import('ip_charm')

  Import('bin_path')
  Import('test_path')

This imports the construction enviroment, the two ways to run tests and the paths needed for storing the files.
  
In the defines section three varibles need to be defined ::
  env['CPIN'] = 'touch parameters.out; mv parameters.out ${TARGET}.in'
  env['RMIN'] = 'rm -f parameters.out'

  date_cmd = 'echo $TARGET > test/STATUS; echo "----------------"; date +"%Y-%m-%d %H:%M:%S";'

This sets ``CPIN`` up to create the parameters.out files and ``RMIN`` to remove them. ``date_cmd`` writes the test run and time to test/STATUS

Then follow the steps above and in the test directory SConscript link the newly made SConscript ::
  SConscript(Balance/SConscript)
