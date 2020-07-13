----------------------
How to create new test
----------------------

Design input parameters that test new function.
===============================================

Create new ``.in`` file that contains these input parameters in relevent input subdirectory.

Defines for new test
====================

In relevent test subdirectory edit SConscript. Create Builder and append to env and env Clone the path for the tests::
  run_balance_none = Builder(action = "$RMIN; " + date_cmd + parallel_run + "$SOURCE $ARGS > $TARGET 2>&1; $CPIN; $COPY")
  env.Append(BUILDERS = { 'RunBalanceNone' : run_balance_none } )
  env_mv_none = env.Clone(COPY = 'mkdir -p ' test_path + '/Balance/None; mv `ls *.png *.h5`' + test_path +'/Balance/None')

Running test
============

Use the created builder in conjunction with the clone to create test run varible. Ths varible will require the test result name(test_WhatsTested.unit), the path followed and the ARGS set to the directory of the previously created ``.in`` file.::
  balance_none = env_mv_none.RunBalanceNone(
     'test_balance_none.unit',
     bin_path + '/enzo-p',
     ARGS='input/Balance/load-balance-4.in')

Clean the test ::
  Clean(balance_none,
        [Glob('#/' + test_path + '/Balance/balance*png'),
	 '#/input/parameters.out',])


If there is not already an appropriate subdirectory create a new subdirectory and in that create a new SConscript, then follow the steps above and in the test directory SConscript link the newly made SConscript ::
  SConscript(Balance/SConscript)
