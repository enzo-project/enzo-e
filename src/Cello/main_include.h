// See LICENSE_CELLO file for license and copyright information

/// @file     main_include.incl
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-08-03
/// @brief    Main entry method include file

  readonly CProxy_Main proxy_main;
  mainchare [migratable] Main {

     entry Main(CkArgMsg *m);

     entry void p_exit (int count_blocks);

     entry void p_checkpoint_output(int count, std::string dir);

     entry void p_initial_exit();
     entry void p_adapt_enter();
     entry void p_adapt_called();
     entry void p_adapt_end();
     entry void p_adapt_update();
     entry void p_adapt_exit();
     entry void p_compute_enter();
     entry void p_compute_continue();
     entry void p_compute_exit();
     entry void p_output_enter ();
     entry void p_output_exit();
     entry void p_stopping_enter();
     entry void p_stopping_balance();
     entry void p_stopping_exit();
     entry void p_text_file_write(int nd, char dir[nd],
     	                          int nf, char file[nf],
     	                          int nl, char line[nl],
				  int count);
     entry void p_exit();

  };
