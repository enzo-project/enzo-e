# Enzo-E/Cello Tools

This directory contains numerous command-line tools of varying degree
of uselessness and uptodatedness. Awk scripts in the awk directory are
typically called from higher-level batch scripts, though can be called
directly.

  - `awk/diff-org.awk` : *Convert "git diff" output to org-mode*

      Converts the output of "git diff" to org-mode format as a
      collapsable list of TODO items. Useful for code reviews of
      changes before commiting and pushing to a remote
      repository. Includes sections for unstaged and staged files.

  - `awk/error-org.awk` : *Convert compiler (gcc) errors to org-mode*

      Converts the error messages in output of a source code build to
      org-mode as a collapsable list of TODO items with links to the
      lines in the files indicated by the error message. Called by
      build.sh.

  - `awk/gdb-org.awk` : *Convert gdb frame info to org-mode*

      Converts the stack trace output of "where" in gdb to org-mode as
      a collapsable list of TODO items with links to the lines in the
      files associated with each call in the stack.

  - `awk/grep-org.awk` : *Convert `"grep -n"` output to org-mode*

      Converts the output of "grep -n <text> <file-glob> to org-mode
      format as a collapsable list of TODO items with links to the
      matching lines in the file. Useful, for example, for finding all
      occurrences of identifiers in source code, and keeping track of
      user operations on them, e.g. renaming.

  - `awk/ls-org.awk`

  - `awk/perf-summary.awk`

      Summarizes performance of a simulation Partes the output of a
      simulation, providing a summary of its performance. May be
      out-of-date with current output format, and computations should
      be double-checked for accuracy.

  - `awk/valgrind-org.awk`

      Convert valgrind output to org-mode Converts the output of e.g.
      "valgrind --leak-check=full ..." to emacs org-mode format as a
      collapsable list of todo items.

  - `awk/warning-org.awk` : *Convert compiler (gcc) warnings to org-made*

      Converts the warning messages in output of a source code build
      to org-mode as a collapsable list of TODO items with links to
      the lines in the files indicated by the warning message. Called
      by build.sh.

  - `cello_parse.py`
  - `check-ppm.sh`

      compare PPM files in Enzo-E with ENZO's requires editing
      directory paths, and files have likely been renamed.

  - `ch-perf.py` : *Create performance plots (broken)*

      Parses output from Enzo-E output to create performance plots of
      time for different phases of a simulation. Old and some phases
      are not supported, breaking the generation of plots, though some
      data files are created in perf/ which may be useful.

  - `ch-perf.sh`
  - `ckpt_restart_test.py`
  - `EZPerf`
  - `EZPerf/data-deriv.sh`
  - `EZPerf/data-integ.sh`
  - `EZPerf/ez-perf.sh`
  - `EZPerf/_plot-perf.py`
  - `field_summary.py`
  - `gen_grackle_testing_file.py`
  - `grep-org.sh` : *Create grep.org org-mode file of grep output*

      Single-line call to awk/grep-org.awk

  - `l1_error_norm.py`
  - `plot_mesh8.py`
  - `plot_mesh.py` : *Plots mesh associated with a list of Block names*

      Given a list of block names (e.g. "B1:10_0:01") plot the
      associated mesh hierarchy.

  - `plot_mesh_short.py`
  - `regenerate_parser.py`
  - `run_cpp_test.py`
  - `test_report.py`


