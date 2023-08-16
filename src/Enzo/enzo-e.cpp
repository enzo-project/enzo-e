// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      enzo-e.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @brief     Cello main
///
/// \mainpage Enzo-E / Cello
///
/// <a href="http://cello-project.org/">Cello</a> is an
/// object-oriented adaptive mesh refinement (AMR) software framework
/// for high performance scientific applications.  The framework is
/// scalable, easy to use, and portable across systems ranging from
/// laptops and PC's to the largest HPC systems available, including
/// Blue Waters, the National Science Foundation's Cray petascale
/// supercomputer at the University of Illinois at Urbana-Champaign.
/// Cello's mesh refinement uses the highly-scalable
/// "array-of-octrees" approach, and the Charm++ parallel programming
/// system enables its high parallel scalability.
///
/// Development of Cello is driven by Enzo, a parallel computational
/// astrophysics and cosmology application. The goal is to efficiently
/// map Enzo's multi-resolution multi-physics capabilities onto large
/// parallel computers with potentially millions of computational
/// units. This "petascale" incarnation of Enzo being built on the
/// Cello framework is called Enzo-E.

//----------------------------------------------------------------------

#define CHARM_ENZO

#include "test.hpp"
#include "enzo.hpp"
#include "main.hpp"
#include "charm_enzo.hpp"

#include "../../auto_config.def"

// The following needs to be included once and only once
// This may not be the perfect place for this, but it is when it is included in
// multiple object files
#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

//----------------------------------------------------------------------

extern CProxy_EnzoSimulation proxy_enzo_simulation;
extern CProxy_Simulation     proxy_simulation;

//----------------------------------------------------------------------

static const char* help_message_ = R"HELP(
USAGE:

    charmrun [charm args] %s [-dryrun] <parameter-file>
    charmrun [charm args] %s [-grackle-version | -help | -precision | -version]

DESCRIPTION:
    Runs a simulation.

OPTIONS:
    For the Charm++ arguments (passed to the charmrun launcher), see the
    "Running Charm++ Programs" section of the the Charm++ documentation at:
        https://charm.readthedocs.io/en/latest/charm++/manual.html#running-charm-programs

    <parameter-file>
        The full path to a parameter file

    -dryrun
        Write parameter file to parameters.[out|libconfig] and exit immediately

    -grackle-version
        When this flag is specified, the program prints out the version of
        Grackle that it was linked against and exits. If the program was not
        linked against Grackle, the program exits with a non-zero exit code.

    -help
        Prints this message

    -precision
        Prints whether the program was compiled with single or double precision

    -version
        Prints the program's version information
)HELP";

//----------------------------------------------------------------------

static void print_usage(const char* argv_0){
  const char* prog;
  if (argv_0 == nullptr) {
    prog = "enzo-e";
  } else {
    int start = std::string(argv_0).rfind("/") + 1;
    prog = (start >= 0) ? argv_0 + start : argv_0;
  }

  printf (help_message_, prog, prog);
}

//----------------------------------------------------------------------

struct Args {

  enum class Mode {Error, Help, Run,
                   Dryrun, GrackleVersion, Precision, Version};

  Mode mode;
  const char* param_fname;

  static Args help() { return { Mode::Help,  nullptr }; }
  static Args err()  { return { Mode::Error, nullptr }; }
};

//----------------------------------------------------------------------

/// Parse the command-line arguments
///
/// In the future, it might be cool to do something similar to Athena++ and let
/// pairs that are prefixed by 2 hyphens specify parameter-value pairs that
/// modify the parameter file.
///
/// For that reason, long and short options are both prefixed by 1 hyphen.
static Args parse_args_(int argc, char** argv) {

  auto eq = [=](int i, const char* ref) {return strcmp(argv[i], ref) == 0;};

  // make a quick pass to scan for help command:
  for (int i = 1; i < argc; i++) {if (eq(i, "-help")) return Args::help(); }

  Args::Mode mode = Args::Mode::Run;
  const char* positional_arg_ptr = nullptr;

  for (int i = 1; i < argc; i++) {
    bool is_positional = argv[i][0] != '-';
    if (is_positional && (positional_arg_ptr == nullptr)) {
      positional_arg_ptr = argv[i];
    } else if (is_positional) {
      CkPrintf("ERR: Only 1 positional argument allowed."); return Args::err();
    } else if (argv[i][1] == '\0') {
      CkPrintf("ERR: invalid argument: \"-\"\n");           return Args::err();
    } else if (argv[i][1] == '-') {
      CkPrintf("ERR: invalid option: \"%s\". All options (even long options) "
               "are prefixed by a single hyphen\n", argv[i]);
      return Args::err();
    } else if (mode != Args::Mode::Run) {
      CkPrintf("ERR: too many flags provided\n");           return Args::err();
    } else {
      if      (eq(i, "-dryrun"))           mode = Args::Mode::Dryrun;
      else if (eq(i, "-grackle-version"))  mode = Args::Mode::GrackleVersion;
      else if (eq(i, "-precision"))        mode = Args::Mode::Precision;
      else if (eq(i, "-version"))          mode = Args::Mode::Version;
      else {
        CkPrintf("ERR: unrecognized option: \"%s\"\n", argv[i]);
        return Args::err();
      }
    }
  }

  return {mode, positional_arg_ptr};
}

//----------------------------------------------------------------------
PARALLEL_MAIN_BEGIN
{

  // Initialize parallelization

  PARALLEL_INIT;

#ifdef PNG_1_2_X
  CkPrintf ("PNG_1_2_X\n");
#endif
#ifdef PNG_1_3_X
  CkPrintf ("PNG_1_3_X\n");
#endif
#ifdef PNG_1_4_X
  CkPrintf ("PNG_1_4_X\n");
#endif
#ifdef PNG_1_5_X
  CkPrintf ("PNG_1_5_X\n");
#endif

  Args parsed_args = parse_args_(PARALLEL_ARGC, PARALLEL_ARGV);

  switch (parsed_args.mode) {
    case Args::Mode::Error: print_usage(PARALLEL_ARGV[0]); CkExit(1);

    // handle easy modes that don't care about the positional argument:

    case Args::Mode::GrackleVersion: {
#ifndef CONFIG_USE_GRACKLE
      CkPrintf("Grackle is not linked to this build.\n");
      CkExit(1);
#else
      grackle_version gversion = get_grackle_version();
      CkPrintf ("Grackle Version: %s\n", gversion.version);
      CkPrintf ("Git Branch:   %s\n", gversion.branch);
      CkPrintf ("Git Revision: %s\n", gversion.revision);
      CkExit(0);
#endif
    }

    case Args::Mode::Help:  print_usage(PARALLEL_ARGV[0]); CkExit(0);

    case Args::Mode::Precision:
      if (sizeof(enzo_float) == 4)         { CkPrintf("single\n"); CkExit(0); }
      else if (sizeof(enzo_float) == 8)    { CkPrintf("double\n"); CkExit(0); }
      ERROR("PARALLEL_MAIN_BEGIN", "error in -precision mode");

    case Args::Mode::Version: {
      CkPrintf("%s\n", CELLO_VERSION);
#ifdef CONFIG_HAVE_VERSION_CONTROL
      const char* changeset = CELLO_CHANGESET;
#else
      const char* changeset = "?";
#endif
      CkPrintf ("changeset: %s\n", changeset);
      CkExit(0);
    }

    // now handle modes that accept parameter file as a positional argument:
    case Args::Mode::Dryrun:
    case Args::Mode::Run:
      if (parsed_args.param_fname != nullptr) {
        break;
      } else {
        CkPrintf("ERR: no parameter file was specified\n");
        print_usage(PARALLEL_ARGV[0]);
        CkExit(1);
      }

    default: ERROR("PARALLEL_MAIN_BEGIN", "unhandled case of Args::Mode");
  }

  // the program only reaches this point if parsed_args.mode is either
  //     Args::Mode::Dryrun or Args::Mode::Dryrun

  // Read parameter file

  g_parameters.read(parsed_args.param_fname);
  g_parameters.write("parameters.out",      param_write_cello);
  g_parameters.write("parameters.libconfig",param_write_libconfig);
  g_parameters.write(stdout,param_write_monitor);
  g_enzo_config.read(&g_parameters);

  // Initialize unit testing

  const int ip = CkMyPe();
  const int np = CkNumPes();

  unit_init(ip,np);

  // Initialize Monitor

  monitor_ = Monitor::instance();
  monitor_->set_mode (monitor_mode_root);
  monitor_->header();
  monitor_->print ("","BEGIN ENZO-E");

  // Exit here if -dryrun
  if (parsed_args.mode == Args::Mode::Dryrun) {
    monitor_->print ("ENZO-E","dryrun == true; exiting.");
    p_exit(0);
  }

  // Print initial baseline memory usage

  Memory * memory = Memory::instance();
  monitor_->print("Memory","bytes %ld bytes_high %ld",
		  memory->bytes(), memory->bytes_high());

#ifdef CONFIG_USE_PAPI
  int retval = PAPI_library_init(PAPI_VER_CURRENT);
  if (retval != PAPI_VER_CURRENT && retval > 0) {
    WARNING("Papi::init","PAPI library version mismatch!");
  } else if (retval < 0) {
    WARNING("Papi::init","PAPI initialization error!");
  }
#endif

 //--------------------------------------------------

  proxy_main     = thishandle;

  // --------------------------------------------------
  // ENTRY: #1 Main::Main() -> EnzoSimulation::EnzoSimulation()
  // ENTRY: create
  // --------------------------------------------------
  proxy_simulation = proxy_enzo_simulation = CProxy_EnzoSimulation::ckNew
    (parsed_args.param_fname, strlen(parsed_args.param_fname)+1);
  // --------------------------------------------------

}

PARALLEL_MAIN_END


//======================================================================
#include "enzo.def.h"
//======================================================================
