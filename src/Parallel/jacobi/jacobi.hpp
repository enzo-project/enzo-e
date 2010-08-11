/* ------------------------------ PARALLEL ------------------------------ */

/* #define CONFIG_USE_MPI */

#define CONFIG_USE_CHARM


/* ------------------------------ DEFINES ------------------------------ */

#ifdef CONFIG_USE_CHARM
#  define PRINTF CkPrintf
#  define ARGC   main->argc
#  define ARGV   main->argv
#else
#  define PRINTF printf
#  define ARGC   argc
#  define ARGV   argc
#endif


