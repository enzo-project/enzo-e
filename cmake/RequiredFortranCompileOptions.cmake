# this module primarily defines a function that adds the minimal-set of
# compiler-specific flags to a target that are required to compile our fortran
# files.
# -> the compiler will also provide some basic warnings about potential
#    compiler-specific issues


# first, we define helper function(s)

function(fc_warn_untested)
  # this accepts an optional extra argument - a suffix to the warning
  if (${ARGC} EQUAL 0)
    set(extra "")
  elseif (${ARGC} EQUAL 1)
    set(extra " ${ARGV0}") # <- the extra space helps with formatting
  else()
    message(FATAL_ERROR "fc_warn_untested receied unexpected number of args")
  endif()

  # compared to a C compiler, it makes a lot of sense to warn about untested
  # fortran compilers because
  # - we adopt a fairly custom fortran dialect
  # - default name-mangling behavior of Fortran symbols can vary a lot between
  #   different compilers

  message(WARNING
    "${CMAKE_Fortran_COMPILER_ID} Fortran compiler is NOT tested.${extra}")
endfunction()



include(CheckFortranCompilerFlag)

function(_test_and_append_FC_flags flagList resultList)
  # for each item in flagList, this function tests whether the compiler
  # supports the flag. If it is supported, it's append to resultList
  set(l "")
  foreach(flag_name IN LISTS flagList)

    # name of the variable where the result of the check will be stored
    set(flag_support_varname "_FC_supports_${flag_name}")

    # check if the variable is supported
    check_fortran_compiler_flag(${flag_name} ${flag_support_varname})

    if(${flag_support_varname})
      list(APPEND l ${flag_name})
    endif()
  endforeach()

  # append the contents of l to resultList
  set(${resultList} ${${resultList}} ${l} PARENT_SCOPE)
endfunction()


# down below, we implement the main function provided by this module

function(get_required_fortran_options outVar)
  # Overview
  # --------
  # identifies the the minimal set of flags required to succesfully compile
  # Enzo-E's fortran dialect and writes them to outVar
  # -> in other words, these flags only address formatting and name-mangling
  # -> it's considered an anti-pattern to hardcode other kinds of flags (it
  #    creates challenges for the end-user to specify desired flags)
  #
  # We use generator expressions to ensure that these flags are ONLY applied to
  # Fortran source files. This lets the caller blindly apply the results to any
  # target.
  #
  # Note
  # ----
  # I think Enzo-E may have 2 sets of fortran dialects:
  #  - dialect for implementing the ppml code; I think this generally includes
  #    the precision-related flags
  #  - dialect for all other fortran code (e.g. stuff that comes from
  #    enzo-classic); (this is mostly just file-formatting)
  # We also set up flags to control name-mangling of fortran files in both
  # cases (but it may be better to use cmake's builtin FortranCInterface
  # module to handle that)
  #
  # While there don't currently appear to be any incompatibilites, we may need
  # to differentiate between the dialects in the future. Since the files
  # written in these dialects generally compiled in separate, self-contained
  # targets, there are 2 options to address this:
  #   1. make this function accept a keyword argument to control the dialect of
  #      the defined flags and then selectively apply them to the targets
  #   2. define a custom property on the targets (maybe call it
  #      USES_PPML_DIALECT) and then have this function enclose this the
  #      relevant flags within a generator expression that conditionally
  #      enables the flags based on the value of the target

  # report to user what we are doing: we follow convention and repeat the
  # description of what we are doing alongside the result
  set(_fc_flag_msg "Identify required compiler-specific Fortran flags")
  message(STATUS "${_fc_flag_msg}")


  # define the required gnu flags (we will try these out as a fallback if we
  # don't know what else to use)
  if (USE_DOUBLE_PREC)
    set(gnu_precision_flags "-fdefault-real-8;-fdefault-double-8")
  else()
    set(gnu_precision_flags "")
  endif()
  set(all_gnu_flags
    "-fno-second-underscore" "-ffixed-line-length-132" ${gnu_precision_flags})

  # do some (non-exhaustive) handling of the fortran-specific flags!
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(FC_FLAGS ${all_gnu_flags})

  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "INTEL")
    set(FC_FLAGS "-nofor-main") # <-- I think only needed for ppml dialect
    if (USE_DOUBLE_PREC)
      list(APPEND FC_FLAGS "SHELL:-real-size 64 -double-size 64")
    endif()

  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "^PGI|NVHPC$")
    # NVHPC is new branding for PGI compilers
    set(FC_FLAGS "-Mnosecond_underscore;-Mextend")
    if (USE_DOUBLE_PREC)
      fc_warn_untested("Testing GNU flags for double precision.")
      _test_and_append_FC_flags(gnu_precision_flags FC_FLAGS)
    else()
      fc_warn_untested()
    endif()

  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "XL") # IBM compiler
    set(FC_FLAGS "-qfixed=132")
    if (USE_DOUBLE_PREC)
      fc_warn_untested("Testing GNU flags for double precision.")
      _test_and_append_FC_flags(gnu_precision_flags FC_FLAGS)
    else()
      fc_warn_untested()
    endif()

  else()
    message(STATUS
      "No explicit logic for ${CMAKE_Fortran_COMPILER_ID} Fortran compiler, "
      "testing GNU flags")

    set(FC_FLAGS "")
    _test_and_append_FC_flags("${all_gnu_flags}" FC_FLAGS)
  endif()

  # report the results of this function
  if (FC_FLAGS STREQUAL "")
    message(STATUS "${_fc_flag_msg} - none needed")
  else()
    message(STATUS "${_fc_flag_msg} - ${FC_FLAGS}")
    set("${outVar}" "$<$<COMPILE_LANGUAGE:Fortran>:${FC_FLAGS}>" PARENT_SCOPE)
  endif()

endfunction()


function(target_required_fortran_compile_options target flag_scope)
  # Adds required fortran compile options to the target.
  #
  #   target_add_required_fortran_options(<target> <INTERFACE|PUBLIC|PRIVATE>)
  #
  # the required options are identified by the get_requried_fortran_flags
  # function (see that function's docstrings for further details)

  # argument checking:
  if (flag_scope STREQUAL "PUBLIC")
    message(FATAL_ERROR
      "you probably don't want the required fortran flags to be PUBLIC. Try "
      "using INTERFACE or PRIVATE instead")
  elseif (NOT (flag_scope MATCHES "^INTERFACE|PRIVATE$"))
    message(FATAL_ERROR
      "target_required_fortran_compile_options's 2nd argument must be "
      "INTERFACE or PRIVATE")
  endif()

  get_required_fortran_options(_REQ_FC_FLAGS)

  if (NOT ("${_REQ_FC_FLAGS}" STREQUAL ""))
    target_compile_options(${target} ${flag_scope} "${_REQ_FC_FLAGS}")
  endif()
endfunction()