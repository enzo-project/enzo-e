# See LICENSE_ENZO file for license and copyright information

#.rst:
# GenericOptionCommand
# --------------------
#
# Defines a command, ``generic_option``, that acts like a more general form of
# the built-in ``option`` command.
#
# In slightly more detail, it provides a variable that the user can optionally
# specify::
#
#   generic_option(<variable> <type> <help_text>
#                  DEFAULT_VALUE <default_value>
#                  [FORBID_EXISTING_NONCACHE_VAR])
#
# Similar to the ``option`` command, the default behavior is to
#   - do nothing if ``<variable>`` is already defined (as a normal or cache
#     variable).
#   - initialize ``<variable>`` as a cache variable of type ``TYPE``
#     with the specified ``<default_value>``. If ``DEFAULT_VALUE`` is not
#     specified, then the default is ``OFF`` for a ``BOOL`` variable and
#     ``""`` in other cases.
#
# When the ``FORBID_EXISTING_NONCACHE_VAR`` option is specified, this command
# will immediately abort cmake with an error if ``<variable>`` is already
# exists as a non-cache variable. NOTE: if it exists as a cache variable and as
# a normal variable that is okay (there isn't any way to identify this case
# anyways)

if(__GenericOptionCommand)
  return()
endif()
set(__GenericOptionCommand YES)


function(generic_option variable type help_text)

  # on the off-chance that the variable-name matches a value used in this
  # function, let's check if it already exists (BEFORE DEFINING ANY VARIABLES)
  if (DEFINED "${variable}")
    set(is_already_defined TRUE)
  else()
    set(is_already_defined FALSE)
  endif()

  # parse args
  set(options FORBID_EXISTING_NONCACHE_VAR)
  set(oneValueArgs DEFAULT_VALUE)
  set(multiValueArgs "")
  cmake_parse_arguments(PARSE_ARGV 3 __GENERIC_OPTION "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" )#${ARGN} )

  # some basic error-handling
  set(_funcname "generic_option")
  if (DEFINED __GENERIC_OPTION_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR
      "${_funcname} recieved invalid arguments: "
      "\"${__GENERIC_OPTION_UNPARSED_ARGUMENTS}\"")
  elseif (DEFINED __GENERIC_OPTION_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR
      "${_funcname} received the ${__GENERIC_OPTION_KEYWORDS_MISSING_VALUES} "
      "keyword(s) without any associated arguments.")
  elseif ("${type}" STREQUAL "INTERNAL")
    message(FATAL_ERROR "${_funcname} does not support the INTERNAL type")
  elseif ("${type}" MATCHES "^PATH|FILEPATH$")
    # these types have some weird handling related to relative-paths
    message(FATAL_ERROR
      "${_funcname} doesn't currently support the ${type} type")
  elseif (NOT "${type}" MATCHES "^BOOL|STRING$")
    message(FATAL_ERROR "invalid type passed to ${_funcname}")
  endif()

  if ("${__GENERIC_OPTION_FORBID_EXISTING_NONCACHE_VAR}" AND
      "${is_already_defined}" AND (NOT DEFINED CACHE{${variable}}))
    # by trial and error, the output format now seems nicer (not sure why we
    # need to escape some line-breaks but not others)
    message(FATAL_ERROR "\
The \"${variable}\" variable is expressly forbidden from being defined as a \
regular variable before the appropriate call to ${_funcname}. \
You probably made the mistake of writing some cmake-code with ``set(...)`` \
somewhere to define the variable.
  -> if you are a user, you should assign a value to the variable on the
     the command-line with the -D<var>=<value> flag or by modifying values
     in the CMakeCache.txt file (by modifying it directly OR using some kind
     of gui-program)
  -> if you are a developer, you are probably trying to assign a default value
     to \"${variable}\". In this case, you should prefer to specify the
     initial-value to ${_func_name} by using the DEFAULT_VALUE keyword. If you
     instead choose to refactor, keep in mind that it's a little tricky to get
     this \"right\" (it can be tricky to allow users to override the value)")
  endif()

  cmake_policy(GET CMP0077 _CMP0077_val)
  if ("${_CMP0077_val}" STREQUAL "OLD")
    message(FATAL_ERROR "${_func_name} needs NEW behavior of CMP0077")
  endif()

  if (NOT "${is_already_defined}")
    if (DEFINED __GENERIC_OPTION_DEFAULT_VALUE)
      set("${variable}" "${__GENERIC_OPTION_DEFAULT_VALUE}" CACHE "${type}"
          "${help_text}")
    elseif("${type}" STREQUAL "BOOL")
      set("${variable}" "OFF" CACHE BOOL "${help_text}")
    else()
      set("${variable}" "" CACHE "${type}" "${help_text}")
    endif()
  endif()

endfunction()

