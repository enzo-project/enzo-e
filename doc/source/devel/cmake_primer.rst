***************
Primer on CMake
***************

There is a lot of excellent documentation for CMake available online.
However, some of that documentation can seem overwhelming.

This primer was written to provide some information at a level between the :ref:`getting_started` tutorial and that documentation.

As an aside, when searching for online CMake documentation and guides, you should generally make sure that the guide was written for version 3.0 or later (sometimes called "Modern CMake").


================================
History of Enzo-E's build system
================================

Historically, Enzo-E made use of the `SCons <https://scons.org/>`_ build system. However, in 2022 Enzo-E migrated to the CMake build system.

While there are a number of viable build systems, there has been a general convergence among C++ developers towards using CMake (this is evidenced by large projects like Boost and Charm++).
Given the very large userbase of CMake, bugs are found and fixed relatively quickly.
Moreover, CMake makes the process of locating and linking against external dependencies relatively easy (especially if they are popular libraries OR were themselves built with CMake).

In contrast, SCons is relatively less popular.
This means that certain areas of the codebase receive less attention (e.g. less widely used compiler toolchains), which can lead to some issues that are more difficult to work around.
Additionally, support staff at clusters are less likely to be familiar with SCons (so they are generally less equipped to help support users)

.. _how_a_cmake_build_works:

============================
How Does a CMake Build Work?
============================

For people mostly familiar with build systems like makefiles, CMake may seem a little unintuitive.
It may be insightful to walk through some commands step-by-step that would be used in a fresh Enzo-E build.

1. Unlike some other build systems, CMake specializes in out-of-source builds (see :ref:`out-of-source_builds` for some additional details/benefits).
   In other words, for a given build CMake only creates/modifies a file in some user-specified directory "build-directory."
   Consequently, the first step in initiating a fresh build is to create this "build-directory" and change the terminal's working directory to your newly created directory.
   Most tutorials suggest that you create a directory called ``build`` in the root directory of the repository. In reality this directory's name is arbitrary and you can pretty much put it wherever you want.

..  code-block:: bash

  # these commands are typically executed from the root directory of the
  # Enzo-E repository
  mkdir build
  cd build

2. The next step is to actually invoke CMake.
   This step involves the configuration of your build (more details about various configuration options for Enzo-E are provided in :ref:`how-to_configure_build`).
   It's important to understand that CMake does not generally compile any code [#f1]_.
   Instead, this step generates build-files (stored within the build-directory) that are natively understood by the user's platform.
   For most Enzo-E user/developers (on unix-like systems), CMake will generate makefiles or (optionally) ninja files.
   When you invoke cmake on the command-line you typically include variables related to the configuration and the final argument specifies the path to the root of the Enzo-E directory.
   If you placed your build file at the root of your Enzo-E directory, this step might look something like:

..  code-block:: bash

  cmake -DCHARM_ROOT=<PATH/TO/charm/build-dir> \
        -DEnzo-E_CONFIG=linux_gcc -DUSE_GRACKLE=OFF ..

3. The last step is to actually invoke the generated build files.
   If you've been following along, this may correlate to a command like:

..  code-block:: bash

  make -j4 # -j4 tells make to execute up to 4 commands in parallel

.. [#f1] Technically, CMake's configuration step may include compilation of some short C-programs.
         This is almost always done in service of testing feature availability of a compiler/platform or testing something about the configuration of a dependency.
         Files that contribute to the main executable/libraries/tests aren't compiled during this stage.


==================================
How are CMake Projects Configured?
==================================

CMake makes use of an interpreted scripting language that is typically found in ``CMakeLists.txt`` text files that are distributed in different directories throughout the project's repository.
The primary ``CMakeLists.txt`` file is found in Enzo-E's root directory.

The scripting language supports many features that include variables, arrays, functions/macros, loops, importing from other files, etc.
Care should be taken when using the scripting features.
Sometimes the scope rules and caching of variables can lead to somewhat unintuitive behavior.

CMake provides built-in commands that are used to provide standardized information.
Examples include ``project`` (define the name of the project), ``add_library`` (define/declare a library that is to be built), ``add_executable`` (define/declare an executable that is to be built), etc.


output targets
~~~~~~~~~~~~~~
The CMake build system revolves around the idea of "targets."

The simplest types of targets are an executable and a static or shared library.
The former can be created with ``add_executable`` and while the latter are created with ``add_library``.
The source files are usually listed when one of these targets is declared.

After targets are declared, then other properties about the targets are defined.
Examples of such properties include:

  * Macro definitions that are passed directly through the compiler (via ``target_compile_definitions``)
  * include directories (via ``target_include_directories``)
  * dependencies on other libraries (via ``target_link_libraries`` -- and yes, this does indeed apply to static libraries)

These properties are typically defined with different "scopes" (``INTERFACE``, ``PUBLIC``, or ``PRIVATE``).
The scope determines if a property affects the target and/or its dependent(s).
When scopes are defined correctly, dependency management becomes straight-forward.

There are also other types of library targets that don't directly correspond to a compiled library.
Examples include ``INTERFACE`` libraries and ``IMPORTED`` libraries.
The former might be used to represent a header-only library or a collection of compiled libraries that you use all at once.
Typically ``IMPORTED`` libraries are used to represent external dependencies.

Note, the usage of Charm++ introduces some additional complexity to this project that is not described here.

.. note:

   At this time, the conventions for organizing Cello/Enzo-E's header files introduce a lot of transitive dependencies.
   Going forward, we may wish to revisit this.
   


``cmake`` directory
~~~~~~~~~~~~~~~~~~~
Following a common convention among CMake projects, we have included a directory called ``cmake`` at the root level of the repository.
This directory holds scripts used for optimizations, building charm++ modules, and locating external dependencies.

``config`` directory
~~~~~~~~~~~~~~~~~~~~
Following a convention from Enzo, the ``config`` directory holds CMake scripts that each define useful variables for the build that are specific to different platforms. 


=========
Questions
=========

.. _out-of-source_builds:

What is an Out-Of-Source Build?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As stated above CMake specializes in out-of-source builds.
In fact, steps have been taken to prevent users from creating in-source builds.

For some context, in-source builds store files created during the build-process (e.g. object files) in the same directories that hold the source files.
It's likely that you encountered this kind of build in a tutorial explaining how to use Make.
For concreteness, Grackle is a software project that employs in-source builds.
Projects that support in-source builds typically need to provide specialized logic for cleaning up from a build (e.g. they usually support ``make clean``).

In contrast, out-of-source builds typically put all files created during the configuration/compilation steps into some build directory.
Because CMake only mutates the contents of the build directory, and it stores all configuration information within the build directory, cleaning up from a build (to start from scratch), is as simple as deleting the build directory.
Out-of-source builds also allow you to maintain multiple different builds at once (that each support incremental recompilation).
For example, you may want to have a separate build directory for each branch that you are actively developing.

The fact that CMake stores a given build's configuration within the build directory also makes it possible to maintain separate build directories dedicated to different configurations.
This might be useful in the following cases:
- while optimizing you might want to have separate build directories dedicated to different compilers (e.g. gcc and icc)
- while implementing a new regression test, you may want to have a build directory supporting single precision and another supporting double precision.
- while working with code related to an optional dependency (e.g. Grackle), it might be helpful to have one build with the dependency and one build without it

.. note:

   As an aside, Enzo-E's build system prior to the migration to CMake also supported out-of-source builds (it did however add some files to the repository's root directory).
   That system had a feature where a new build-directory was automatically created for different git branches (and the system would automatically switch between them).
   That feature is not currently supported.
   Note however, that it was more difficult to maintain builds for different configurations with the old system because of the way that configuration was tracked.


How do I start a new build from scratch?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Just delete your current build directory and make a new one.

As explained in :ref:`out-of-source_builds`, you don't actually need to delete an older build to make a new one from scratch.
You can just make new build directory with a different name.

Creating a fresh build is tedious. How do I make this easier?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It can be a little tedious to start a fresh build if you have defined a number of parameters on the command line.
If it has been a while since your last fresh build it may be difficult to remember all of the configurations.
One way to make this easier is to create a custom machine configuration file.

As is detailed below, :ref:`here <how_do_I_add_a_new_source_file>`, a small oversight previously increased how frequently fresh builds were required.
At the time of writing this page, this should now be less of an issue.


.. _how_do_I_add_a_new_source_file:

What do I need to change when I add a new file to the source code?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The answer to this question is: "usually nothing".
The CMake system for Enzo-E is configured to use globbing expressions to locate source files that are included in the build.
As long your new files are added to one of the main source directories and follows standard naming conventions, the build system should automatically find it.

As an aside, our usage of globbing expressions to locate sources needed by CMake is considered a bad practice by the CMake developers.
Their recommendation is that all files used in the build are explicitly listed.
This is because CMake generally only adjusts the build system if a ``CMakeLists.txt`` file has been updated (e.g. if you append a new entry to the list of source files).
Historically, if you globbed for source files, CMake wouldn't know that it needed to rebuild the list of source files when you added a new one (since the ``CMakeLists.txt`` file wasn't touched).

To mitigate some of these issues we make use of a relatively new flag called the ``CONFIGURE_DEPENDS`` Flag (it's specified with the globbing expressions).
When this flag is present, CMake tells the generator (e.g. makefile, ninja) to rerun the glob expression at build time.
If there are any differences in the recovered set of files, CMake reconstructs the build system.
The CMake developers generally view this unfavorably because:

 * globbing expressions can be expensive (this may be less of an issue on Linux compared to windows)

 * "[t]he ``CONFIGURE_DEPENDS`` flag may not work reliably on all generators"

For that reason, we may want to reconsider our approach in the future.

.. note:
   In the interval of time between Enzo-E's transition to using CMake and the patch introducing this documentation, the command to locate source files with globbing expressions did not include the ``CONFIGURE_DEPENDS`` Flag.
   For reasons described above, steps generally needed to be taken to force the regeneration of the build system (e.g. deleting the build and starting again).
   This should now be less of an issue.


What is Ninja?
~~~~~~~~~~~~~~

`Ninja <https://ninja-build.org/>`_ is an alternative build-system to something like Make.
Ninja was designed to be a smaller, less-feature-rich alternative to Make that is intended to be used with tools (like CMake) that generate its input files.
These design goals can facilitate faster build-times in large projects.

To use ``Ninja``, you need to make sure it is installed on your machine and you need to specify ``-G Ninja`` as one of the configuration arguments in the second step described in :ref:`how_a_cmake_build_works`.
Then in the build step, replace ``make`` with ``ninja``.

Does CMake save a log of compiler outputs to any files?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The short answer is no.
Unlike the old build system, the current system does not save the compilation outputs to disk.
If this is something that would be broadly useful, we might be able to hack something together that accomplishes this.
