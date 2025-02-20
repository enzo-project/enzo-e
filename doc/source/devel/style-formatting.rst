.. include:: ../roles.incl

***************
Code Formatting
***************

Enzo-E is configured with multiple automated tools to enforce control formatting (and other properties) of various files in the repository.
At this time, there files fall into 2 broad categories: (i) C/C++ formatting and (ii) miscellaneous formatting.
Tools for formatting more types of files may also be introduced in the future (e.g. formatting python files).

This page is organized as follows:

1. First, we provide a :ref:`high-level overview <high-level-code-formatting>` that describes the simplest way to interact with the formatting tools.
   We describe the simplest possible workflow for invoking the tools.

2. Next, we provide a :ref:`slightly more detailed description <intermediate-level-code-formatting>` of the framework that is used to manage/invoke the tools.
   We also describe the easiest way to invoke the tools on your machine.

3. Finally, we describe some of the :ref:`individual formatting and enforcement tools <individual-enforcement-checks-tools>`.


.. _high-level-code-formatting:

High Level Overview (no local installation)
===========================================

Because we realize that these tools may seem overwhelming or complicated, we have configured the repository to try to simplify the experience to the greatest extent possible.

At a basic level, you don't have to worry about manually invoking any of these tools on your machine.
In fact, you are free to entirely ignore the existence of these tools until it comes time to submit a Pull Request.
When you submit a Pull Request (and whenever you update it), various forms of continuous integration are triggered.

For the present discussion, the `pre-commit.ci <https://pre-commit.ci/>`__ continuous integration tool is of primary interest.
This tool executes all of the formatting tools and if your submission doesn't satisfy all of the requirements, it will fail.
To address formatting issues, you can take 2 approaches:

 1. *(Not-Recommended)*: You can fix the changes locally (and then push your changes).
    You would either read the errors displayed by pre-commit.ci and try to manually fix the problems **OR** you can invoke the linting tools on your local machine (more on that below).

 2. *(Recommended)* You can leave a comment on the Pull Request that simply says ``pre-commit.ci autofix`` and pre-commit.ci will contribute push a commit to your branch that tries to fix the errors.
    While this should address all of the issues in most cases, there are rare cases where you may need manually fix one or two remaining problems.

.. important::

    This second approach modifies the remote-copy of your git-branch that holds your Pull Request's proposed changes (in a similar way to how pressing the "Update branch" button on the webpage for your Pull Request page will modify the branch).

    If you don't know what this means, then you should try to follow the following procedure to avoid any issues (or you can ask for help).
    For concreteness, let's assume that you previously made a Pull Request that proposed changes held in a git-branch called ``my-new-feature``:

     1. First, you need to make sure that the local copy of your ``my-new-feature`` git-branch (the version of the branch on your machine) is synchronized with your remote copy (i.e. the version of the branch on GitHub).
        In other words, you need to ensure that you synchronize any changes that you made since you originally issued the Pull Request.
        If you are using the terminal, invoke ``git checkout my-new-feature`` from anywhere within the Enzo-E repository and then immediately invoke ``git push``.
        If you get a failure message, then you should get help.
        (It's good practice even if you don't think this is necessary)

     2. Navigate to the webpage on GitHub for your pull request (you need to make sure you are logged in).
        If the previous step introduced new changes, then Continuous Integration Checks will be running (if the precommit.ci check is running, you should probably wait until it's done)
        Leave a comment that just says ``pre-commit.ci autofix``.
        Shortly thereafter, you should see a new commit appear on your screen.

     3. Finally, you need to "pull" the new commit from your GitHub repository back to your machine.
        If you are using the terminal, invoke ``git checkout my-new-feature`` from anywhere within the Enzo-E repository and then immediately invoke ``git pull``.

.. _intermediate-level-code-formatting:

At an Intermediate Level (Running Locally)
==========================================

About the precommit software
----------------------------

We rely upon the `pre-commit <https://pre-commit.com/>`__ software to actually manage these enforcement tools and checks.

In more detail, `pre-commit <https://pre-commit.com/>`__ software is a framework that uses a plugin-system (configured through a file called ``.pre-commit-config.yaml``) to manage the invocation and enforcement tools.
It was originally written to help manage when particular software is invoked by one (or more) of `git's hooks <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>__` .
The software's name reflects the fact that it was primarily intended to be used with git's "pre-commit" hook.

The power and utility of the pre-commit software has made it useful outside the context of git-hooks.
In fact, we totally ignore its usage in the context of git-hooks for the remainder of the section (you are free to configure the software that way, but that choice is left up to you).
We primarily use the pre-commit software because of its ability to automatically install locally cached copies of arbitrary versions of automated enforcement tools (this is extremely useful for tools like :ref:`clang-format <clang-format-details>`) and to easily apply all enforcement checks at once.

.. note::

   Because it's confusing, it's worth emphasizing that there are essentially 3 distinct entities named "pre-commit":

   1. the `pre-commit <https://pre-commit.com/>`__ software.
      Enzo-E contributors **only** need to know about this if they want to apply the enforcement tools locally.

   2. the `pre-commit.ci <https://pre-commit.ci/>`__ continuous integration service.
      This is named because the service simply executes the pre-commit software.
      Just about all Enzo-E contributors will encounter this at one time or another.

   3. the "pre-commit" `git hook <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>__`.
      This is one of multiple different "hooks" offered by git.
      The pre-commit software is named after this hook because it was originally designed to be used with this hook.
      We do **NOT** currently recommend using the pre-commit software with the pre-commit hook.


Running the Checks Locally
--------------------------

To run the checks locally, we encourage you need to install the pre-commit software.
This software is written in python and can be installed with ``pip``.
The `installation instructions <https://pre-commit.com/#installation>`__ also mention an alternative approach where you can run download and run pre-commit without fully installing it (as a "zipapp").

Once you have installed ``pre-commit``, you can enforce the checks by invoking the following command from the root of your Enzo-E repository:

.. code-block:: shell-session

   ~/enzo-e $ pre-commit run --all-files

The above command does two things:

 1. it ensures that local copies of the correct versions of the required enforcement tools are installed.
    These local copies are only accessed by pre-commit and won't affect other parts of your system.
    These copies are also cached (so that the tools don't need to be reinstalled on every invocation).

 2. it applies the enforcement tools on the files in your repository (tool-specific exclusions, like files listed ``.clang-format-ignore`` are obviously respected).

.. caution::

   The above command will modify the files in your repository (after all, that's the whole point of the command).
   The pre-commit software does not provide a way to reverse this change.

.. _individual-enforcement-checks-tools:

Enforcement Checks and Enforcement Tools
========================================

Here we briefly discuss some of the individual checks and enforcement tools managed by pre-commit.

.. _clang-format-details:

``clang-format`` (C/C++ Formatting)
-----------------------------------
C/C++ code is formatted by `clang-format <https://releases.llvm.org/18.1.8/tools/clang/docs/ClangFormat.html>`__.

 * At the time of writing, ``clang-format`` enforces formatting rules that are largely derived from the google-style (with a handful of tweaks that derive from the llvm style guide).
   Details about the enforced style are configured in the ``.clang-format`` file at the root of the Enzo-E repository.

 * files that are formatted this way will generally have far fewer merge conflicts.

 * Sometimes, you may need to disable clang-format to disable formatting for individual pieces of code.
   You can do that with a pair of special comments, `described here <https://releases.llvm.org/18.1.8/tools/clang/docs/ClangFormatStyleOptions.html#disabling-formatting-on-a-piece-of-code>`__.

 * Unless you have a very good reason, [#f1]_ **all** newly introduced C/C++ files will be formatted by ``clang-format``.
   With that said, we have temporarily disabled clang-format from applying to most older files (listed in ``.clang-format-ignore``, that existed before we adopted clang-format in order to minimize merge conflicts (we plan to enable ``clang-format`` for these files in the future).
   We plan to enable formatting on these files in the future.

 * **NOTE:** Trying to manually learn all of the style-rules is an exercise in frustration.
   Instead, we recommend that you rely upon clang-format.


.. important::

   The ``clang-format`` version number is tied to the version number of the entire LLVM project, which has a fairly rapid release cadence and  clang-format is **NOT** forward or backwards compatible.
   If you want to manually install and invoke clang-format on your machine (outside of the pre-commit framework), you **NEED** to make sure that you use the exact same version of clang-format that is used by the pre-commit framework.
   If you use a different version, differences **will** arise. [#f2]_
   This version number is stored in the ``.pre-commit-config.yaml`` file by the ``rev`` parameter in the section for the "clang-format plugin".

   For this reason, we strongly encourage you to invoke the ``pre-commit`` tool for local formatting.

miscellaneous checks
--------------------
Some miscellaneous checks are also implemented by a set of miscellaneous enforcement scripts provided by the authors of pre-commit.


.. rubric:: Footnotes

.. [#f1] For example, we might disable ``clang-format`` on source code files that have been generated by external tools like flex/bison.

.. [#f2] While this hasn't been tested, it's *probably* okay to use a different version of clang-format as long as the Major and Minor version numbers match (e.g. if the pre-commit plugin was configured to use version 18.1.2 of clang-format, you could probably use version 18.1.6).


