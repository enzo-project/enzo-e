This directory exists as a centralized location for defining functionallity, implemented in python that may be useful for implementing tests of Enzo-E and implementing tools provided in the root-directory.

To make use of code in this directory from the ``answer_tests`` subdirectory, you can simply import things as you would with any python module.

From other parts of the codebase, you will need to add code like the following before you try importing anything:

    ```python
    import sys
    sys.path.append('path/to/test_utils')
    ```