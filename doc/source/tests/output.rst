------------
Output Tests
------------

Tests checking for correct output


output-header
=============



output-stride-1
===============

Runs a simple 2D implosion. Tests output stride  where ``stride_write = 1``. This means testing where specifying that exactly one physical processor is used to output data.

output-stride-2
===============

Similar to above test but on two physical processors where ``stride_write = 2``

output-stride-4
===============

Similar to above test but on four physical processors where ``stride_write = 4``


output-data
===========

Output test of 2D implosion problem
