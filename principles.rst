Enzo-P / Cello Guiding Principles
=================================

Key guiding principles in Enzo-P / Cello design and construction are
performance, scalability, usability, reliability, and flexibility.

Performance
-----------

Cello enables high single-node performance through the design of its
field classes for representing arrays of field variables (density,
velocity, etc.) that are associated with blocks of the AMR hierarchy.
The size of these block-local arrays is specified by the user in
Enzo's configuration file, allowing flexibility for optimizing for
cache size and parallel task granularity.  Cello's field data also
supports improving performance through data alignment and data
structure padding.  *Data alignment*---aligning array elements along
128-bit boundaries---can greatly improve performance of a CPU's SIMD
instructions.  *Data structure padding*---separating field arrays by
some small multiple of the cache line length---can reduce cache
conflict and improve cache reuse.  Block size, data alignment, and
data structure padding are all controllable through Enzo-P's parameter
file.

Scalability
-----------

Cello enables high scalability in part through the innovative
combination of Charm++ and Cello's array-of-octree AMR.  Charm++, the
externally-developed parallel programming system on which Cello is
built, is a data-driven, asynchronous, programming language based on
C++.  Charm++ incorporates state-of-the-art dynamic load balancing,
its asynchronous execution model is naturally latency-tolerant, and it
provides automatic overlap of computation and communication.  Emergent
scalability issues such as fault-tolerance and energy efficiency are
also being researched and incorporated into Charm++, allowing Enzo-P /
Cello to benefit automatically as the parallel programming system
continues to improve.

The *array-of-octree* AMR approach that is implemented in Cello has
been shown to be among the most scalable known AMR approaches to date.
The specific AMR remeshing algorithm is based on that developed and
prototyped by the Charm++ group (see `Langer, et al
<http://charm.cs.illinois.edu/papers/12-35>`_).  As in their approach,
our mesh refinement algorithm proceeds mostly asynchronously, using
only a lightweight barrier-like operation provided by the Charm++
runtime system to detect consensus on refinement levels between
blocks.

Usability
---------

.. |NSF000| image:: nsf-000.png
   :scale: 50

.. |NSF300| image:: nsf-300.png
   :scale: 50

.. table:: A simple hydrodynamics problem

   ===========  ===========
     |NSF000|      |NSF300| 
   ===========  ===========

Cello is designed to be both *user-friendly* and *developer-friendly*.
The integrated Enzo-P / Cello structured configuration files are easy
to read and write, and allow for much more powerful yet
easy-to-understand problem initialization than in Enzo.  As an
example, the simple test problem whose output is shown above used the
combination of a PNG image file created using a paint program,
together with the following configuration to initialize the density
field:

.. code:: python

   Initial {
      density { value = [ 1.0 + x, "input/density_mask.png", 0.125 ]; }
   }

.. role:: config(code)
   :language: python
 
This reads naturally as "initialize density with the value 1.0
wherever the input image mask is set (i.e. not black), and 0.125
elsewhere."

Logical expressions such as :config:`x + sin(y) < 0.5` can be used in
place of image masks, and multiple masks and logical expressions can
be combined to generate arbitrarily complex initial conditions.  This
capability greatly simplifies the user experience of setting up a test
problem, making Enzo-P not only useful for researchers to run
scientifically viable simulations, but also for students and educators
to learn and experiment with concepts in computational fluid dynamics
and astrophysics.

Migrating new functionality from Enzo or other science applications
to Enzo-P / Cello is also very developer-friendly.  Only the serial
code implementing a numerical method on a single grid patch in Enzo
is required, and frequently no modifications of the method code are
even necessary.  The bulk of the migration effort required is writing
an EnzoBlock class that contains static variables to replace
Enzo's global variables, but this has already been done.  Even Enzo's
troublesome redefinition of `float` to `double` can be
easily handled, by globally replacing `float` with the macro
`enzo_float` in the code being migrated.

Reliability
-----------

Cello currently contains a large and growing collection of test
problems, together with built-in support for both running tests and
visualizing results.  The running of tests is incorporated into the
software build system and invoked by executing the provided
`./regression.sh` script, and test results are visualized through a
supplied PHP hypertext preprocessor web page. Cello also has support
for testing invariants during code execution, helping to identify
numerical anomalies and reporting them to the user before they may
otherwise make themselves known, which reduces time debugging and
leads to overall higher software quality.  Defect tracking is
coordinated through `Bugzilla <http://http://www.bugzilla.org/>`_ to
ensure that users have a means to report bugs, and that developers
have the means to obtain the information required to fix them.

Flexibility
-----------

We anticipate new and possibly disruptive developments in hardware, as
well as new developments in numerical methods and parallel data
structures, so flexibility is a high priority in our design and
development.  Object-oriented programming provides a natural
flexibility through polymorphism, encapsulation, and separation of
interface and implementation.  Also, we have adopted an agile
iterative test-driven development approach with continuous
refactoring, which improves flexibility by reducing code complexity
and improving code readability.
