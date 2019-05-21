
======================================
Enzo-E / Cello Status Report: May 2019
======================================


  "There is a theory called the *Broken Windows Theory*. This
  theory states that if there is a disorder in a neighbourhood for a
  substantial period of time, the mess will only get worse. The
  authors give the example of a broken window in an abandoned
  building. If it isn't replaced within some time, the chances are
  that vandals will break another window. And another. Gradually
  making the mess worse, and therefore making the issue bigger."

  "This theory applies to quality in software development as well. Over
  time the quality of our software may decrease for a variety of
  reasons. For instance, when under the pressure of a deadline we tend
  to do concessions to the quality of our software. Or when some tests
  fail, and no-one acts, it doesn't motivate to write good tests and
  more tests will fail over time. Effectively making the testsuite
  useless. This disorder is also known as *technical debt*."
  
  --Richard Tuin, *Software Development and the Broken Windows Theory*


The Enzo-E / Cello project requires improvement in several project
components.  We organize these into two sections: improving existing
components, and adding new features.

-----------------------------
Improving existing components
-----------------------------

Documentation needs to be restructured, reviewed against existing
code, missing pieces filled in, and updated for correctness.  Testing
needs to be better-organized, new unit and feature tests need to be
added, automated testing and viewable results for multiple
architectures and multiple configurations along multiple dimensions
need to be added, and we cannot continue to ignore failed tests.  Bug
reports need to be recreated to replace those lost in the bugzilla
tracking site hardware crash, and fixing known bugs need to be
actively worked on.  The current codebase needs improvement along
several dimensions, including using modern C++ features, improving
code readability and reliability, reviewing and improving existing
API's, and refactoring the more complex and error-prone code sections.
Performance also needs to be more in the forefront, with new
performance and scalability tests added to the automated testing
framework, and identifying and creating bug reports for performance
and scaling issues found.  User parameters should be reviewed for
organization, naming, usefulness, and usability.  The problems in the
input directory need to be reviewed, renaming existing tests in a more
structured and organized manner, outdated parameters brought
up-to-date, and commented.  Lastly, the project itself needs to be
better organized, with project documents covering different components
of the project, and design documents for proposed added features and
physics capabilities.

We propose the following documents for organizing these different
aspects of the project.  These documents are concerned with 
more general organization rather than specific details.

  - ``project-plan``: Project planning
  - ``project-doc``: Documentation
  - ``project-test``: Automated testing
  - ``project-params``: User parameters
  - ``project-input``: Input files
  - ``project-code``: Code quality improvements
  - ``project-perf``: Performance optimization
  - ``project-bug``: Bug reporting and debugging

-------------------    
Adding new features
-------------------    

Many new features are required and desired, some related to physics
capabilities, some related to performance, and some to basic
infrastructure for highly-scalable AMR astrophysics.

New features involving basic physics and underlying technical
capabilities are required, including flux-correction for AMR,
scalable gravity for AMR, scalable and reliable I/O.

Performance features include support for accelerators, which will
likely involve other dependent capabilities, such as
dynamic load balancing that keeps data locality, agglomerating
neighboring Block's into larger components to reduce ghost cells
and improve accelerator efficiency.

We propose the following documents for documenting the
design of these proposed updates.  These documents are concerned
more with mid- and low-level design details and issues than
the higher-level project documents.

  - ``design-accel``: Accelerator usage
  - ``design-flux``: Flux correction
  - ``design-interp``: Conservative interpolation
  - ``design-units``: 
  - ``design-io``:
  - ``design-particle``: Complete missing particle capabilities
  - ``design-array``
  - ``design-gravity``:
  - ``design-active``:
  - ``design-ray``:
  - ``design-timestep``
  - ``design-grackle``:
  - ``design-balance``:
  - ``design-restart``:
  - ``design-smp``


