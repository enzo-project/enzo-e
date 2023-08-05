***********************************
Parameter Parsing and Configuration
***********************************

*[ This page is under development ]*

========
Overview
========

Most of the parsing work is handled by the :cpp:class:`!Parameters` class from the Cello layer and stored internally.
All parameters are parsed after startup on a single process and then the parsed information is communicated between processes.

We are currently in the process of migrating between approaches for accessing parsed parameters (and configuring.
Given the large-scope of this transition, we decided to gradually migrate between the approaches, rather than try to do it one shot.

* At the time of writing, we have migrated a large number of :cpp:class:`!Method` classes.
  
* Our intention is to eventually migrate as much as possible to this new approach (e.g. all :cpp:class:`!Method` classes, all :cpp:class:`!Physics` classes, all :cpp:class:`!Initial` classes, etc.), to stop ultimately get rid of the :cpp:class:`!Config` and :cpp:class:`!EnzoConfig` classes.

First, we briefly review properties of an Enzo-E/Cello parameters and talk about the idea of a "parameter-path"
Next, we talk through an example of how you might add a new parameter using the new approach.
Then, we provide some more details about the design of the new-approach.
Finally, we briefly discuss the older approach (and some of its problems shortcomings).

===================================
Parameter Files and Parameter-Paths
===================================

.. todo::

   Consider merging the content of this section into :ref:`parameter-file`\.

As explained in :ref:`parameter-file`\, Enzo-E/Cello makes use of a hierarichal parameter file (this documentation assumes you are familiar with the basics from that section).

Essentially, parameters are organized within groups and (possibly subgroups).
In other words, the parameters can be thought of as leaf nodes in a tree-hierarchy of "groups".
The following snippet illustrates the organization of parameters from a hypothetical parameter file:

 ::

  Method {
     list = [ "mhd_vlct", "grackle"];

     mhd_vlct {
        # other assorted parameters...
        courant = 0.3;
     }

     grackle {
        # other assorted parameters...
        courant = 0.4;
     }
  }

The organization of parameters in a group hierarchy is analogous to the organization of files in a directory hierarchy.
Continuing this analogy, we have devised shorthand for naming parameters in the documentation (and throughout the codebase) that is similar to a file path.
One can think of these names as a "parameter-path".

- documentation refering to the parameters in the above snippet would mention :par:paramfmt:`Method:list`, :par:paramfmt:`Method:mhd_vlct:courant`, or :par:paramfmt:`Method:grackle:courant`.

- in the codebase, one could see strings refering to ``"Method:list"``, ``"Method:mhd_vlct:courant"``, ``"Method:grackle:courant"`` (in practice, you may not see these particular strings).

- to be clear, a parameter-path for a given parameter lists the names of ancestor "groups", separated by colons, and lists the name of the parameter at the end (i.e. the string that directly precedes an assignment).

=========================================
How to add a new parameter (new approach)
=========================================

As mentioned above, the nitty-gritty details of parsing are handled by Enzo-E automatically.
Values associated with parameter names can be queried by invoking methods on an instance of the :cpp:class:`!ParameterAccessor` class.
Instances of this class are commonly passed to the constructors of classes that inherit from Cello class-hierarchy.


Let’s walk through an example where we want to introduce a new parameter to :cpp:class:`!EnzoMethodHeat`. 
Suppose we want to add a new parameter to the called :par:param:`!my_param`.
The full name of this parameter would be :par:param:`!Method:heat:my_param`.

The steps are as-follows:

1. Introduce a new member-variable (aka an attribute) to :cpp:class:`!EnzoMethodHeat` (in the ``EnzoMethodHeat.hpp`` file).
   For the sake of example, let's imagine that we want to directly store the value specified in the parameter-file in a member-variable (a.k.a. an attribute) named  ``my_param_``.

   - The convention is to declare all member-variables as ``private`` or ``protected`` (if the value of that attribute is needed outside of the class, you should define a ``public`` accessor-function).

   - Relatedly, the names of all ``private`` & ``protected`` member-variables or member-functions should generally be suffixed with an underscore.
     An undersore should **NEVER** be the first character in the name of a member-variable or member-function.

   - **NOTE:** the value of the parameter doesn't necessarily need to initialize a variable with a matching name (or type), this is just a convenience in this example (although, it does make the code a little easier to follow)

2. Modify the pup routine of :cpp:class:`!EnzoMethodHeat` and the ``PUP::able`` migration constructor to properly handle the newly added member-variable

3. Modify the main constructor of :cpp:class:`!EnzoMethodHeat` to initialize ``my_param_`` based on the value parsed from the parameter file.
   The constructor of :cpp:class:`!EnzoMethodHeat` is passed a copy of an instance of :cpp:class:`!ParameterAccessor`, in an argument ``p``.
   To access the value specified in the parameter-file, you would probably invoke one of the following expressions (based on the expected type of the parameter):

   - ``p.value_logical("my_param", false)`` if the parameter is expected to specify a boolean value and if it defaults to a value of ``false`` when the parameter is not specified.

   - ``p.value_integer("my_param", 7)`` if the parameter is expected to be an integer and if it defaults to a value of ``7`` when the parameter is not specified.

   - ``p.value_float("my_param", 2.0)`` if the parameter is expected to be a floating-point value and if it defaults to a value of ``2.0`` when the parameter is not specified.

   - ``p.value_string("my_param", "N/A")`` if the parameter is expected to be a string and if it defaults to a value of ``"N/A"`` when the parameter is not specified.

   In all of these cases, ``"my_param"`` is internally understood to be the last part of a full parameter-name.
   We will discuss this down below in more detail, but for the sake of this example, it's generally expanded to ``"Method:heat:my_param"``.
   Alternative logic may be needed to the above expressions in slightly more sophisticated cases (for example if the parameter expects a list of values or if you want to abort the program if the parameter can't be found).

.. note::

    In some cases, parameter-values may be retrieved from the parameter-file in a static factory method, typically called ``from_parameters``, rather than in a constructor.

    There are sometimes cases where the use of factory-method could be very useful (that should usually justified in a comment).
    But factory-methods are most commonly used in cases where it simplified the transition from the older parameter-parsing approach to this newer approach.


==============================
Design Overview (new approach)
==============================

Our new approach revolves around the usage of the :cpp:class:`!ParameterAccessor` class for accessing/querying parameters stored in an instance of the :cpp:class:`Parameters` class.
Instances of the :cpp:class:`!ParameterAccessor` class are light-weight and are expected to have a short-lifetime (akin to :cpp:class:`!Field` or :cpp:class:`!Particle`).

**As illustrated above, instances of the** :cpp:class:`!ParameterAccessor` **class are expected to be passed to the constructor (or factory-method) of classes that inherit from Cello class-hierarchy.**

The main feature of the :cpp:class:`!ParameterAccessor` class is that it provides methods for querying/accessing parameters with parameter-paths that share a common root.

- The root parameter-path is specified during the construction of a :cpp:class:`!ParameterAccessor` instance and cannot be changed over the lifetime of the instance.

  - The immutable nature of the root parameter-path this property is a feature: whenever a :cpp:class:`!ParameterAccessor` instance is passed to a function, you ALWAYS know that the root parameter-path is unchanged (without needing to check the helper function's implementation).

  - If a developer is ever tempted to mutate the root-path, they should just initialize a new :cpp:class:`!ParameterAccessor` (since the instances are lightweight)

  - The root path can be queried with :cpp:expr:`ParameterAccessor::get_root_parpath()`

- When a string is passed to one of the accessor methods, that string is internally appended to the root pararameter-path and the result represents the full name of the queried parameter.
  (You can think of this string as specifying the relative path to the parameter).
  You can use :cpp:expr:`ParameterAccessor::full_name(s)` to see the full parameter name that a string, ``s``, maps to.


Why do we even need :cpp:class:`!ParameterAccessor`?
----------------------------------------------------

To motivate the existence of the :cpp:class:`!ParameterAccessor` class, it's useful to consider alternative approaches.
The most obvious option is to simply pass instances of the :cpp:class:`!Parameters` class to constructors and static-factory methods of other classes.

To flesh out this alternative case more, let's consider the following snippet of a hypothetical parameter file.

.. code-block::

       Method {
         list = [
            "output", # ...
         ];

         output {
           file_name = # ...
           all_fields = true;
           all_particles = true;
           # ...

         }
       }

This particular snippet can easily be parsed if we pass a reference to the :cpp:class:`!Parameters` object to the static factory method of :cpp:class:`!MethodOutput`.
An example code block is included here, to show what that the initialization might look like:

.. code-block:: c++

    // NOTE: MethodOutput is a special case. At least historically, its factory
    // method needed to accept an argument other than just the parameters

    MethodOutput* MethodOutput::from_parameters(/* ... */, Parameters &p)
    {
      std::string file_name = p.value_string("Method:output:file_name", "");
      bool all_fields = p.value_logical("Method:output:all_fields", false);
      bool all_particles = p.value_logical("Method:output:all_particles", false);

      // parse other parameters...

      return new MethodOutput( /* pass file_name, all_fields, all_particles,
                                  and other values */ );
    }

There is nothing wrong with the above snippet, and it will work in a lot of cases.
However, we will encounter issues, when we want to setup a simulation that makes use of multiple :cpp:class:`!MethodOutput` instances.
To illustrate how this is done in Enzo-E, see the following snippet from a hypothetical parameter file:

.. COMMENT-BLOCK

    I'm a little tempted to show a schedule grouping, but I'm not sure
    it adds much (since parsing of schedule objects is handled
    separately)...

    Additionally, should we be using a real method in these examples?
    What if parameter choices change in the future?

.. code-block::

        Method {
          list = [
             "output_field", "output_particle", # ...
          ];

          output_field {
            file_name = # ...
            all_fields = true;
            all_particles = false;
            # parse other parameters ...

            type = "output";
          }

          output_particle {
            file_name = # ...
            all_fields = false;
            all_particles = true;
            # parse other parameters ...

            type = "output";
          }

        }

As you can see from the above snippet, a parameter-subgroup carrying the configuration of the :cpp:class:`!MethodOutput` instance is no longer called ``"output"``.

- now 2 :cpp:class:`!MethodOutput` instances should be initialized, using the configuration from the parameter-subgroups called ``"output_field"`` and ``"output_particle"``.

- Cello/Enzo-E determines the :cpp:class:`!Method` subclass that a given parameter-subgroup, :par:param:`!Method:<subgroup>`, is meant to describe based on the value of :par:param:`!Method:<subgroup>:type`.
  In both above subgroups, we have specified the type as ``"output"``.
  (In the common case where :par:param:`!Method:<subgroup>:type` is omitted, the type parameter defaults to the string-value of :par:param:`!<subgroup>`)

**Importantly,** the absolute paths of the parameters that are used to initialize the :cpp:class:`!MethodOutput` instances are different in the second parameter file compared to the first.
The main difference is in the the root-path to the subgroup.

To gracefully handle both scenarios, we now make use of the of the :cpp:class:`!ParameterAccessor` class.
A code snippet using our new approach is shown below:

.. code-block:: c++

    // NOTE: MethodOutput is a special case. At least historically, its factory
    // method needed to accept an argument other than just the parameters

    MethodOutput* MethodOutput::from_parameters(/* ... */, ParameterAccessor &p)
    {
      std::string file_name = p.value_string("file_name", "");
      bool all_fields = p.value_logical("all_fields", false);
      bool all_particles = p.value_logical("all_particles", false);

      // parse other parameters...

      return new MethodOutput( /* pass file_name, all_fields, all_particles,
                                  and other values */ );
    }

.. note::

   Historically, the :cpp:class:`!Parameters` class has also had the capability to track a common root-path.
   However, the code was not very explicit about whether that capability was being used or not (although, most of the time you could safely assume that the feature wasn't being)

   It's our intention to eventually remove this capability from the :cpp:class:`!Parameters` class, since the :cpp:class:`!ParameterAccessor` class can be used for the same purpose (and it's more explict)


.. note::

   The main disadvantage of this approach is that we no longer specify the full, absolute parameter names, when accessing the values.
   However, this is mostly unavoidable if we want to gracefully accomodate initialization of multiple instances of the same :cpp:class:`!Method` subclass.
   Hopefully, this page of documentation will help to offset this disadvantage.

   The *only* other alternative is have :cpp:class:`!ParameterAccessor` instances "auto-magically" redirect absolute parameter-paths, but I think that will generally be more confusing. 

Other Benefits
--------------

Another benefit to the :cpp:expr:`ParameterAccessor::wrapped_Parameters_ref()` class is that it restricts access to parameters within the associated root-path.
This discourages the design of classes that are configured by parameters scattered throughout the parameter file.

In rare cases (e.g. during refactoring when we convert a previously Method-specific parameter to a Physics parameter and want to retain backwards compatability), exceptions need to be made.
Thus, an "escape-hatch" is provided to directly access the wrapped Parameters object: call the :cpp:expr:`ParameterAccessor::wrapped_Parameters_ref()` method.
Please, avoid using this "escape-hatch" unless it's truly necessary.

.. todo::

   Consider extending the analogy between a parameter-path and a file path.
   For example, one could imagine interpretting a path that begins with a ``:`` as an absolute parameter-path and all other strings as relative parameter-paths.

   This would probably streamline the documentation to some degree.

   If we were to do that, we would need to modify the code to recognize this convention.
   We would probably also want to modify the various parameter-accessor methods of the :cpp:class:`!ParameterAccessor` to continue to restrict access to parameters within the common root-path that a :cpp:class:`!ParameterAccessor` is configured with.

=================
Historic Approach
=================

Historically, all parameters were parsed shortly after startup and then the results were stored as variables in the :cpp:class:`!Config` and :cpp:class:`!EnzoConfig` classes.
However, this approach had a number of warts:

- Adding a new parameter "properly" was laborious. Let’s imagine that we want to add  a parameter, ``<param>``, to class ``<Klass>``.
  This class might be a subclass of :cpp:class:`Method`, :cpp:class:`!Initial`, :cpp:class:`!Physics`, etc.
  To add this new parameter, we need to

  (i) define a new member-variable (aka attribute) to :cpp:class:`!EnzoConfig` to hold the value of the parameter

  (ii) ensure that the new member-variable of :cpp:class:`!EnzoConfig` is properly serialized

  (iii) add the logic to retrieve the value associated with the parameter from the :cpp:class:`!Parameters` object and store that value in the newly defined member variable of :cpp:class:`!EnzoConfig`

  (iv) modify the line of code where :cpp:class:`!EnzoProblem` calls the constructor of ``<Klass>``, in order to pass the parameter-value stored in the newly-defined member-variable of :cpp:class:`!EnzoConfig`.

  (v) add a newly-defined member variable on ``<Klass>`` in order to store the value of the parsed parameter

  (vi) ensure that the new member-variable of :cpp:class:`!EnzoConfig` is properly serialized

  (vii) modify the primary constructor of ``<Klass>`` to actually initialize the new member-variable

- Because of how laborious this is, developers have a tendency to just skip the last for steps and access the attributes of the global :cpp:class:`!EnzoConfig` instance.
  This has all the short-comings of global variables (it makes things hard to refactor)


- If you want to do error-checking of the parameter-values, it’s not always clear where to do that (within :cpp:class:`!EnzoConfig` vs within the constructor of the class that uses the parameter)


- complications arise if multiple instances of a class can be initialized with different configurations.


Our new practice takes inspiration from Athena++.
Essentially, the new approach's intention is to have every :cpp:class:`!Method`\/:cpp:class:`!Initial`\/:cpp:class:`!Physics` class just directly access the needed values from the parameter file.
We skip the whole step of storing the parsed values in an :cpp:class:`!Config` or :cpp:class:`!EnzoConfig` instance and then forwarding those values.
We essentially "cut out the middleman".
