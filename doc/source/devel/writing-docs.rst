.. include:: ../roles.incl

=====================
Writing Documentation
=====================

The Enzo-E/Cello documentation is written in a plaintext markup language called reStructuredText (reST).
The `Sphinx <https://www.sphinx-doc.org/>`_ documentation generator is used to convert the plaintext files into HTML.

The Sphinx documentation provides an excellent primer on reStructuredText `here <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.
This is a great place to start if you are interesting in contributing to Cello/Enzo-E's documentation.

This page describes some additional extensions that introduce special features for Enzo-E/Cello's documentation.

The ``breathe`` extension
=========================

To reduce some duplication, we make use of `Breathe <https://www.breathe-doc.org/>`_.
It is a Sphinx plugin that we use to selectively integrate documentation generated with `Doxygen <https://www.doxygen.nl/>`_ from docstrings in Enzo-E/Cello's source code.

The ``par`` extension
=====================

We have also written an extension for Sphinx to assist with formatting Cello/Enzo-E's parameters, called ``par``.
At present, this extension provides 3 primary features:

    1. A directive called :rst:dir:`par:parameter` that is to be used when defining a new parameter in the reference section.
       Parameters defined with this directive automatically provide an anchor point to facillitate cross-referencing (i.e. using the :rst:role:`par:param` role).

    2. Pretty-formatting of parameter names.
       This is done by the :rst:dir:`par:parameter` directive, as well as the :rst:role:`par:param` and :rst:role:`par:paramfmt` roles.
       Note, this functionallity supports special formatting of components of parameter names enclosed by angle brackets (e.g. ``:par:paramfmt:`Boundary:<condition>:type``` and ``:par:paramfmt:`Initial:value:<field>``` are respectively rendered as :par:paramfmt:`Boundary:<condition>:type` and :par:paramfmt:`Initial:value:<field>`)

    3. Pretty-formatting of parameter types with the :rst:role:`par:typefmt` role.
       For example, ``:par:typefmt:`list ( float )``` gets rendered as :par:typefmt:`list ( float )`.



Defining a Parameter
--------------------

When defining a new parameter in the reference section, you should use the following directive:

.. rst:directive:: par:parameter

   This directive should be used when defining a new Cello/Enzo-E parameter.

   The argument passed to the directive is the name of the parameter being documented.
   After the line ``.. par:parameter:: ParamName``, you should leave a blank line, and then you should define a list of fields (e.g. ``Summary``, ``Type``, ``Default``, ``Scope``) that provide an overview of the parameter.
   Finally, you generally provide an extended description after the field-list.

   .. rubric:: Example

   It is insightful to consider a concrete example.
   Below we show how to do this for the :par:paramfmt:`Method:heat:alpha` parameter:

   .. code-block:: rst

       .. par:parameter:: Method:heat:alpha

          :Summary:    :s:`Parameter for the forward euler heat equation solver`
          :Type:       :par:typefmt:`float`
          :Default:    :d:`1.0`
          :Scope:      :z:`Enzo`

          :e:`Thermal diffusivity parameter for the heat equation.`

   This snippet will be rendered as the following (cross-referencing to this output has been disabled):

   ..
      Note that below we are passing the ``:noindex:`` option to the directive
      in order to disable the cross-referencing.

   .. par:parameter:: Method:heat:alpha
      :noindex:

      :Summary:    :s:`Parameter for the forward euler heat equation solver`
      :Type:       :par:typefmt:`float`
      :Default:    :d:`1.0`
      :Scope:      :z:`Enzo`

      :e:`Thermal diffusivity parameter for the heat equation.`

   .. note:: Other potential benefits of :rst:dir:`par:parameter`

      Defining parameters with this directive could also facillitate other benefits in the future such as:

         * the automated populating of parameter tables displayed in the User Guide section
         * the automated formatting of the text in the parameter fields
         * the displayed format can be easily updated


Formatting parameters names in text
-----------------------------------

The extension also defines two roles that can be used when mentioning parameters in text.

.. rst:role:: par:param

   This formats the name of a parameter nicely and creates a cross reference to an existing parameter definition.
   For example, ``:par:param:`Mesh:root_size``` renders as :par:param:`Mesh:root_size`.

   You can modify this role's behavior by prefixing the content with a ``~`` or ``!``:

      * when the role's content is prefixed with a ``~``, the rendered link text only shows the final component of the parameter.
        For example, ``:par:param:`~Mesh:root_size``` renders as :par:param:`~Mesh:root_size` (it links to the same place as :par:param:`Mesh:root_size`).

      * when the role's content is prefixed with a ``!``, no link is constructed (e.g. ``:par:param:`!Mesh:root_size``` renders as :par:param:`!Mesh:root_size`).
        This does the same thing as :rst:role:`par:paramfmt`

   .. note::

      Parameter names that also serve as links are intentionally styled different from non-linked parameter names (to provide readers with a visual cue to the difference).
      Additionally, unlinked parameter names rendered with this role are more compactly than the titular parameter name in :rst:dir:`par:parameter` (since this role is commonly used inline)

.. rst:role:: par:paramfmt

   This role is provided as a convenience.
   It effectively aliases the behavior of the :rst:role:`par:param` role when the content is prefixed by ``!`` (i.e. the name is formatted nicely, but no link is created).
   For example, ``:par:paramfmt:`Mesh:root_size``` renders as :par:paramfmt:`Mesh:root_size`.

Formatting parameter types in text
----------------------------------

We have defined the following role to nicely format the names of parameter types.

.. rst:role:: par:typefmt

   Some examples include:

     * ``:par:typefmt:`float``` renders as :par:typefmt:`float`
     * ``:par:typefmt:`list ( float )``` renders as :par:typefmt:`list ( float )`
     * ``:par:typefmt:`list ( float-expr, [ logical-expr, float-expr, [ ... ] ] )``` renders as :par:typefmt:`list ( float-expr, [ logical-expr, float-expr, [ ... ] ] )`
