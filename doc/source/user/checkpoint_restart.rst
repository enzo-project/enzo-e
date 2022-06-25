********************
Checkpoint / restart
********************

Enzo-E can perform checkpoint dumps to save the current state, and can
read them at initialization to continue a simulation where it left
off. Setting up checkpointing involves including methods in the
`Method` parameter group, and provided a schedule for them. Restarting
requires adding a small number of parameters to a parameter file to
turn on restart, and where the checkpoint files are located.


==========
Checkpoint
==========

.. code-block::

   Method {
     ...
     list = [ "order_morton", "check", ... ];
     order_morton {
        schedule { var = "cycle";  start = 5;  step = 5;  }
     }
     check {
        schedule { var = "cycle";  start = 5;  step = 5;  }
        dir = [ "Check-%02d", "cycle" ];
        num_files = 4;
        ordering = "order_morton";
     }
     ...
   }

This example writes checkpoint dumps every 5 cycles, starting with the 5th
cycle. See the input/Checkpoint/test_cosmo-checkpoint.in parameter
file for a working example of writing checkpoint directories.

=======
Restart
=======

To restart from a checkpoint directory, edit a copy of the original parameter
file to add two parameters, `Initial:restart` and `Initial:restart_dir`.
For example, to restart the problem in the previous "Checkpoint" section
at cycle 50, one would add

.. code-block::

   Initial {

      restart = true;
      restart_dir = "Check-50";
   }

