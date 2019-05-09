

Parameter File Example
======================

Below we describe the different sections of the "HelloWorld.in" input
file we ran in the previous section:


Comments
--------

The example input file begins with a comment header.  Comments in
Cello parameter files begin with ``#`` :

  ::

    ### Problem: Hello World
    ### Summary: Evolve high-pressure region shaped as an "E"
    ###    Date: 2011-10-27
    ###  Author: James Bordner (jobordner@ucsd.edu)

Parameters are organized into "groups", which may be nested.  Below
we look in turn at each group in the HelloWorld.in sample parameter file.

Domain Group
------------

  ::

    Domain { 
       lower = [ 0.0, 0.0 ];
       upper = [ 1.0, 1.0 ];
    } 

Mesh Group
----------

  ::

   Mesh { 
      root_size      = [160, 160];
      root_blocks    = [  1,   1 ];
   }
    
Field Group
-----------

  ::

    Field {
    
       ghost_depth  = 3;
       courant = 0.8;
    
       fields = [ 
          "density",
          "velocity_x",
          "velocity_y",
          "total_energy",
          "internal_energy"
        ];
    }

Method Group
------------

  ::

    Method {
    
       sequence = [ "ppm" ];

       ppm {
          diffusion   = true;
          flattening  = 3;
          steepening  = true;
          dual_energy = false;
       }
    }

Physics Group
-------------

  ::

    Physics {
    
       dimensions = 2;
       gamma = 1.4;
    
    }

Initial Group
-------------

  ::

   Initial {
       density       { value = [ 1.0, 
       (0.35 < x && x < 0.40 && 0.25 < y && y < 0.70) ||
       (0.35 < x && x < 0.60 && 0.25 < y && y < 0.30) ||
       (0.35 < x && x < 0.60 && 0.45 < y && y < 0.50) ||
       (0.35 < x && x < 0.60 && 0.65 < y && y < 0.70),
                                 0.125 ]; };
       total_energy  { value = [ 1.0  / (0.4 * 1.0) ,
       (0.35 < x && x < 0.40 && 0.25 < y && y < 0.75) ||
       (0.35 < x && x < 0.60 && 0.25 < y && y < 0.30) ||
       (0.35 < x && x < 0.60 && 0.45 < y && y < 0.55) ||
       (0.35 < x && x < 0.60 && 0.65 < y && y < 0.75),
	                         0.14 / (0.4 * 0.1)]; };
       velocity_x    { value = [0.0]; };
       velocity_y    { value = [0.0]; }
   }

Boundary Group
--------------

  ::

    Boundary { type = "reflecting" }


Stopping Group
--------------

  ::

    Stopping {
       cycle = 500;
    }

Output Group
------------

  ::

   Output { 

      file_groups = ["cycle_step"];

      cycle_step {
         field_list = ["density"];
         type     = "image";
         name       = ["E-%03d.png","cycle"];
         schedule = ["cycle","interval", 10];
         colormap_alpha = [0.00, 0.00, 1.00, 0.0,
                           0.00, 1.00, 0.00, 0.33,
                           1.00, 1.00, 0.00, 0.66, 
                           1.00, 0.00, 0.00, 1.0];
      }
   }
