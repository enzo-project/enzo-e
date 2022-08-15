#define RHOUNIT 1.674e-24
#define TUNIT 3.154e13
#define LUNIT (6.6e3*3.086e18)
#define VOLUNIT LUNIT*LUNIT*LUNIT
#define FUNIT (1.0/(TUNIT*LUNIT*LUNIT))
#define NUNIT (1.0/(VOLUNIT))
#define CLIGHT 0.1*(29979245800.0)
#define DX_CGS LUNIT/128.0
#define DT_CGS DX_CGS/(3.0*CLIGHT)
#define CELL_VOLUME_CGS (DX_CGS*DX_CGS*DX_CGS)

Boundary {
     list = ["N0", "N1", "N2", "Fx0", "Fx1", "Fx2"];
     Fx0 {
       type = "inflow";
       field_list = "flux_x_0";
       value = [0.447*1e6 / FUNIT, x <= 0.0, 1e-16/FUNIT];
     }
     Fx1 {
       type = "inflow";
       field_list = "flux_x_1";
       value = [0.494*1e6 / FUNIT, x <= 0.0, 1e-16 / FUNIT];
     }
     Fx2 {
       type = "inflow";
       field_list = "flux_x_2";
       value = [0.59*1e6 / FUNIT, x <= 0.0, 1e-16 / FUNIT];
     }
     N0 {
       type = "inflow";
       field_list = "photon_density_0";
      value = [2.235e48*DT_CGS/CELL_VOLUME_CGS / NUNIT, x <= 0.0, 1e-16 / NUNIT];
     }
     N1 {
       type = "inflow";
       field_list = "photon_density_1";
       value = [2.47e48*DT_CGS/CELL_VOLUME_CGS / NUNIT, x <= 0.0, 1e-16 / NUNIT];
     }
     N2 {
       type = "inflow";
       field_list = "photon_density_2";
       value = [2.95e48*DT_CGS/CELL_VOLUME_CGS / NUNIT, x <= 0.0, 1e-16 / NUNIT];
     }
 }

 Domain {
     lower = [ 0.000000000000000, 0.000000000000000, 0.000000000000000 ];
     upper = [ 1.000000000000000, 1.000000000000000, 1.000000000000000 ];
 }

 Field {
     alignment = 8;
     gamma = 1.400000000000000;
     ghost_depth = 4;
     list = [ "density", "internal_energy", "total_energy", "velocity_x", "velocity_y", "velocity_z", "pressure", "temperature", "photon_density", "flux_x", "flux_y", "flux_z", "photon_density_0", "flux_x_0", "flux_y_0", "flux_z_0", "photon_density_1", "flux_x_1", "flux_y_1", "flux_z_1", "photon_density_2", "flux_x_2", "flux_y_2", "flux_z_2", "metal_density", "HI_density", "HII_density", "HeI_density", "HeII_density", "HeIII_density", "e_density", "RT_heating_rate", "RT_HI_ionization_rate", "RT_HeI_ionization_rate", "RT_HeII_ionization_rate" ];
 }

 Group {
     color {
         field_list = [ "HI_density", "HII_density", "HeI_density", "HeII_density", "HeIII_density", "e_density" ];
     };
     derived {
         field_list = [ "temperature", "pressure" ];
     };
     list = [ "color", "derived" ];
 }

 Initial {
     list = [ "value" ];
     value {  
         density = [ 200.0*2e-4 * 1.674e-24/RHOUNIT,
                         sqrt((x - 0.75757575756)*(x - 0.757575757576) + (y - 0.5)*(y - 0.5) + (z - 0.5)*(z - 0.5)) <= 0.1212121212121, 
                   2e-4 * 1.674e-24/RHOUNIT];
         HI_density = [ 200.0 * 2e-4 * 1.674e-24/RHOUNIT,
                         sqrt((x - 0.75757575756)*(x - 0.757575757576) + (y - 0.5)*(y - 0.5) + (z - 0.5)*(z - 0.5)) <= 0.1212121212121, 
                   2e-4 * 1.674e-24/RHOUNIT];
         temperature = [ 40.0,
                         sqrt((x - 0.75757575756)*(x - 0.757575757576) + (y - 0.5)*(y - 0.5) + (z - 0.5)*(z - 0.5)) <= 0.1212121212121, 
                       8000.0];

         HII_density = 0.0;
         HeIII_density = 0.0;
         HeII_density = 0.0;
         HeI_density = 0.0;
         e_density = 0.0;
  
         photon_density_0 = 1e-16/NUNIT;
         photon_density_1 = 1e-16/NUNIT;
         photon_density_2 = 1e-16/NUNIT;

         flux_x_0 = 1e-16/FUNIT;
         flux_x_1 = 1e-16/FUNIT;
         flux_x_2 = 1e-16/FUNIT;

         flux_y_0 = 1e-16/FUNIT;
         flux_y_1 = 1e-16/FUNIT;
         flux_y_2 = 1e-16/FUNIT;

         flux_z_0 = 1e-16/FUNIT;
         flux_z_1 = 1e-16/FUNIT;
         flux_z_2 = 1e-16/FUNIT;
     }
 }

 Mesh {
     root_blocks = [4,4,4];
     root_rank = 3;
     root_size = [128, 128, 128]; 
 }

 Method {
     list = [ "null", "ramses_rt", "grackle", "ppm" ];
     grackle {
         CaseBRecombination = 1;
         UVbackground = 0;
         data_file = "fuck";
         metal_cooling = 0;
         primordial_chemistry = 1;
         self_shielding_method = 0;
         use_cooling_timestep = false;
         use_radiative_transfer = 1;
         with_radiative_cooling = 1;
     };
     ppm {
       courant = 0.8;
       diffusion = true;
       flattening = 3;
       steepening = true;
       dual_energy = true;
     };
     ramses_rt {
         N_groups = 3;
         bin_lower = [ 13.60000000000000, 24.59000000000000, 54.42000000000000 ];
         bin_upper = [ 24.59000000000000, 54.42000000000000, 100.0000000000000 ];
         clight_frac = 0.100000000000000;
         radiation_spectrum = "none";
         recombination_radiation = false;
         average_global_quantities = false;
     };
     null {dt = 1e-4;}
 }

 Output {
     list = [ "N", "Fx", "Fy", "Fz", "de", "vx", "vy", "e", "HI", "HII", "hdf5", "check" ];
     de {
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "density" ];
         image_log = false;
         image_size = [ 512, 512 ];
         name = [ "de-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     vx {
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "velocity_x" ];
         image_log = false;
         image_size = [ 512, 512 ];
         name = [ "vx-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     vy {
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "velocity_y" ];
         image_log = false;
         image_size = [ 512, 512 ];
         name = [ "vy-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     Fx {
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "flux_x" ];
         image_log = false;
         image_size = [ 512, 512 ];
         name = [ "RT_Fx-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     Fy {
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "flux_y" ];
         image_log = false;
         image_size = [ 512, 512 ];
         name = [ "RT_Fy-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     Fz {
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "flux_z" ];
         image_log = false;
         image_size = [ 512, 512 ];
         name = [ "RT_Fz-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     HI {
         colormap = [ 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 ];
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "HI_density" ];
         image_log = false;
         image_size = [ 512, 512 ];
         image_type = "data";
         name = [ "RT_HI-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     HII {
         colormap = [ 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 ];
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "HII_density" ];
         image_log = false;
         image_size = [ 512, 512 ];
         image_type = "data";
         name = [ "RT_HII-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     N {
         colormap = [ 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 ];
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "photon_density" ];
         image_log = false;
         image_size = [ 512, 512 ];
         image_type = "data";
         name = [ "RT_N-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     check {
         dir = [ "checkpoint-%06d", "count" ];
         schedule {
             start = 10100;
             step = 10000;
             var = "cycle";
         };
         type = "checkpoint";
     };
     e {
         colormap = [ 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 ];
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "e_density" ];
         image_log = false;
         image_size = [ 512, 512 ];
         image_type = "data";
         name = [ "RT_e-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
     hdf5 {
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "photon_density", "flux_x", "flux_y", "flux_z", "density", "velocity_x", "velocity_y", "velocity_z", "HI_density", "HII_density", "HeI_density", "HeII_density", "HeIII_density", "e_density", "temperature", "pressure", "internal_energy", "total_energy", "RT_heating_rate", "RT_HI_ionization_rate", "RT_HeI_ionization_rate", "RT_HeII_ionization_rate", "photon_density_0", "photon_density_1", "photon_density_2" ];
         name = [ "data-%03d-%02d.h5", "count", "proc" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "data";
     };
     mesh {
         colormap = [ 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 ];
         dir = [ "ILIEV7-%06d", "cycle" ];
         field_list = [ "photon_density" ];
         image_log = false;
         image_size = [ 512, 512 ];
         image_type = "mesh";
         name = [ "RT_N_mesh-%06d.png", "cycle" ];
         schedule {
             start = 0;
             step = 100;
             var = "cycle";
         };
         type = "image";
     };
 }

 Particle {
     list = [ "star" ];
     star {
         attributes = [ "x", "double", "y", "double", "z", "double", "vx", "double", "vy", "double", "vz", "double", "ax", "double", "ay", "double", "az", "double", "id", "double", "mass", "double", "is_copy", "int64", "creation_time", "double", "lifetime", "double", "metal_fraction", "double" ];
         groups = [ "is_gravitating" ];
         mass_is_mass = true;
         position = [ "x", "y", "z" ];
         velocity = [ "vx", "vy", "vz" ];
     };
 }

 Stopping {
     cycle = 100000;
     time = 50.0000000000000;
 }

 Units {
     density = RHOUNIT;
     length = LUNIT;
     time = TUNIT;
 }
