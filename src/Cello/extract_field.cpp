// See LICENSE_CELLO file for license and copyright information

/// @file     extract_field.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-05-21
/// @brief    Program for extracting fields from checkpoint or data dumps
 
#include "main.hpp" 
#include "test.hpp"
#include "disk.hpp"
#include <fstream>
#include <iostream>

Main::Main(CkArgMsg* m)
  : count_exit_(0), count_checkpoint_(0),
    monitor_(NULL), fp_text_(), sync_text_()
{
  PARALLEL_INIT;

  // Read in list of h5 files "file_list" from provided directory
  std::vector<std::string> file_list;
  int n_file;
  const int n_field=m->argc - 2;
  if (m->argc>1) {
    std::string file_name = std::string(m->argv[1]) + "/check.file_list";
    std::ifstream file_stream;
    file_stream.open(file_name);
    file_stream >> n_file;
    file_list.resize(n_file);
    for (int i_f=0; i_f<n_file; i_f++) {
      file_stream >> file_list[i_f];
    }
    file_stream.close();
  } else {
    CkPrintf ("\n\n  Usage: %s <directory> <field_1> [ <field_2> ... ]\n\n",m->argv[0]);
    CkExit();
  }

  const std::string path = std::string(m->argv[1]) + "/";
  int type_data=-1, mox=1,moy=1,moz=1;
  for (int i_field=0; i_field<n_field; i_field++) {
    std::string field = m->argv[i_field+2];
    std::string dataset_name = std::string("field_") + field;
    CkPrintf ("Extracting field %s\n\n",field.c_str());
    // Open output file for field
    FileHdf5 hdf5_out(".",field + ".h5");
    hdf5_out.file_create();

    // Loop over data files, extracting current field and writing
    cello_float * field_out = nullptr;
    for (int i_file=0; i_file<n_file; i_file++) {

      const std::string file_base = file_list[i_file];
      const std::string file_block_list = path + "/" + file_base + ".block_list";
      const std::string file_block_data = file_base+".h5";
      std::ifstream block_stream;
      block_stream.open(file_block_list);
      FileHdf5 hdf5_in(path,file_block_data);
      hdf5_in.file_open();
      int blocking[3];
      int root_size[3];
      int type_scalar;
      int tx,ty,tz;
      hdf5_in.group_chdir("/");
      hdf5_in.group_open();
      hdf5_in.group_read_meta(blocking,"blocking",&type_scalar,&tx,&ty,&tz);
      hdf5_in.group_read_meta(root_size,"root_size",&type_scalar,&tx,&ty,&tz);
      mox = root_size[0];
      moy = root_size[1];
      moz = root_size[2];
      const int nx = root_size[0]/blocking[0];
      const int ny = root_size[1]/blocking[1];
      const int nz = root_size[2]/blocking[2];
      hdf5_in.group_close();
      if (field_out == nullptr) {
        field_out = new cello_float [mox*moy*moz];
      }
      while (block_stream.peek() != EOF) {
        std::string name_block;
        int level;
        block_stream >> name_block;
        block_stream >> level;
        if (name_block.size() > 0 && level == 0) {
          std::string group_name = "/" + name_block;
          hdf5_in.group_chdir(group_name);
          hdf5_in.group_open();

          int is_leaf;
          hdf5_in.group_read_meta(&is_leaf,"is_leaf",&type_scalar,&tx,&ty,&tz);

          if (is_leaf && level >= 0) {
            int array[3];
            hdf5_in.group_read_meta(array,"array",&type_scalar,&tx,&ty,&tz);
            int m4[4];

            // Open input dataset
            hdf5_in.data_open (dataset_name, &type_data,
                               m4,m4+1,m4+2,m4+3);
            int n4[4] = {1,1,1,1},o4[4] = {0,0,0,0};
            n4[0] = m4[0];
            n4[1] = m4[1];
            n4[2] = m4[2];
            const int gx=(m4[0]-nx)/2;
            const int gy=(m4[1]-ny)/2;
            const int gz=(m4[2]-nz)/2;
            const int ox=array[0]*nx;
            const int oy=array[1]*ny;
            const int oz=array[2]*nz;
            hdf5_in. data_slice
              (m4[0],m4[1],m4[2],m4[3],
               n4[0],n4[1],n4[2],n4[3],
               o4[0],o4[1],o4[2],o4[3]);

            hdf5_in.mem_create (m4[0],m4[1],m4[2],m4[0],m4[1],m4[2],0,0,0);
            cello_float * field_in = new cello_float [m4[0]*m4[1]*m4[2]];
            CkPrintf ("   reading %s%s block %s\n",
                      path.c_str(),file_block_data.c_str(),
                      name_block.c_str());
            hdf5_in.data_read (field_in);
            for (int iz=0; iz<nz; iz++) {
              for (int iy=0; iy<ny; iy++) {
                for (int ix=0; ix<nx; ix++) {
                  const int i_in = (ix+gx) + m4[0]*((iy+gy) + m4[1]*(iz+gz));
                  const int i_out = (ix+ox) + mox*((iy+oy) + moy*(iz+oz));
                  field_out[i_out] = field_in[i_in];
                }
              }
            }
            hdf5_in.data_close();
          } // is_leaf
        }
        hdf5_in.group_close();
      }
      hdf5_in.file_close();
    }

    hdf5_out.group_chdir (std::string("/") + field);
    hdf5_out.group_create();
    hdf5_out.mem_create(mox,moy,moz,mox,moy,moz,0,0,0);
    hdf5_out.data_create(field,type_data,moz,moy,mox,1,moz,moy,mox,1);
    CkPrintf ("\n   writing %d %d %d field %s in %s.h5\n\n",
              mox,moy,moz,field.c_str(),field.c_str());
    hdf5_out.data_write(field_out);
    hdf5_out.file_close();
    delete [] field_out;
    field_out = nullptr;

  }

 exit_();

}
PARALLEL_MAIN_END
