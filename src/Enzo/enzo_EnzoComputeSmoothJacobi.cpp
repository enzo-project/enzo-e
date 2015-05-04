// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeSmoothJacobi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoComputeSmoothJacobi class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputeSmoothJacobi::EnzoComputeSmoothJacobi
(
 std::string x_field,
 std::string r_field,
 std::string d_field,
 double weight,
 FieldDescr * field_descr) throw()
  : i_x_ (field_descr->field_id(x_field)),
    i_r_ (field_descr->field_id(r_field)),
    i_d_ (field_descr->field_id(d_field)),
    w_(weight)
{
}

//----------------------------------------------------------------------

void EnzoComputeSmoothJacobi::compute ( Block * block) throw()
{
  Field field = block->data()->field();
  if (field.precision(0) == precision_single) {
    compute_<float>(block);
  } else if (field.precision(0) == precision_double) {
    compute_<double>(block);
  } else if (field.precision(0) == precision_quadruple) {
    compute_<long double>(block);
  }
}

//----------------------------------------------------------------------

template <typename T>
void EnzoComputeSmoothJacobi::compute_(Block * block)
{
  
}

