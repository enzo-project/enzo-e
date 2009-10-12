//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      test_parameters.cpp
 * @brief     Program implementing unit tests for the Parameters class
 * @author    James Bordner
 * @date      Thu Feb 21 16:04:03 PST 2008
 * 
 * $Id: test_parameters.cpp 715 2009-07-08 23:48:09Z bordner $
 * 
 *********************************************************************
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "test.hpp"
#include "error.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{

  unit_class ("Parameters");

  unit_open();

  //----------------------------------------------------------------------
  // test parameter
  //----------------------------------------------------------------------

  Parameters * parameters = new Parameters;

  unit_func("set_group()");

  printf ("parameters->get_group() = %s\n",parameters->get_group().c_str());
  parameters->set_group("Group 1");

  unit_assert(parameters->get_group() == "Group 1");
  unit_assert(parameters->get_subgroup() == "");

  unit_func("read()");
  FILE * fp = fopen ("test.in","r");

  parameters->read ( fp );

  parameters->set_group("Group4");

  parameters->set_subgroup("subgroup3");

  double value;

  value = parameters->value_scalar("param_scalar_expr1");
  unit_assert(value == 30.625);

  value = parameters->value_scalar("param_scalar_expr2");
  unit_assert(value == 18.375);

  value = parameters->value_scalar("param_scalar_expr3");
  unit_assert(value == 150.0625);

  value = parameters->value_scalar("param_scalar_expr4");
  unit_assert(value == 4.0000000000);

  parameters->set_subgroup("subgroup4");

  value = parameters->value_scalar("param_scalar_expr1");
  unit_assert(value == 36.750);

  value = parameters->value_scalar("param_scalar_expr2");
  unit_assert(value == 67.375);

  value = parameters->value_scalar("param_scalar_expr3");
  unit_assert(value == -30.625);



  double x[] = { 1, 2, 3};
  double y[] = {5 , 4, 3};
  double z[] = {8, 9, 10};
  double t[] = {-1, 2, -7};
  double values[] = {0,0,0};
  double deflts[] = {-1,-2,-3};

  parameters->set_group("Group5");
  parameters->set_subgroup("subgroup5");
  
  //  param_scalar_expr1 = x;
  parameters->evaluate_scalar("param_scalar_expr1",
			      3,values,deflts,x,y,z,t);
  unit_assert (values[0]==x[0]);
  unit_assert (values[1]==x[1]);
  unit_assert (values[2]==x[2]);

  
  //  param_scalar_expr2 = x+3.0;
  parameters->evaluate_scalar("param_scalar_expr2",
			      3,values,deflts,x,y,z,t);
  unit_assert (values[0]==x[0]+3.0);
  unit_assert (values[1]==x[1]+3.0);
  unit_assert (values[2]==x[2]+3.0);

  //  param_scalar_expr3 = x+y+z+t;
  parameters->evaluate_scalar("param_scalar_expr3",
			      3,values,deflts,x,y,z,t);
  unit_assert (values[0]==x[0]+y[0]+z[0]+t[0]);
  unit_assert (values[1]==x[1]+y[1]+z[1]+t[1]);
  unit_assert (values[2]==x[2]+y[2]+z[2]+t[2]);

  parameters->set_group("Group6");
  parameters->set_subgroup("subgroup6");

  //  param_scalar_expr1 = sin(x);
  parameters->evaluate_scalar("param_scalar_expr1",
			      3,values,deflts,x,y,z,t);
  unit_assert (values[0]==sin(x[0]));
  unit_assert (values[1]==sin(x[1]));
  unit_assert (values[2]==sin(x[2]));

  //  param_scalar_expr2 = atan(y/3.0+2*t);
  parameters->evaluate_scalar("param_scalar_expr2",
			      3,values,deflts,x,y,z,t);
  unit_assert (values[0]==atan(y[0]/3.0+2*t[0]));
  unit_assert (values[1]==atan(y[1]/3.0+2*t[1]));
  unit_assert (values[2]==atan(y[2]/3.0+2*t[2]));

  unit_close();

  delete parameters;

}
