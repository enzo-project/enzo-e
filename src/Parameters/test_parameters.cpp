/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * SYNOPSIS:
 *
 *    
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

/** 
 *********************************************************************
 *
 * @file      test_parameters.cpp
 * @brief     Program implementing unit tests for the Parameters class
 * @author    James Bordner
 * @date      Thu Feb 21 16:04:03 PST 2008
 * 
 * $Id$
 * 
 *********************************************************************
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error.hpp"
#include "test.hpp"
#include "parameters.hpp"

#define CLOSE(a,b) ((((a) - (b)) / (fabs(a) + fabs(b))) < 1e-16)

int main(int argc, char **argv)
{

  unit_class ("Parameters");

  unit_open();

  //----------------------------------------------------------------------
  // test parameter
  //----------------------------------------------------------------------

  Parameters * parameters = new Parameters;

  // Read

  unit_func("read");
  FILE * fpin = fopen ("test.in","r");
  parameters->read ( fpin );

  // Write

  unit_func("write");
  FILE * fpout = fopen ("test.out","w");
  parameters->write ( fpout );
  unit_assert(0); //FAILS

  // set_group()

  unit_func("set_group");
  parameters->set_group("Group");
  unit_assert(parameters->get_group() == "Group");

  // set_subgroup()
  unit_func("set_subgroup");
  parameters->set_subgroup("subgroup_1");
  unit_assert(parameters->get_subgroup() == "subgroup_1");
  parameters->set_subgroup("subgroup_2");
  unit_assert(parameters->get_subgroup() == "subgroup_2");

  // set_group() without set_subgroup() clears subgroup to ""

  unit_func("set_group");
  parameters->set_group("Group2");
  unit_assert(parameters->get_group() == "Group2");
  unit_assert(parameters->get_subgroup() == "");

  // value_logical()

  unit_func("value_logical");

  parameters->set_group("Logical");
  
  unit_assert (parameters->value_logical("test_true")  == true);
  unit_assert (parameters->value_logical("test_false") == false);
  unit_assert (parameters->value_logical("test_none",true) == true);
  unit_assert (parameters->value_logical("test_none",false) == false);

  // value_integer()

  unit_func("value_integer");

  parameters->set_group("Integer");
  
  unit_assert (parameters->value_integer("test_1")  == 1);
  unit_assert (parameters->value_integer("test_37") == 37);
  unit_assert (parameters->value_integer("test_none",58) == 58);

  // value_scalar()
  
  unit_func("value_scalar");

  parameters->set_group("Scalar");
  
  unit_assert (parameters->value_scalar("test_1_5")  == 1.5);
  unit_assert (parameters->value_scalar("test_37_25") == 37.25);
  unit_assert (parameters->value_scalar("test_none",58.75) == 58.75);

  // Constant scalar expressions
  // subgroups

  unit_func("value_scalar");

  parameters->set_group("Scalar");

  parameters->set_subgroup("const_scalar_1");

  unit_assert(parameters->value_scalar("num1") == 30.625);
  unit_assert(parameters->value_scalar("num2") == 18.375);
  unit_assert(parameters->value_scalar("num3") == 150.0625);
  unit_assert(parameters->value_scalar("num4") == 4.0000000000);

  parameters->set_subgroup("const_scalar_2");

  unit_assert(parameters->value_scalar("num1") == 36.750);
  unit_assert(parameters->value_scalar("num2") == 67.375);
  unit_assert(parameters->value_scalar("num3") == -30.625);

  // Strings

  parameters->set_group("String");
  unit_assert(parameters->value_string("str1") == "testing");
  unit_assert(parameters->value_string("str1","blah") == "testing");
  unit_assert(parameters->value_string("none","blah") == "blah");

  // Variable scalar expressions

  unit_func("evaluate_scalar");

  double x[] = { 1, 2, 3};
  double y[] = {5 , 4, 3};
  double z[] = {8, 9, 10};
  double t[] = {-1, 2, -7};
  double values_scalar[] = {0,0,0};
  double deflts_scalar[] = {-1,-2,-3};

  parameters->set_group("Scalar_expr");

  parameters->set_subgroup("var_scalar_1");

  parameters->evaluate_scalar("num1",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (values_scalar[0]==x[0]);
  unit_assert (values_scalar[1]==x[1]);
  unit_assert (values_scalar[2]==x[2]);

  
  parameters->evaluate_scalar("num2",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (values_scalar[0]==x[0]+3.0);
  unit_assert (values_scalar[1]==x[1]+3.0);
  unit_assert (values_scalar[2]==x[2]+3.0);

  parameters->evaluate_scalar("num3",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (values_scalar[0]==x[0]+y[0]+z[0]+t[0]);
  unit_assert (values_scalar[1]==x[1]+y[1]+z[1]+t[1]);
  unit_assert (values_scalar[2]==x[2]+y[2]+z[2]+t[2]);

  parameters->set_subgroup("var_scalar_2");

  parameters->evaluate_scalar("num1",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (CLOSE(values_scalar[0],sin(x[0])));
  unit_assert (CLOSE(values_scalar[1],sin(x[1])));
  unit_assert (CLOSE(values_scalar[2],sin(x[2])));

  parameters->evaluate_scalar("num2",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (CLOSE(values_scalar[0],atan(y[0]/3.0+3*t[0])));
  unit_assert (CLOSE(values_scalar[1],atan(y[1]/3.0+3*t[1])));
  unit_assert (CLOSE(values_scalar[2],atan(y[2]/3.0+3*t[2])));

  // Logical expressions

  unit_func("evaluate_logical");

  bool values_logical[] = {false, false, false};
  bool deflt_logical[] = {true, false,true};

  parameters->set_group("Logical_expr");
  parameters->set_subgroup("var_logical_1");

  parameters->evaluate_logical("num1",3,values_logical,deflt_logical,x,y,z,t);
  unit_assert (values_logical[0] == (x[0] < y[0]));
  unit_assert (values_logical[1] == (x[1] < y[1]));
  unit_assert (values_logical[2] == (x[2] < y[2]));

  parameters->evaluate_logical("num2",3,values_logical,deflt_logical,x,y,z,t);
  unit_assert (values_logical[0] == (x[0] + y[0] >= t[0] + 3.0));
  unit_assert (values_logical[1] == (x[1] + y[1] >= t[1] + 3.0));
  unit_assert (values_logical[2] == (x[2] + y[2] >= t[2] + 3.0));

  parameters->evaluate_logical("num3",3,values_logical,deflt_logical,x,y,z,t);
  unit_assert (values_logical[0] == (x[0] == y[0]));
  unit_assert (values_logical[1] == (x[1] == y[1]));
  unit_assert (values_logical[2] == (x[2] == y[2]));

  // Lists

  unit_func("value_list");

  parameters->set_group("List");

  unit_assert(parameters->list_length("num1") == 5);
  unit_assert(parameters->list_value_scalar(0,"num1") == 1.0);
  unit_assert(parameters->list_value_logical(1,"num1") == true);
  unit_assert(parameters->list_value_integer(2,"num1") == 37);
  unit_assert(parameters->list_value_string(3,"num1") == "string");
  printf ("%s\n",parameters->list_value_string(3,"num1").c_str());
  parameters->list_evaluate_scalar(4,"num1",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (values_scalar[0] == (x[0]-y[0]+2.0*z[0]));
  unit_assert (values_scalar[1] == (x[1]-y[1]+2.0*z[1]));
  unit_assert (values_scalar[2] == (x[2]-y[2]+2.0*z[2]));

  delete parameters;

}
