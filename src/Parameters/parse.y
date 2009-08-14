%{
#include <stdio.h>
#include "parse.tab.h"
%}

%token GROUP
%token PARAMETER

%token SCALAR
%token STRING
%token CONSTANT
%token LIST_BEGIN
%token LIST_END
%token GROUP_BEGIN
%token GROUP_END

%token TRUE
%token FALSE
%token LE
%token GE
%token NE
%token EQ
%token AND
%token OR

/* double foo (double) */

%token ACOS
%token ACOSH
%token ASIN
%token ASINH
%token ATAN
%token ATANH
%token CBRT
%token CEIL
%token COS
%token COSH
%token ERFC
%token ERF
%token EXP
%token EXPM1
%token FABS
%token FLOOR
%token GAMMA
%token J0
%token J1
%token LGAMMA
%token LOG10
%token LOG1P
%token LOGB
%token LOG
%token SIN
%token SINH
%token SQRT
%token TAN
%token TANH
%token Y0
%token Y1
%token RINT

/* int foo (double) */

/* int    ilogb(double); */
/* int    isnan(double); */

/*  */
/* double foo (double,double) */

%token ATAN2
%token FMOD
%token HYPOT
%token NEXTAFTER
%token POW
%token REMAINDER
%token SCALB

/* double jn(int, double); */
/* double ldexp(double, int); */
/* double yn(int, double); */

%%

file : /* nothing */
 | file group { printf ("file group EOF\n") }
 ;

group : group_name '{' parameter_list '}' 
 { printf ("group { ... }\n"); }
 ;

group_name : GROUP 
{ printf ("group_name\n"); }
 ;

parameter_list : /* nothing */
 | parameter_list parameter_assignment ';'
{ printf ("parameter_list\n") }
 ;

parameter_assignment : parameter_name '=' parameter_value
  { printf ("parameter_assignment\n"); }
 ;

parameter_name : PARAMETER 
   { printf ("parameter_name\n"); }
 ;

parameter_value : STRING | scalar_expression
{ printf ("parameter_value\n"); }
;

scalar_expression: scalar_factor
 | scalar_expression '+' scalar_factor { printf (" + \n"); }
 | scalar_expression '-' scalar_factor { printf (" - \n"); }
 ;

scalar_factor: scalar_term
 | scalar_factor      '*' scalar_term   { printf (" *\n"); }
 | scalar_factor      '/' scalar_term   { printf (" /\n"); }
 ;

scalar_term: SCALAR { printf ("SCALAR\n"); };

%%
main(int argc, char **argv)
{
  yyparse();
}

yyerror(char *s)
{
  fprintf(stderr, "error: %s\n", s);
}
