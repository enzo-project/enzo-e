%{
/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

#include <stdio.h>
#include "parse.tab.h"
  /*  int yydebug;     */
  /* #define YYDEBUG 1 */
%}

%token GROUP_NAME
%token IDENTIFIER

%union {int logical_type; double scalar_type; char * string_type;}

%token <string_type> STRING
%token <scalar_type> SCALAR
%token <logical_type> LOGICAL
%type <scalar_type>  scalar_expression
%type <logical_type> logical_expression

%token CONSTANT
%token LIST_BEGIN
%token LIST_END
%token GROUP_BEGIN
%token GROUP_END

%token LE
%token GE
%token NE
%token EQ
%token AND
%token OR

%left OR
%left AND
%left EQ NE 
%left  LE GE '<' '>'
%left '+' '-'
%left '*' '/'

/* double foo (double) */

%token ACOS ACOSH ASIN ASINH ATAN ATANH CBRT CEIL COS COSH ERFC ERF
%token EXP EXPM1 FABS FLOOR GAMMA J0 J1 LGAMMA LOG10 LOG1P LOGB
%token LOG SIN SINH SQRT TAN TANH Y0 Y1 RINT

/* double foo (double,double) */

%token ATAN2 FMOD HYPOT NEXTAFTER POW REMAINDER SCALB

%%

file : /* nothing */
 | file named_group { }
 ;

named_group: GROUP_NAME group           { printf ("Group\n"); }
          |  GROUP_NAME IDENTIFIER group { printf ("Group member\n"); }
;

group :  '{' parameter_list '}'  { }
 ;

parameter_list : 
   parameter_assignment 
 | parameter_assignment  ';' parameter_list
{  }
 ;

parameter_assignment : IDENTIFIER '=' parameter_value
  {  }
 ;

parameter_value : 
   STRING              { printf ("string %s\n",$1); free $1;}
 | scalar_expression   { printf ("scalar expression\n"); }
 | logical_expression  { printf ("logical expression\n"); }
 | list                { printf ("list\n"); }
 | group               { printf ("group\n"); }
 | IDENTIFIER          { printf ("variable\n"); }
{  }
;

list: '[' list_elements ']'


list_elements:
  parameter_value 
  | list_elements  ',' parameter_value
{ }
;


logical_expression: 
 '(' logical_expression ')' { }
 | scalar_expression  LE scalar_expression { }
 | scalar_expression  GE scalar_expression { }
 | scalar_expression  '<' scalar_expression { }
 | scalar_expression  '>' scalar_expression { }
 | scalar_expression  EQ scalar_expression { }
 | scalar_expression  NE scalar_expression { }
 | logical_expression OR logical_expression {  }
 | logical_expression AND logical_expression {  }
 | LOGICAL { }
;


scalar_expression: 
  '(' scalar_expression ')' { }
 | scalar_expression '+' scalar_expression {  }
 | scalar_expression '-' scalar_expression {  }
 | scalar_expression '*' scalar_expression {  }
 | scalar_expression '/' scalar_expression {  }
 | SCALAR { }
 | CONSTANT { }
 ;


%%
main(int argc, char **argv)
{
  /*  yydebug=1; */
  yyparse();
}

