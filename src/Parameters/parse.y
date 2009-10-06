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
%}

%token GROUP_NAME
%token IDENTIFIER

%union {int logical_type; double scalar_type; char * string_type;}

%token <string_type> STRING
%token <scalar_type> SCALAR
%token <logical_type> LOGICAL
%type <scalar_type>  constant_scalar_expression
%type <logical_type> constant_logical_expression
%type <string_type>  variable_scalar_expression
%type <string_type>  variable_logical_expression

%token VARIABLE
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

%token SIN
 /* %token ACOS ACOSH ASIN ASINH ATAN ATANH CBRT CEIL COS COSH ERFC ERF */
 /* %token EXP EXPM1 FABS FLOOR GAMMA J0 J1 LGAMMA LOG10 LOG1P LOGB */
 /* %token LOG SIN SINH SQRT TAN TANH Y0 Y1 RINT */

/* double foo (double,double) */

 /* %token ATAN2 FMOD HYPOT NEXTAFTER POW REMAINDER SCALB */

%%

file : /* nothing */
 | file group { }
 ;

group: 
    GROUP_NAME '{' parameter_list '}'                 { printf ("group \n"); }
 |  GROUP_NAME '{' parameter_list ';' '}'             { printf ("group;\n"); }
 |  GROUP_NAME IDENTIFIER '{' parameter_list '}'      { printf ("group subgroup\n"); }
 |  GROUP_NAME IDENTIFIER '{' parameter_list ';' '}'  { printf ("group subgroup;\n"); }

;
subgroup : 
    IDENTIFIER  '{' parameter_list '}'      { printf ("subgroup \n"); }
 |  IDENTIFIER  '{' parameter_list ';' '}'  { printf ("subgroup; \n"); }

parameter_list : 
   parameter_assignment                     {  printf ("parameter assignment\n"); }
 | subgroup                                 {  }
 | parameter_list ';' parameter_assignment  {  }
 | parameter_list ';' subgroup              {  } 
{  }
 ;

parameter_assignment : IDENTIFIER '=' parameter_value
  {  }
 ;

parameter_value : 
   STRING                     { printf ("string %s\n",$1); free $1;}
 | constant_scalar_expression  { printf ("constant scalar expression\n"); }
 | constant_logical_expression { printf ("constant logical expression\n"); }
 | variable_scalar_expression  { printf ("variable scalar expression\n"); }
 | variable_logical_expression { printf ("variable logical expression\n"); }
 | list                { printf ("list\n"); }
 | IDENTIFIER          { printf ("identifier\n"); }
{  }
;

list: '[' list_elements ']'


list_elements:
  parameter_value 
  | list_elements  ',' parameter_value
{ }
;


constant_logical_expression: 
 '(' constant_logical_expression ')' { }
 | constant_scalar_expression  LE constant_scalar_expression { }
 | constant_scalar_expression  GE constant_scalar_expression { }
 | constant_scalar_expression  '<' constant_scalar_expression { }
 | constant_scalar_expression  '>' constant_scalar_expression { }
 | constant_scalar_expression  EQ constant_scalar_expression { }
 | constant_scalar_expression  NE constant_scalar_expression { }
 | constant_logical_expression OR constant_logical_expression {  }
 | constant_logical_expression AND constant_logical_expression {  }
 | LOGICAL { }
;

variable_logical_expression: 
 '(' variable_logical_expression ')' { }
 | variable_scalar_expression  LE constant_scalar_expression  { printf (" VSE <= CSE\n"); }
 | constant_scalar_expression  LE variable_scalar_expression { printf (" CSE <= VSE\n"); }
 | variable_scalar_expression  LE variable_scalar_expression { printf (" VSE <= VSE\n"); }
 | variable_scalar_expression  GE constant_scalar_expression { printf (" VSE >= CSE\n"); }
 | constant_scalar_expression  GE variable_scalar_expression { printf (" CSE >= VSE\n"); }
 | variable_scalar_expression  GE variable_scalar_expression { printf (" VSE >= VSE\n"); }
 | variable_scalar_expression  '<' constant_scalar_expression { printf (" VSE < CSE\n"); }
 | constant_scalar_expression  '<' variable_scalar_expression { printf (" CSE < VSE\n"); }
 | variable_scalar_expression  '<' variable_scalar_expression { printf (" VSE < VSE\n"); }
 | variable_scalar_expression  '>' constant_scalar_expression { printf (" VSE > CSE\n"); }
 | constant_scalar_expression  '>' variable_scalar_expression { printf (" CSE > VSE\n"); }
 | variable_scalar_expression  '>' variable_scalar_expression { printf (" VSE > VSE\n"); }
 | variable_scalar_expression  EQ constant_scalar_expression { printf (" VSE == CSE\n"); }
 | constant_scalar_expression  EQ variable_scalar_expression { printf (" CSE == VSE\n"); }
 | variable_scalar_expression  EQ variable_scalar_expression { printf (" VSE == VSE\n"); }
 | variable_scalar_expression  NE constant_scalar_expression { printf (" VSE != CSE\n"); }
 | constant_scalar_expression  NE variable_scalar_expression { printf (" CSE != VSE\n"); }
 | variable_scalar_expression  NE variable_scalar_expression { printf (" VSE != VSE\n"); }
 | variable_logical_expression OR constant_logical_expression {  printf (" VSE || CSE\n"); }
 | constant_logical_expression OR variable_logical_expression {  printf (" CSE || VSE\n"); }
 | variable_logical_expression OR variable_logical_expression {  printf (" VSE || VSE\n"); }
 | variable_logical_expression AND constant_logical_expression {  printf (" VSE && CSE\n"); }
 | constant_logical_expression AND variable_logical_expression {  printf (" CSE && VSE\n"); }
 | variable_logical_expression AND variable_logical_expression {  printf (" VSE && VSE\n"); }
;


constant_scalar_expression: 
   '(' constant_scalar_expression ')'               { printf ("(CSE)\n"); }
 | constant_scalar_expression '+' constant_scalar_expression { printf ("CSE + CSE\n");  }
 | constant_scalar_expression '-' constant_scalar_expression { printf ("CSE - CSE\n");  }
 | constant_scalar_expression '*' constant_scalar_expression { printf ("CSE * CSE\n");  }
 | constant_scalar_expression '/' constant_scalar_expression { printf ("CSE / CSE\n");  }
 | SIN '(' constant_scalar_expression ')' { printf ("sin(CSE)\n"); }
 | SCALAR { printf ("S\n"); }
 ;

variable_scalar_expression: 
   '(' variable_scalar_expression ')'      { printf ("(VSE)"); }
 | variable_scalar_expression '+' constant_scalar_expression { printf (" VSE + CSE \n");  }
 | constant_scalar_expression '+' variable_scalar_expression { printf (" CSE + VSE \n");  }
 | variable_scalar_expression '+' variable_scalar_expression { printf (" CSE + VSE \n");  }
 | variable_scalar_expression '-' constant_scalar_expression { printf (" VSE - CSE \n");  }
 | constant_scalar_expression '-' variable_scalar_expression { printf (" CSE - VSE \n");  }
 | variable_scalar_expression '-' variable_scalar_expression { printf (" CSE - VSE \n");  }
 | variable_scalar_expression '*' constant_scalar_expression { printf (" VSE * CSE \n");  }
 | constant_scalar_expression '*' variable_scalar_expression { printf (" CSE * VSE \n");  }
 | variable_scalar_expression '*' variable_scalar_expression { printf (" CSE * VSE \n");  }
 | variable_scalar_expression '/' constant_scalar_expression { printf (" VSE / CSE \n");  }
 | constant_scalar_expression '/' variable_scalar_expression { printf (" CSE / VSE \n");  }
 | variable_scalar_expression '/' variable_scalar_expression { printf (" CSE / VSE \n");  }
 | SIN '(' variable_scalar_expression ')' { printf ("sin(VSE)\n"); }
 | VARIABLE { printf ("V\n")}
 ;


%%
main(int argc, char **argv)
{
  /*  yydebug=1; */
  yyparse();
}

