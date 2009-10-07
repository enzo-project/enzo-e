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
#include <math.h>
#include <string.h>
#include <malloc.h>

#include "parse.tab.h"

#define TRACE  printf ("%s:%d\n",__FILE__,__LINE__); fflush(stdout);

#include "type_parameter.h"

  /* Structure for storing a single parameter / value pair in a linked list */
   
  struct param_type {
    char * group;
    char * subgroup;
    char * parameter;
    enum type_parameter type;
    union  {
       int                 logical_value; 
       double              scalar_value; 
       char *              string_value;
       struct param_type * list_value;
    };
    struct param_type *   next;
  } ;

  /* The head of the linked list of parameter / value pairs */

  struct param_type * param_head = NULL;

  /* The current group, subgroup, and parameter type */

  char *              current_parameter = 0;
  enum type_parameter current_type      = type_unknown;

  /* Function for creating and inserting a new parameter / value pair */
  /* in the linked list */


  void update_group (char * group)
    {
      printf ("update_group (%s)\n",group);
      struct param_type * p = param_head;
      while (p && p->group == NULL) {
	p->group = strdup(group);
        p  = p -> next;
      }
    }

  void update_subgroup (char * subgroup)
    {
      printf ("update_subgroup (%s)\n",subgroup);
      struct param_type * p = param_head;
      while (p && p->subgroup == NULL) {
	p->subgroup = strdup(subgroup);
        p  = p -> next;
      }
    }

  struct param_type * new_param ()
  {
    /* Create the new node */

    TRACE;
     struct param_type * p = 
       (struct param_type *) malloc (sizeof (struct param_type));

   /* Fill in the non-type-specific values for the new node */

    TRACE;
    
     p->parameter    = strdup(current_parameter);

   /* Update the linked list pointers */

    TRACE;
     p->next         = param_head;
 
     param_head      = p;

   /* Clear variables for the next assignment */

    TRACE;
     current_type = type_unknown;

    TRACE;
     return p;
  }

  /* New string parameter assignment */
  void new_string_param (char * value)
  {
    /* Create the new node */
 
     struct param_type * p = new_param();

     p->type         = type_string;
     p->string_value = strdup(value);

  }

  /* New scalar parameter assignment */
  void new_scalar_param (double value)
  {
    struct param_type * p = new_param();

    p->type         = type_scalar;
    p->scalar_value = value;
  }

%}


%union { int integer_type;  double scalar_type;  char * string_type; }

%token <string_type> GROUP_NAME
%token <string_type> STRING
%token <scalar_type> SCALAR
%token <integer_type> LOGICAL
%token <string_type> IDENTIFIER
%type <scalar_type>  constant_scalar_expression
%type <integer_type> constant_logical_expression
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
   GROUP_NAME '{' parameter_list '}'               { update_group($1);
                                                     update_subgroup(""); }
 | GROUP_NAME '{' parameter_list ';' '}'           { update_group($1);
                                                     update_subgroup(""); }
 | GROUP_NAME IDENTIFIER '{' parameter_list '}'    { update_group($1); 
                                                     update_subgroup($2); }
 | GROUP_NAME IDENTIFIER '{' parameter_list ';' '}'{ update_group($1); 
                                                     update_subgroup($2); }

subgroup : 
    IDENTIFIER  '{' parameter_list '}'      { update_subgroup($1);
                                              printf ("subgroup\n"); }
 |  IDENTIFIER  '{' parameter_list ';' '}'  { update_subgroup($1);
                                              printf ("subgroup\n"); }

parameter_list : 
   parameter_assignment                     { }
 | subgroup                                 {  }
 | parameter_list ';' parameter_assignment  {  }
 | parameter_list ';' subgroup              {  } 
{  }
 ;

parameter :
  IDENTIFIER { current_parameter = strdup($1);} 

parameter_assignment : 
   parameter '=' parameter_value { 
     switch (current_type) {
     case type_string: 
       new_string_param(yylval.string_type);
       break;
     case type_scalar:
       new_scalar_param(yylval.scalar_type);
       break;
     case type_integer:
       break;
     case type_scalar_expr:
       break;
     case type_logical_expr:
       break;
     case type_list:
       break;
     case type_identifier:
       break;
    default:
       break;
     }
     printf ("assignment\n"); }
 ;

parameter_value : 
   STRING                      { current_type = type_string; 
                                 yylval.string_type = strdup($1);}
 | constant_scalar_expression  { current_type = type_scalar; 
                                 yylval.scalar_type = $1;}
 | constant_logical_expression { current_type = type_integer; }
 | variable_scalar_expression  { current_type = type_scalar_expr; }
 | variable_logical_expression { current_type = type_logical_expr; }
 | list                        { current_type = type_list; }
 | IDENTIFIER                  { current_type = type_identifier; }
{  }
;

list: '[' list_elements ']'


list_elements:
   parameter_value { printf ("begin list\n"); }
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
 | variable_scalar_expression  LE constant_scalar_expression  { }
 | constant_scalar_expression  LE variable_scalar_expression { }
 | variable_scalar_expression  LE variable_scalar_expression { }
 | variable_scalar_expression  GE constant_scalar_expression { }
 | constant_scalar_expression  GE variable_scalar_expression { }
 | variable_scalar_expression  GE variable_scalar_expression { }
 | variable_scalar_expression  '<' constant_scalar_expression { }
 | constant_scalar_expression  '<' variable_scalar_expression { }
 | variable_scalar_expression  '<' variable_scalar_expression { }
 | variable_scalar_expression  '>' constant_scalar_expression { }
 | constant_scalar_expression  '>' variable_scalar_expression { }
 | variable_scalar_expression  '>' variable_scalar_expression { }
 | variable_scalar_expression  EQ constant_scalar_expression { }
 | constant_scalar_expression  EQ variable_scalar_expression { }
 | variable_scalar_expression  EQ variable_scalar_expression { }
 | variable_scalar_expression  NE constant_scalar_expression { }
 | constant_scalar_expression  NE variable_scalar_expression { }
 | variable_scalar_expression  NE variable_scalar_expression { }
 | variable_logical_expression OR constant_logical_expression {  }
 | constant_logical_expression OR variable_logical_expression {  }
 | variable_logical_expression OR variable_logical_expression {  }
 | variable_logical_expression AND constant_logical_expression {  }
 | constant_logical_expression AND variable_logical_expression {  }
 | variable_logical_expression AND variable_logical_expression {  }
;


constant_scalar_expression: 
 '(' constant_scalar_expression ')'               { $$ = $2; }
 | constant_scalar_expression '+' constant_scalar_expression { $$ = $1 + $3;}
 | constant_scalar_expression '-' constant_scalar_expression { $$ = $1 - $3;}
 | constant_scalar_expression '*' constant_scalar_expression { $$ = $1 * $3;}
 | constant_scalar_expression '/' constant_scalar_expression { $$ = $1 / $3;}
 | SIN '(' constant_scalar_expression ')' { $$ = sin($3); }
 | SCALAR { printf ("S = %g\n",$1); $$ = $1;}
 ;

variable_scalar_expression: 
   '(' variable_scalar_expression ')'      { }
 | variable_scalar_expression '+' constant_scalar_expression { }
 | constant_scalar_expression '+' variable_scalar_expression { }
 | variable_scalar_expression '+' variable_scalar_expression { }
 | variable_scalar_expression '-' constant_scalar_expression { }
 | constant_scalar_expression '-' variable_scalar_expression { }
 | variable_scalar_expression '-' variable_scalar_expression { }
 | variable_scalar_expression '*' constant_scalar_expression { }
 | constant_scalar_expression '*' variable_scalar_expression { }
 | variable_scalar_expression '*' variable_scalar_expression { }
 | variable_scalar_expression '/' constant_scalar_expression { }
 | constant_scalar_expression '/' variable_scalar_expression { }
 | variable_scalar_expression '/' variable_scalar_expression { }
 | SIN '(' variable_scalar_expression ')' { }
 | VARIABLE { }
 ;


%%

void cello_parameters_read(FILE * fp)
{
  yyrestart(fp);
  yyparse();
}

void cello_parameters_print()
{
  struct param_type * p = param_head;
  while (p != NULL) {
    printf ("group = %s subgroup = %s  parameter = %s value = ", 
	    p->group, p->subgroup, p->parameter);
    switch (p->type) {
    case type_unknown: printf ("???\n"); break;
    case type_scalar: printf ("%g\n",p->scalar_value); break;
    case type_string: printf ("%s\n",p->string_value); break;
    case type_logical: printf ("%s\n",p->logical_value ? "true" : "false"); break;
    }
    p = p->next;
  }
}

