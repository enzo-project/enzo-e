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
#include "type_parameter.h"
const char * type_name[]  = {
  "empty",
  "integer",
  "scalar",
  "string",
  "identifier",
  "logical",
  "list",
  "scalar_expr",
  "logical_expr" };



  /* Structure for storing a single parameter / value pair in a linked list */
   
  struct param_type {
    char * group;
    char * subgroup;
    char * parameter;
    enum type_parameter type;
    union  {
       int                 logical_value; 
       int                 integer_value; 
       double              scalar_value; 
       char *              string_value;
       struct param_type * list_value;
    };
    struct param_type *   next;
  } ;

  /* The head of the linked list of parameter / value pairs */

  struct param_type * param_head = NULL;
  struct param_type * param_curr = NULL;

  /* The current group, subgroup, and parameter type */

  char *              current_parameter = NULL;
  enum type_parameter current_type      = type_empty;

  /* Function to update parameter's groups once the group is known */

  void update_group (char * group)
    {
      struct param_type * p = param_curr;
      printf ("update_group %p\n",p);
      while (p && p->group == NULL) {
	p->group = strdup(group);
        p  = p -> next;
      }
    }

  /* Function to update parameter's subgroups once the subgroup is known */

  void update_subgroup (char * subgroup)
    {
      struct param_type * p = param_curr;
      printf ("update_subgroup %p\n",p);
      while (p && p->subgroup == NULL) {
	p->subgroup = strdup(subgroup);
        p  = p -> next;
      }
    }

  /* Function for creating and inserting a new parameter / value pair */
  /* in the linked list */

  struct param_type * new_param ()
  {
    /* Create the new node */

     struct param_type * p = 
       (struct param_type *) malloc (sizeof (struct param_type));

   /* Fill in the non-type-specific values for the new node */

     p->parameter = strdup(current_parameter);

   /* Update the linked list pointers */

     printf ("%p->next = %p\n",p,param_curr);
     p->next = param_curr->next;
     param_curr->next = p;
     param_head = param_curr;

   /* Clear variables for the next assignment */

     current_type = type_empty;

     return p;
  }

  /* New string parameter assignment */
  void new_string_param (char * value)
  {
    printf ("new_string_param(%s)\n",value);
    struct param_type * p = new_param();
    p->type         = type_string;
    p->string_value = strdup(value);

  }

  /* New empty parameter assignment: FIRST NODE IN LIST IS A SENTINEL  */
  void new_empty_param ()
  {
    if (param_head != NULL || param_curr != NULL) {
      printf ("Error: %s:%d Calling new_empty_param with non-empty list\n",
	      __FILE__,__LINE__);
    } else {
      struct param_type * p = 
	(struct param_type *) malloc (sizeof (struct param_type));
      p->group     = NULL;
      p->subgroup  = NULL;
      p->parameter = NULL;
      p->type = type_empty;
      p->next = NULL;

      param_head = p;
      param_curr = p;
    }
  }

  /* New list parameter assignment */
  void new_list_param (struct param_type * curr)
  {
    printf ("new_list_param()\n");
    struct param_type * p = new_param();
    p->type       = type_list;
    p->list_value = curr;
  }

  /* New scalar parameter assignment */
  void new_scalar_param (double value)
  {
    printf ("new_scalar_param(%g)\n",value);
    struct param_type * p = new_param();
    p->type         = type_scalar;
    p->scalar_value = value;
  }

  /* New logical parameter assignment */
  void new_logical_param (int value)
  {
    printf ("new_logical_param(%d)\n",value);
    struct param_type * p = new_param();
    p->type          = type_logical;
    p->logical_value = value;
  }

  /* New integer parameter assignment */
  void new_integer_param (int value)
  {
    printf ("new_integer_param(%d)\n",value);
    struct param_type * p = new_param();
    p->type          = type_integer;
    p->integer_value = value;
  }

%}


%union { int logical_type;  int integer_type; double scalar_type;  char * string_type; }

%token <string_type> GROUP_NAME
%token <string_type> STRING
%token <scalar_type> SCALAR
%token <integer_type> INTEGER
%token <logical_type> LOGICAL
%token <string_type> IDENTIFIER
%type <integer_type>  constant_integer_expression
%type <scalar_type>  constant_scalar_expression
%type <logical_type> constant_logical_expression
%type <string_type>  variable_scalar_expression
%type <string_type>  variable_logical_expression

%token VARIABLE

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
    IDENTIFIER  '{' parameter_list '}'      { update_subgroup($1); }
 |  IDENTIFIER  '{' parameter_list ';' '}'  { update_subgroup($1); }

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
     printf ("current_type = %s\n",type_name[current_type]);
     switch (current_type) {
     case type_integer:
       new_integer_param(yylval.integer_type);
       break;
     case type_scalar:
       new_scalar_param(yylval.scalar_type);
       break;
     case type_string: 
       new_string_param(yylval.string_type);
       break;
     case type_identifier:
       printf ("IDENTIFIER\n");
       break;
     case type_logical:
       new_logical_param(yylval.logical_type);
       break;
     case type_list:
       printf ("LIST\n");
       break;
     case type_scalar_expr:
       printf ("SCALAR EXPRESSION\n");
       break;
     case type_logical_expr:
       printf ("LOGICAL EXPRESSION\n");
       break;
    default:
       printf ("%s:%d Unknown type %d\n",__FILE__,__LINE__,current_type);
       break;
     }
   }
 ;

parameter_value : 
   STRING                      { current_type = type_string; 
                                 printf ("S=%s\n",$1);
                                 yylval.string_type = strdup($1);}
 | constant_integer_expression  { current_type = type_integer; 
                                 yylval.integer_type = $1;}
 | constant_scalar_expression  { current_type = type_scalar; 
                                 yylval.scalar_type = $1;}
 | constant_logical_expression { current_type = type_logical;
                                 yylval.logical_type = $1; }
 | variable_scalar_expression  { current_type = type_scalar_expr; }
 | variable_logical_expression { current_type = type_logical_expr; }
 | list                        { current_type = type_list; }
 | IDENTIFIER                  { current_type = type_identifier; }
{  }
;

list: LIST_BEGIN list_elements LIST_END {  }

LIST_BEGIN:
 '[' { printf ("[\n"); }
LIST_END:
 ']' { printf ("]\n"); }


list_elements:
   parameter_value { }
  | list_elements  ',' parameter_value
{ }
;


constant_logical_expression: 
'(' constant_logical_expression ')' { $$ = $2; }
 | constant_scalar_expression  LE constant_scalar_expression { $$ = $1 <= $3; }
 | constant_scalar_expression  GE constant_scalar_expression { $$ = $1 >= $3; }
 | constant_scalar_expression  '<' constant_scalar_expression { $$ = $1 < $3; }
 | constant_scalar_expression  '>' constant_scalar_expression { $$ = $1 > $3; }
 | constant_scalar_expression  EQ constant_scalar_expression { $$ = $1 == $3; }
 | constant_scalar_expression  NE constant_scalar_expression { $$ = $1 != $3; }
 | constant_logical_expression OR constant_logical_expression {  $$ = $1 || $3; }
 | constant_logical_expression AND constant_logical_expression {  $$ = $1 && $3; }
 | LOGICAL { printf ("L=%d\n",$1); $$ = $1; }
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
 | SCALAR { printf ("S=%g\n",$1); $$ = $1;}
 ;

constant_integer_expression: 
 '(' constant_integer_expression ')'               { $$ = $2; }
 | constant_integer_expression '+' constant_integer_expression { $$ = $1 + $3;}
 | constant_integer_expression '-' constant_integer_expression { $$ = $1 - $3;}
 | constant_integer_expression '*' constant_integer_expression { $$ = $1 * $3;}
 | constant_integer_expression '/' constant_integer_expression { $$ = $1 / $3;}
 | INTEGER { printf ("I=%d\n",$1); $$ = $1;}
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

struct param_type * 
cello_parameters_read(FILE * fp)
{
  /* initialize the linked list with an initial empty (sentinel) node */
  new_empty_param();
  
  yyrestart(fp);
  yyparse();
  return param_head;
}

void cello_parameters_print()
{
  printf ("cello_parameters_print()\n");
  struct param_type * p = param_head;
  printf ("param_head = %p\n",param_head);
  if (param_head != NULL) printf ("param_head->next = %p\n",param_head->next);
  while (p != NULL) {
    printf ("%s %s:%s:%s value = ", 
	    type_name[p->type],p->group, p->subgroup, p->parameter);
    switch (p->type) {
    case type_empty:   printf ("empty\n"); break;
    case type_scalar:  printf ("%g\n",p->scalar_value);  break;
    case type_integer: printf ("%d\n",p->integer_value); break;
    case type_string:  printf ("%s\n",p->string_value);  break;
    case type_logical: printf ("%s\n",p->logical_value ? "true" : "false"); break;
    default: printf ("unknown type\n"); break;
    }
    p = p->next;
  }
}

