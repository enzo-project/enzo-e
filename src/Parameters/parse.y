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
#include "type_op.h" 

  const char * type_name[]  = {
    "unknown",
    "sentinel",
    "integer",
    "scalar",
    "string",
    "identifier",
    "logical",
    "list",
    "scalar_expr",
    "logical_expr" };

  /* Structure for storing a single parameter / value pair in a linked list */
   
  struct op_node {
    enum type_parameter type;
     union {
       enum op_type oper;   /* operator */
       double scalar_value; /* number */ 
       char var_value;      /* variable x,y,z,t */
     };
    struct op_node * left;
    struct op_node * right;
  };

  struct op_node * new_oper
    (enum type_parameter type,
     enum op_type oper,
     struct op_node * left, 
     struct op_node * right)
  {
    struct op_node * node = malloc (sizeof (struct op_node));
    node->oper = oper;
    node->type = type;
    node->left = left;
    node->right = right;
  }

  struct op_node * new_scalar
    (enum type_parameter type,
     double value,
     struct op_node * left, 
     struct op_node * right)
  {
    struct op_node * node = malloc (sizeof (struct op_node));
    node->type = type;
    node->scalar_value = value;
    node->left = NULL;
    node->right = NULL;
  }
  struct op_node * new_var
    (enum type_parameter type,
     char value,
     struct op_node * left, 
     struct op_node * right)
  {
    struct op_node * node = malloc (sizeof (struct op_node));
    node->type = type;
    node->var_value = value;
    node->left = NULL;
    node->right = NULL;
  }

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
      struct op_node    * op_value;    /* expression tree */
    };
    struct param_type *   next;
  } ;

  /* The head of the linked list of parameter / value pairs */

  struct param_type * param_head = NULL; /* head of entire list */
  struct param_type * param_curr = NULL; /* head of current list */

  /* The current group, subgroup, and parameter type */

  char *              current_parameter = NULL;
  enum type_parameter current_type      = type_sentinel;

  /* Function to update parameter's groups once the group is known */

  void update_group (char * group)
    {
      struct param_type * p = param_curr;
      while (p->next->type != type_sentinel && p->next->group == NULL) {
	p->next->group = strdup(group);
        p  = p -> next;
      }
    }

  /* Function to update parameter's subgroups once the subgroup is known */

  void update_subgroup (char * subgroup)
    {
      struct param_type * p = param_curr;
      while (p->next->type != type_sentinel && p->next->subgroup == NULL) {
	p->next->subgroup = strdup(subgroup);
        p  = p -> next;
      }
    }

  void insert_param(struct param_type * head, struct param_type * new)
  {
     new->next  = head->next;
     head->next = new;
  }

  /* Function for creating and inserting a new parameter / value pair */
  /* in the linked list */

  struct param_type * new_param ()
  {
    /* Create the new node */

     struct param_type * p = 
       (struct param_type *) malloc (sizeof (struct param_type));

   /* Fill in the non-type-specific values for the new node */

     p->group     = NULL;
     p->subgroup  = NULL;
     p->parameter = strdup(current_parameter);
     current_type = type_unknown;

     insert_param(param_curr,p);

   /* Clear variables for the next assignment */

     return p;
  }

  /* New string parameter assignment */
  void new_string_param (char * value)
  {
    struct param_type * p = new_param();
    p->type         = type_string;
    p->string_value = strdup(value);

  }

  /* New empty parameter assignment: FIRST NODE IN LIST IS A SENTINEL  */
  struct param_type * new_sentinel_param ()
  {
    struct param_type * p = 
      (struct param_type *) malloc (sizeof (struct param_type));
    p->group     = NULL;
    p->subgroup  = NULL;
    p->parameter = NULL;
    p->type      = type_sentinel;
    p->next       = p;
    p->list_value = NULL;

    return p;
  }

  /* New list parameter assignment */
  void new_list_param (struct param_type * curr)
  {
    struct param_type * p = new_param();
    p->type       = type_list;
    p->list_value = curr;
  }

  /* New scalar parameter assignment */
  void new_scalar_param (double value)
  {
    struct param_type * p = new_param();
    p->type         = type_scalar;
    p->scalar_value = value;
  }

  /* New logical parameter assignment */
  void new_logical_param (int value)
  {
    struct param_type * p = new_param();
    p->type          = type_logical;
    p->logical_value = value;
  }

  /* New integer parameter assignment */
  void new_integer_param (int value)
  {
    struct param_type * p = new_param();
    p->type          = type_integer;
    p->integer_value = value;
  }

  void new_parameter()
  {
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
       break;
     case type_scalar_expr:
       new_string_param(yylval.string_type);
       param_curr->next->type = type_scalar_expr;
       break;
     case type_logical_expr:
       new_string_param(yylval.string_type);
       param_curr->next->type = type_logical_expr;
       break;
    default:
       printf ("%s:%d Unknown type %d\n",__FILE__,__LINE__,current_type);
       break;
     }
  }

  char * strcat3 (const char * s1,const char * s2,const char * s3)
  {
    char * s = malloc (strlen(s1) + strlen(s2) + strlen(s3) + 1);
    strcpy(s,s1);
    strcpy(s+strlen(s1),s2);
    strcpy(s+strlen(s1)+strlen(s2),s3);
    return s;
  }

  char * ftoa (double f)
    { 
      char * a = malloc(25); 
      sprintf (a,"%24.16e",f);
      return a;
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

%token <string_type> VARIABLE

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
 parameter '=' parameter_value { new_parameter(); }
 ;

parameter_value : 
   STRING                      { current_type = type_string;
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
 '[' { 
   struct param_type * p = new_sentinel_param();
   p->list_value = param_curr;
   new_list_param(p);
   param_curr = p;
 }
LIST_END:
 ']' { param_curr = param_curr->list_value; }


list_elements:
  parameter_value { new_parameter();  }
 | list_elements  ',' parameter_value    { new_parameter(); }

{ }
;


constant_logical_expression: 
'(' constant_logical_expression ')' { $$ = $2; }
 | constant_scalar_expression  LE constant_scalar_expression  { $$ = $1 <= $3; }
 | constant_scalar_expression  GE constant_scalar_expression  { $$ = $1 >= $3; }
 | constant_scalar_expression  '<' constant_scalar_expression { $$ = $1 <  $3; }
 | constant_scalar_expression  '>' constant_scalar_expression { $$ = $1 >  $3; }
 | constant_scalar_expression  EQ constant_scalar_expression  { $$ = $1 == $3; }
 | constant_scalar_expression  NE constant_scalar_expression  { $$ = $1 != $3; }
 | constant_logical_expression OR constant_logical_expression { $$ = $1 || $3; }
 | constant_logical_expression AND constant_logical_expression {$$ = $1 && $3; }
 | LOGICAL { printf ("%d\n",$1);$$ = $1; }
;

constant_scalar_expression: 
 '(' constant_scalar_expression ')'               { $$ = $2; }
 | constant_scalar_expression '+' constant_scalar_expression { $$ = $1 + $3;}
 | constant_scalar_expression '-' constant_scalar_expression { $$ = $1 - $3;}
 | constant_scalar_expression '*' constant_scalar_expression { $$ = $1 * $3;}
 | constant_scalar_expression '/' constant_scalar_expression { $$ = $1 / $3;}
 | SIN '(' constant_scalar_expression ')' { $$ = sin($3); }
| SCALAR { printf ("%g\n",$1); $$ = $1;}
 ;


constant_integer_expression: 
 '(' constant_integer_expression ')'               { $$ = $2; }
 | constant_integer_expression '+' constant_integer_expression { $$ = $1 + $3;}
 | constant_integer_expression '-' constant_integer_expression { $$ = $1 - $3;}
 | constant_integer_expression '*' constant_integer_expression { $$ = $1 * $3;}
 | constant_integer_expression '/' constant_integer_expression { $$ = $1 / $3;}
 | INTEGER { printf ("%d\n",$1); $$ = $1;}
 ;


variable_scalar_expression: 
'(' variable_scalar_expression ')' 
 { yylval.string_type = strcat3 ("(",$2,")");}
 | variable_scalar_expression '+' constant_scalar_expression 
 { yylval.string_type = strcat3 (strdup($1)," + ", ftoa($3)); }
 | constant_scalar_expression '+' variable_scalar_expression
 { yylval.string_type = strcat3 (ftoa($1)," + ", strdup($3)); }
 | variable_scalar_expression '+' variable_scalar_expression
 { yylval.string_type = strcat3 (strdup($1)," + ", strdup($3)); 
   printf ("(%s)\n",yylval.string_type); }
 | variable_scalar_expression '-' constant_scalar_expression
 { yylval.string_type = strcat3 (strdup($1)," - ", ftoa($3)); }
 | constant_scalar_expression '-' variable_scalar_expression
 { yylval.string_type = strcat3 (ftoa($1)," - ", strdup($3)); }
 | variable_scalar_expression '-' variable_scalar_expression
 { yylval.string_type = strcat3 (strdup($1)," - ", strdup($3)); }
 | variable_scalar_expression '*' constant_scalar_expression
 { yylval.string_type = strcat3 (strdup($1)," * ", ftoa($3)); }
 | constant_scalar_expression '*' variable_scalar_expression
 { yylval.string_type = strcat3 (ftoa($1)," * ", strdup($3)); }
 | variable_scalar_expression '*' variable_scalar_expression
 { yylval.string_type = strcat3 (strdup($1)," * ", strdup($3)); }
 | variable_scalar_expression '/' constant_scalar_expression
 { yylval.string_type = strcat3 (strdup($1)," / ", ftoa($3)); }
 | constant_scalar_expression '/' variable_scalar_expression
 { yylval.string_type = strcat3 (ftoa($1)," / ", strdup($3)); }
 | variable_scalar_expression '/' variable_scalar_expression
 { yylval.string_type = strcat3 (strdup($1)," / ", strdup($3)); }
 | SIN '(' variable_scalar_expression ')'
 { yylval.string_type = strcat3 ("sin(",$3,")"); }
 | VARIABLE { printf ("%s\n",$1); }
 ;


variable_logical_expression: 
 '(' variable_logical_expression ')' { }
 | variable_scalar_expression  LE constant_scalar_expression  
{ yylval.string_type = strcat3 (strdup($1)," <= ", ftoa($3)); }
 | constant_scalar_expression  LE variable_scalar_expression  
{ yylval.string_type = strcat3 (ftoa($1)," <= ", strdup($3)); }
 | variable_scalar_expression  LE variable_scalar_expression  
{ yylval.string_type = strcat3 (strdup($1)," <= ", strdup($3)); }
 | variable_scalar_expression  GE constant_scalar_expression 
{ yylval.string_type = strcat3 (strdup($1)," >= ", ftoa($3)); }
 | constant_scalar_expression  GE variable_scalar_expression 
{  yylval.string_type = strcat3 (ftoa($1)," >= ", strdup($3)); }
 | variable_scalar_expression  GE variable_scalar_expression 
{  yylval.string_type = strcat3 (strdup($1)," >= ", strdup($3)); }
 | variable_scalar_expression  '<' constant_scalar_expression 
{   yylval.string_type =  strdup ("7");}
 | constant_scalar_expression  '<' variable_scalar_expression 
{   yylval.string_type =  strdup ("7");}
 | variable_scalar_expression  '<' variable_scalar_expression 
{   yylval.string_type =  strdup ("7");}
 | variable_scalar_expression  '>' constant_scalar_expression 
{   yylval.string_type =  strdup ("10");}
 | constant_scalar_expression  '>' variable_scalar_expression 
{   yylval.string_type =  strdup ("10");}
 | variable_scalar_expression  '>' variable_scalar_expression 
{   yylval.string_type =  strdup ("10");}
 | variable_scalar_expression  EQ constant_scalar_expression 
{   yylval.string_type =  strdup ("13");}
 | constant_scalar_expression  EQ variable_scalar_expression 
{   yylval.string_type =  strdup ("13");}
 | variable_scalar_expression  EQ variable_scalar_expression 
{   yylval.string_type =  strdup ("13");}
 | variable_scalar_expression  NE constant_scalar_expression 
{   yylval.string_type =  strdup ("16");}
 | constant_scalar_expression  NE variable_scalar_expression 
{   yylval.string_type =  strdup ("16");}
 | variable_scalar_expression  NE variable_scalar_expression 
{   yylval.string_type =  strdup ("16");}
 | variable_logical_expression OR constant_logical_expression 
{   yylval.string_type =  strdup ("19");}
 | constant_logical_expression OR variable_logical_expression 
{   yylval.string_type =  strdup ("19");}
 | variable_logical_expression OR variable_logical_expression 
{   yylval.string_type =  strdup ("19");}
 | variable_logical_expression AND constant_logical_expression
{   yylval.string_type =  strdup ("22");}
 | constant_logical_expression AND variable_logical_expression
{   yylval.string_type =  strdup ("22");}
 | variable_logical_expression AND variable_logical_expression
{   yylval.string_type =  strdup ("22");}
;



%%

struct param_type * 
cello_parameters_read(FILE * fp)
{
  /* initialize the linked list with an initial sentinel (sentinel) node */
  param_head = param_curr = new_sentinel_param();
  
  yyrestart(fp);
  yyparse();
  return param_head;
}

void cello_parameters_print_list(struct param_type * head)
{
/*   if (head==NULL) return; */
  struct param_type * p = head->next;
  int count = 0;
  while (p && p->type != type_sentinel && count++ < 100) {
    if (p->group != NULL) {
      printf ("%s %s:%s:%s = ", 
	      type_name[p->type],p->group, p->subgroup, p->parameter);
    } else {
      /* list element */
      printf ("   %s %s = ", 
	      type_name[p->type], p->parameter);
    }
    switch (p->type) {
    case type_scalar:  printf ("%g\n",p->scalar_value);  break;
    case type_integer: printf ("%d\n",p->integer_value); break;
    case type_string:  printf ("%s\n",p->string_value);  break;
    case type_logical: printf ("%s\n",p->logical_value ? "true" : "false"); break;
    case type_list:    
      printf ("[\n"); 
      cello_parameters_print_list(p->list_value);
      printf ("]\n"); 
      break;
    case type_logical_expr:
      printf ("%s\n",p->string_value);  
      break;
    case type_scalar_expr:
      printf ("%s\n",p->string_value);  
      break;
    default: printf ("unknown type\n"); break;
    }
    p = p->next;
  }
}

void cello_parameters_print()
{
  cello_parameters_print_list(param_head);
}

