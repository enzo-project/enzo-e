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
#include <assert.h>
#include <string.h>
#include <malloc.h>

#include "parse.tab.h"

#include "parse.h"

const char * node_name[] = {
  "node_unknown",
  "node_operation",
  "node_scalar",
  "node_integer",
  "node_variable",
  "node_function"
  };

const char * op_name[] = {
    "+",
    "-",
    "*",
    "/",
    "<=",
    "<",
    ">=",
    ">",
    "==",
    "!=",
    "&&",
    "||"};

  const char * parameter_name[]  = {
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


  struct node_expr * new_node_operation
    (struct node_expr * left, 
     enum enum_op oper,
     struct node_expr * right)
  {
    
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type     = enum_node_operation;
    node->op_value = oper;
    node->left     = left;
    node->right    = right;
    return node;
  }

  struct node_expr * new_node_scalar (double value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type = enum_node_scalar;
    node->scalar_value = value;
    node->left = NULL;
    node->right = NULL;
    return node;
  }
  struct node_expr * new_node_logical (int value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type = enum_node_integer;
    node->integer_value = value;
    node->left = NULL;
    node->right = NULL;
    return node;
  }
  struct node_expr * new_node_variable (char value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type = enum_node_variable;
    node->var_value = value;
    node->left = NULL;
    node->right = NULL;
    return node;
  }
  struct node_expr * new_node_function
    (double (*function)(double),
     char * function_name,
     struct node_expr * argument)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type = enum_node_function;
    node->fun_value = function;
    node->left = argument;
    node->right = NULL;
    node->function_name = function_name;
    return node;
  }


  /* The head of the linked list of parameter / value pairs */

  struct param_type * param_head = NULL; /* head of entire list */
  struct param_type * param_curr = NULL; /* head of current list */

  /* The current group, subgroup, and parameter type */

  char *              current_parameter = NULL;
  enum enum_parameter current_type      = enum_parameter_sentinel;

  /* Function to update parameter's groups once the group is known */

  void update_group (char * group)
    {
      struct param_type * p = param_curr;
      while (p->next->type != enum_parameter_sentinel && p->next->group == NULL) {
	p->next->group = strdup(group);
        p  = p -> next;
      }
    }

  /* Function to update parameter's subgroups once the subgroup is known */

  void update_subgroup (char * subgroup)
    {
      struct param_type * p = param_curr;
      while (p->next->type != enum_parameter_sentinel && p->next->subgroup == NULL) {
	p->next->subgroup = strdup(subgroup);
        p  = p -> next;
      }
    }

  void insert_param(struct param_type * head, struct param_type * new)
  {
     new->next  = head->next;
     head->next = new;
  }

  struct param_type * reverse_param(struct param_type * old_head)
  {
    /* Keep sentinel the same */

    struct param_type * new_head = old_head;

    struct param_type * p = old_head;
    struct param_type * c = p->next;
    struct param_type * n = c->next;

    do {
      /* If parameter is a list, recursively reverse it as well */
      if (c->type == enum_parameter_list) {
	c->list_value = reverse_param(c->list_value);
      }
      c->next = p;

      p=c;
      c=n;
      n=n->next;
    } while (p->type != enum_parameter_sentinel) ;

    new_head = p;
    return new_head;
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
     current_type = enum_parameter_unknown;

     insert_param(param_curr,p);

   /* Clear variables for the next assignment */

     return p;
  }

  /* New string parameter assignment */
  void new_param_string (char * value)
  {
    struct param_type * p = new_param();
    p->type         = enum_parameter_string;
    p->string_value = strdup(value);

  }

  /* New empty parameter assignment: FIRST NODE IN LIST IS A SENTINEL  */
  struct param_type * new_param_sentinel ()
  {
    struct param_type * p = 
      (struct param_type *) malloc (sizeof (struct param_type));

    p->group     = NULL;
    p->subgroup  = NULL;
    p->parameter = NULL;
    p->type      = enum_parameter_sentinel;
    p->next       = p;
    p->list_value = NULL;

    return p;
  }

  /* New list parameter assignment */
  void new_param_list (struct param_type * curr)
  {
    struct param_type * p = new_param();
    p->type       = enum_parameter_list;
    p->list_value = curr;
  }

  /* New scalar parameter assignment */
  void new_param_scalar (double value)
  {
    struct param_type * p = new_param();
    p->type         = enum_parameter_scalar;
    p->scalar_value = value;
  }

  /* New logical parameter assignment */
  void new_param_logical (int value)
  {
    struct param_type * p = new_param();
    p->type          = enum_parameter_logical;
    p->logical_value = value;
  }

  /* New integer parameter assignment */
  void new_param_integer (int value)
  {
    struct param_type * p = new_param();
    p->type          = enum_parameter_integer;
    p->integer_value = value;
  }

  /* New string parameter assignment */
  void new_param_expr (enum enum_parameter type,
		       struct node_expr * value)
  {
    struct param_type * p = new_param();
    p->type     = type;
    p->op_value = value;

  }

  void new_parameter()
  {
     switch (current_type) {
     case enum_parameter_integer:
       new_param_integer(yylval.integer_type);
       break;
     case enum_parameter_scalar:
       new_param_scalar(yylval.scalar_type);
       break;
     case enum_parameter_string: 
       new_param_string(yylval.string_type);
       break;
/*      case enum_parameter_identifier: */
/*        printf ("IDENTIFIER\n"); */
/*        break; */
     case enum_parameter_logical:
       new_param_logical(yylval.logical_type);
       break;
     case enum_parameter_list:
       break;
     case enum_parameter_scalar_expr:
       new_param_expr(enum_parameter_scalar_expr,yylval.node_type);
       break;
     case enum_parameter_logical_expr:
       new_param_expr(enum_parameter_logical_expr,yylval.node_type);
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


%union { 
  int logical_type;  
  int integer_type; 
  double scalar_type;  
  char * string_type; 
  struct node_expr * node_type;
  }

%token <string_type>  GROUP_NAME
%token <string_type>  STRING
%token <scalar_type>  SCALAR
%token <integer_type> INTEGER
%token <logical_type> LOGICAL
%token <string_type>  IDENTIFIER
%token <string_type> VARIABLE

%type <integer_type> cie
%type <scalar_type>  cse
%type <logical_type> cle

%type <node_type>  vse
%type <node_type>  vle

%token LE
%token GE
%token NE
%token EQ
%token AND
%token OR

%left OR
%left AND
%left EQ NE 
%left LE GE '<' '>'
%left '+' '-'
%left '*' '/'

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
/* %token GAMMA */
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
 STRING { current_type = enum_parameter_string;       yylval.string_type = strdup($1);}
 | cie  { current_type = enum_parameter_integer;      yylval.integer_type = $1;}
 | cse  { current_type = enum_parameter_scalar;       yylval.scalar_type = $1;}
 | cle  { current_type = enum_parameter_logical;      yylval.logical_type = $1; }
 | vse  { current_type = enum_parameter_scalar_expr;  yylval.node_type = $1; }
 | vle  { current_type = enum_parameter_logical_expr; yylval.node_type = $1; }
 | list { current_type = enum_parameter_list; }
/*  | IDENTIFIER  { current_type = enum_parameter_identifier; } */
/* {  } */
;

list: LIST_BEGIN list_elements LIST_END {  }

LIST_BEGIN:
 '[' { 
   struct param_type * p = new_param_sentinel();
   p->list_value = param_curr;
   new_param_list(p);
   param_curr = p;
 }
LIST_END:
 ']' { param_curr = param_curr->list_value; }


list_elements:
  parameter_value { new_parameter();  }
 | list_elements  ',' parameter_value    { new_parameter(); }

{ }
;


cle: 
'(' cle ')' { $$ = $2; }
 | cse LE  cse { $$ = $1 <= $3; }
 | cse GE  cse { $$ = $1 >= $3; }
 | cse '<' cse { $$ = $1 <  $3; }
 | cse '>' cse { $$ = $1 >  $3; }
 | cse EQ  cse { $$ = $1 == $3; }
 | cse NE  cse { $$ = $1 != $3; }
 | cle OR  cle { $$ = $1 || $3; }
 | cle AND cle { $$ = $1 && $3; }
 | LOGICAL { $$ = $1; }
;

cse: 
 '(' cse ')' { $$ = $2; }
 | cse '+' cse { $$ = $1 + $3;}
 | cse '-' cse { $$ = $1 - $3;}
 | cse '*' cse { $$ = $1 * $3;}
 | cse '/' cse { $$ = $1 / $3;}
| ACOS '(' cse ')' { $$ = acos($3); }
| ACOSH '(' cse ')' { $$ = acosh($3); }
| ASIN '(' cse ')' { $$ = asin($3); }
| ASINH '(' cse ')' { $$ = asinh($3); }
| ATAN '(' cse ')' { $$ = atan($3); }
| ATANH '(' cse ')' { $$ = atanh($3); }
| CBRT '(' cse ')' { $$ = cbrt($3); }
| CEIL '(' cse ')' { $$ = ceil($3); }
| COS '(' cse ')' { $$ = cos($3); }
| COSH '(' cse ')' { $$ = cosh($3); }
| ERFC '(' cse ')' { $$ = erfc($3); }
| ERF '(' cse ')' { $$ = erf($3); }
| EXP '(' cse ')' { $$ = exp($3); }
| EXPM1 '(' cse ')' { $$ = expm1($3); }
| FABS '(' cse ')' { $$ = fabs($3); }
| FLOOR '(' cse ')' { $$ = floor($3); }
/* | GAMMA '(' cse ')' { $$ = gamma($3); } */
| J0 '(' cse ')' { $$ = j0($3); }
| J1 '(' cse ')' { $$ = j1($3); }
| LGAMMA '(' cse ')' { $$ = lgamma($3); }
| LOG10 '(' cse ')' { $$ = log10($3); }
| LOG1P '(' cse ')' { $$ = log1p($3); }
| LOGB '(' cse ')' { $$ = logb($3); }
| LOG '(' cse ')' { $$ = log($3); }
| SIN '(' cse ')' { $$ = sin($3); }
| SINH '(' cse ')' { $$ = sinh($3); }
| SQRT '(' cse ')' { $$ = sqrt($3); }
| TAN '(' cse ')' { $$ = tan($3); }
| TANH '(' cse ')' { $$ = tanh($3); }
| Y0 '(' cse ')' { $$ = y0($3); }
| Y1 '(' cse ')' { $$ = y1($3); }
| RINT '(' cse ')' { $$ = rint($3); }
| SCALAR { $$ = $1;}
 ;

cie: 
 '(' cie ')' { $$ = $2; }
 | cie '+' cie { $$ = $1 + $3;}
 | cie '-' cie { $$ = $1 - $3;}
 | cie '*' cie { $$ = $1 * $3;}
 | cie '/' cie { $$ = $1 / $3;}
 | INTEGER { $$ = $1;}
 ;

vse: 
'(' vse ')'    { $2; }
 | vse '+' cse { $$ = new_node_operation ($1, enum_op_add,new_node_scalar($3)); }
 | cse '+' vse { $$ = new_node_operation (new_node_scalar($1), enum_op_add,$3); }
 | vse '+' vse { $$ = new_node_operation ($1, enum_op_add,$3); }
 | vse '-' cse { $$ = new_node_operation ($1, enum_op_sub,new_node_scalar($3)); }
 | cse '-' vse { $$ = new_node_operation (new_node_scalar($1), enum_op_sub,$3); }
 | vse '-' vse { $$ = new_node_operation ($1, enum_op_sub,$3); }
 | vse '*' cse { $$ = new_node_operation ($1, enum_op_mul,new_node_scalar($3)); }
 | cse '*' vse { $$ = new_node_operation (new_node_scalar($1), enum_op_mul,$3); }
 | vse '*' vse { $$ = new_node_operation ($1, enum_op_mul,$3); }
 | vse '/' cse { $$ = new_node_operation ($1, enum_op_div,new_node_scalar($3)); }
 | cse '/' vse { $$ = new_node_operation (new_node_scalar($1), enum_op_div,$3); }
 | vse '/' vse { $$ = new_node_operation ($1, enum_op_div,$3); }
 | ACOS   '(' vse ')' { $$ = new_node_function ( acos, "acos", $3); }
 | ACOSH  '(' vse ')' { $$ = new_node_function ( acosh, "acosh", $3); }
 | ASIN   '(' vse ')' { $$ = new_node_function ( asin, "asin", $3); }
 | ASINH  '(' vse ')' { $$ = new_node_function ( asinh, "asinh", $3); }
 | ATAN   '(' vse ')' { $$ = new_node_function ( atan, "atan", $3); }
 | ATANH  '(' vse ')' { $$ = new_node_function ( atanh, "atanh", $3); }
 | CBRT   '(' vse ')' { $$ = new_node_function ( cbrt, "cbrt", $3); }
 | CEIL   '(' vse ')' { $$ = new_node_function ( ceil, "ceil", $3); }
 | COS    '(' vse ')' { $$ = new_node_function ( cos, "cos", $3); }
 | COSH   '(' vse ')' { $$ = new_node_function ( cosh, "cosh", $3); }
 | ERFC   '(' vse ')' { $$ = new_node_function ( erfc, "erfc", $3); }
 | ERF    '(' vse ')' { $$ = new_node_function ( erf, "erf", $3); }
 | EXP    '(' vse ')' { $$ = new_node_function ( exp, "exp", $3); }
 | EXPM1  '(' vse ')' { $$ = new_node_function ( expm1, "expm1", $3); }
 | FABS   '(' vse ')' { $$ = new_node_function ( fabs, "fabs", $3); }
 | FLOOR  '(' vse ')' { $$ = new_node_function ( floor, "floor", $3); }
/* | GAMMA '(' vse ')' { $$ = new_node_function ( gamma, "gamma", $3); } */
 | J0     '(' vse ')' { $$ = new_node_function ( j0, "j0", $3); }
 | J1     '(' vse ')' { $$ = new_node_function ( j1, "j1", $3); }
 | LGAMMA '(' vse ')' { $$ = new_node_function ( lgamma, "lgamma", $3); }
 | LOG10  '(' vse ')' { $$ = new_node_function ( log10, "log10", $3); }
 | LOG1P  '(' vse ')' { $$ = new_node_function ( log1p, "log1p", $3); }
 | LOGB   '(' vse ')' { $$ = new_node_function ( logb, "logb", $3); }
 | LOG    '(' vse ')' { $$ = new_node_function ( log, "log", $3); }
 | SIN    '(' vse ')' { $$ = new_node_function ( sin, "sin", $3); }
 | SINH   '(' vse ')' { $$ = new_node_function ( sinh, "sinh", $3); }
 | SQRT   '(' vse ')' { $$ = new_node_function ( sqrt, "sqrt", $3); }
 | TAN    '(' vse ')' { $$ = new_node_function ( tan, "tan", $3); }
 | TANH   '(' vse ')' { $$ = new_node_function ( tanh, "tanh", $3); }
 | Y0     '(' vse ')' { $$ = new_node_function ( y0, "y0", $3); }
 | Y1     '(' vse ')' { $$ = new_node_function ( y1, "y1", $3); }
 | RINT   '(' vse ')' { $$ = new_node_function ( rint, "rint", $3); }
 | VARIABLE { $$ = new_node_variable ($1[0]);  }
 ;


vle: 
 '(' vle ')' { }
 | vse LE cse  { $$ = new_node_operation ($1, enum_op_le,new_node_scalar($3)); }
 | cse LE vse  { $$ = new_node_operation (new_node_scalar($1), enum_op_le,$3); }
 | vse LE vse  { $$ = new_node_operation ($1, enum_op_le,$3); }
 | vse GE cse  { $$ = new_node_operation ($1, enum_op_ge,new_node_scalar($3)); }
 | cse GE vse  { $$ = new_node_operation (new_node_scalar($1), enum_op_ge,$3); }
 | vse GE vse  { $$ = new_node_operation ($1, enum_op_ge,$3); }
 | vse '<' cse { $$ = new_node_operation ($1, enum_op_lt,new_node_scalar($3)); }
 | cse '<' vse { $$ = new_node_operation (new_node_scalar($1), enum_op_lt,$3); }
 | vse '<' vse { $$ = new_node_operation ($1, enum_op_lt,$3); }
 | vse '>' cse { $$ = new_node_operation ($1, enum_op_gt,new_node_scalar($3)); }
 | cse '>' vse { $$ = new_node_operation (new_node_scalar($1), enum_op_gt,$3); }
 | vse '>' vse { $$ = new_node_operation ($1, enum_op_gt,$3); }
 | vse EQ cse  { $$ = new_node_operation ($1, enum_op_eq,new_node_scalar($3)); }
 | cse EQ vse  { $$ = new_node_operation (new_node_scalar($1), enum_op_eq,$3); }
 | vse EQ vse  { $$ = new_node_operation ($1, enum_op_eq,$3); }
 | vse NE cse  { $$ = new_node_operation ($1, enum_op_ne,new_node_scalar($3)); }
 | cse NE vse  { $$ = new_node_operation (new_node_scalar($1), enum_op_ne,$3); }
 | vse NE vse  { $$ = new_node_operation ($1, enum_op_ne,$3); }
 | vle OR cle  { $$ = new_node_operation ($1, enum_op_or,new_node_logical($3)); }
 | cle OR vle  { $$ = new_node_operation (new_node_logical($1), enum_op_or,$3); }
 | vle OR vle  { $$ = new_node_operation ($1, enum_op_or,$3); }
 | vle AND cle { $$ = new_node_operation ($1, enum_op_and,new_node_logical($3)); }
 | cle AND vle { $$ = new_node_operation (new_node_logical($1), enum_op_and,$3); }
 | vle AND vle { $$ = new_node_operation ($1, enum_op_and,$3); }
;



%%

struct param_type * 
cello_parameters_read(FILE * fp)
{
  /* initialize the linked list with an initial sentinel (sentinel) node */
  param_head = param_curr = new_param_sentinel();
  
  yyrestart(fp);
  yyparse();
  param_head = reverse_param(param_head);
  return param_head;
}

void indent (int level)
{
  int i;
  for (i=0; i<level; i++) {
    printf ("  "); 
  }
}

void print_expression (struct node_expr * node)
{
  if (node == NULL) {
    printf ("NULL");
  } else {
    switch (node->type) {
    case enum_node_integer:
      printf ("%d",node->integer_value);
      break;
    case enum_node_scalar:
      printf ("%g",node->scalar_value);
      break;
    case enum_node_variable:
      printf ("%c",node->var_value);
      break;
    case enum_node_function:
      printf ("%s(",node->function_name); fflush(stdout);
      print_expression(node->left);
      printf (")"); fflush(stdout);
      break;
    case enum_node_operation:
      printf ("(");  fflush(stdout);
      print_expression(node->left);
      printf (") %s (",op_name[node->op_value]); fflush(stdout);
      print_expression(node->right);
      printf (")"); fflush(stdout);
      break;
    default:
      break;
    }
  }

}

void cello_parameters_print_list(struct param_type * head, int level)
{
  struct param_type * p = head->next;
  int count = 0;
  while (p && p->type != enum_parameter_sentinel && count++ < 100) {
/*     printf ("%p %s\n",p,type_name[p->type]); */
    if (p->group != NULL) {
      indent(level);
      printf ("%s %s:%s:%s = ", 
	      parameter_name[p->type],p->group, p->subgroup, p->parameter);
    } else {
      /* list element */
      indent(level);
      printf ("%s %s = ", 
	      parameter_name[p->type], p->parameter);
    }
    switch (p->type) {
    case enum_parameter_scalar:  
      printf ("%g\n",p->scalar_value);  
      break;
    case enum_parameter_integer: 
      printf ("%d\n",p->integer_value); 
      break;
    case enum_parameter_string:  
      printf ("%s\n",p->string_value); 
      break;
    case enum_parameter_logical:
      printf ("%s\n",p->logical_value ? "true" : "false");
      break;
    case enum_parameter_list:    
      indent(level);
      printf ("[\n"); 
      cello_parameters_print_list(p->list_value, level + 1);
      indent(level);
      printf ("]\n"); 
      break;
    case enum_parameter_logical_expr:
      indent(level);
      print_expression(p->op_value); printf ("\n");
      break;
    case enum_parameter_scalar_expr:
      indent(level);
      print_expression(p->op_value); printf ("\n");
      break;
    default: 
      indent(level);
      printf ("unknown type\n"); 
      break;
    }
    p = p->next;
  }
}

void cello_parameters_print()
{
  cello_parameters_print_list(param_head,0);
}

