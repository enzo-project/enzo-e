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
#include <stdlib.h>

#ifndef __APPLE__
	#include <malloc.h>
#endif
	
#define YYDEBUG 1

  /* Quiet a few -Wall errors */

int yylex (void);
void yyrestart  (FILE * input_file );
void yyerror(char *s);
void yylex_destroy();

#include "parse.tab.h"

#include "parse.h"

const char * node_name[] = {
  "node_unknown",
  "node_operation",
  "node_float",
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

  /* ANY CHANGES HERE MUST BE REFLECTED IN parse.h enum_parameter[] */
  const char * parameter_name[]  = {
    "unknown",
    "sentinel",
    "group",
    "integer",
    "float",
    "string",
    "identifier",
    "logical",
    "list",
    "float_expr",
    "logical_expr",
    "function" };

  /* Structure for storing a single parameter / value pair in a linked list */


  struct node_expr * new_node_operation
    (struct node_expr * left, 
     enum enum_op oper,
     struct node_expr * right)
  {
    
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_operation;
    node->op_value      = oper;
    node->left          = left;
    node->right         = right;
    node->function_name = NULL;
    return node;
  }

  struct node_expr * new_node_float (double value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_float;
    node->float_value  = value;
    node->left          = NULL;
    node->right         = NULL;
    node->function_name = NULL;
    return node;
  }
  struct node_expr * new_node_logical (int value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_integer;
    node->integer_value = value;
    node->left          = NULL;
    node->right         = NULL;
    node->function_name = NULL;
    return node;
  }
  struct node_expr * new_node_variable (char * value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_variable;
    node->var_value     = value[0];
    node->left          = NULL;
    node->right         = NULL;
    node->function_name = NULL;
    free (value);
    return node;
  }
  struct node_expr * new_node_function
    (double (*function)(double),
     char * function_name,
     struct node_expr * argument)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_function;
    node->fun_value     = function;
    node->left          = argument;
    node->right         = NULL;
    node->function_name = strdup(function_name);
    return node;
  }


  /* The head of the linked list of parameter / value pairs */

  struct param_struct * param_head = NULL; /* head of entire list */
  struct param_struct * param_curr = NULL; /* head of current list */

  /* The current groups and parameter type */


  char *              current_parameter = NULL;
  char *              current_group[MAX_GROUP_DEPTH];
  int                 current_group_level = 0;
  enum enum_parameter current_type      = enum_parameter_sentinel;

  void clear_groups (char * groups[]) {
    int i;
    for (i=0; i<MAX_GROUP_DEPTH; i++) {
      groups[i] = 0; 
    }
  };


  void copy_groups (char * group_dest[], char * group_src[]) {
    int i;
    for (i=0; i<MAX_GROUP_DEPTH; i++) {
      /* MEMORY LEAK */
      group_dest[i] = (group_src[i]) ? strdup(group_src[i]) : 0;
    }
  };

  /* Function to update parameter's groups once the group is known */

/*   void update_group (char * group) */
/*     { */
/*       struct param_struct * p = param_curr; */
/*       while (p->next->type  != enum_parameter_sentinel &&  */
/* 	     p->next->group == NULL) { */
/* 	p->next->group = strdup(group); */
/*         p = p -> next; */
/*       } */
/*     } */

  /* Insert a parameter into the list */

  void insert_param(struct param_struct * head, struct param_struct * new)
  {
     new->next  = head->next;
     head->next = new;
  }

  /* Delete a parameter from the list given a pointer to the previous element */

  void delete_param(struct param_struct * previous)
  {
    struct param_struct * item = previous->next;
    previous->next = item->next;
    free (item);     
  }

  /* Function to update parameter's subgroups once the subgroup is known */

/*   void update_subgroup (char * subgroup) */
/*     { */
/*       struct param_struct * p = param_curr; */
/*       int inside_subgroup = 1; */
/*       while (p->next->type     != enum_parameter_sentinel &&  */
/* 	     p->next->subgroup == NULL) { */
/* 	if (p->next->type == enum_parameter_subgroup) { */
/* 	  inside_subgroup = 0; */
/*           delete_param(p); */
/*         } else if (inside_subgroup) { */
/*           p->next->subgroup = strdup(subgroup); */
/*           p = p -> next; */
/*         } */
/*       } */
/*     } */

  struct param_struct * reverse_param(struct param_struct * old_head)
  {
    /* Keep sentinel the same */

    struct param_struct * new_head = old_head;

    struct param_struct * p = old_head;
    struct param_struct * c = p->next;
    struct param_struct * n = c->next;

    do {
      /* If parameter is a list, recursively reverse it as well */
      if (c->type == enum_parameter_list) {
	c->list_value = reverse_param(c->list_value);
      }
      c->next = p;
      p       = c;
      c       = n;
      n = n->next;
    } while (p->type != enum_parameter_sentinel) ;

    new_head = p;
    return new_head;
  }

  /* Function for creating and inserting a new parameter / value pair */
  /* in the linked list */

  struct param_struct * new_param ()
  {
    /* Create the new node */

    /* MEMORY LEAK */
     struct param_struct * p = 
       (struct param_struct *) malloc (sizeof (struct param_struct));

   /* Fill in the non-type-specific values for the new node */

     /* MEMORY LEAK */
     
     copy_groups(p->group,current_group);

     p->parameter = (current_parameter) ? strdup(current_parameter) : 0;

     free (current_parameter);
     current_parameter = 0;

     current_type = enum_parameter_unknown;

     insert_param(param_curr,p);

   /* Clear variables for the next assignment */

     return p;
  }

  /* New integer parameter assignment */

  void new_param_integer (int value)
  {
    struct param_struct * p = new_param();
    p->type          = enum_parameter_integer;
    p->integer_value = value;
  }


  /* New floating-point parameter assignment */

  void new_param_float (double value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_float;
    p->float_value = value;
  }

  /* New logical parameter assignment */

  void new_param_logical (int value)
  {
    struct param_struct * p = new_param();
    p->type          = enum_parameter_logical;
    p->logical_value = value;
  }

  /* New string parameter assignment */
  void new_param_string (char * value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_string;
    p->string_value = value;
  }

  /* New subgroup  */
  void new_param_group (char * value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_group;
    p->string_value = value;
  }

  /* New empty parameter assignment: FIRST NODE IN LIST IS A SENTINEL  */
  struct param_struct * new_param_sentinel ()
  {

    /* MEMORY LEAK */
    struct param_struct * p = 
      (struct param_struct *) malloc (sizeof (struct param_struct));

    clear_groups(p->group);
    p->parameter = NULL;
    p->type      = enum_parameter_sentinel;
    p->next       = p;
    p->list_value = NULL;

    return p;
  }

  /* New list parameter assignment */

  void new_param_list (struct param_struct * curr)
  {
    struct param_struct * p = new_param();
    p->type       = enum_parameter_list;
    p->list_value = curr;
  }

  /* New string parameter assignment */

  void new_param_expr (enum enum_parameter type,
		       struct node_expr * value)
  {
    struct param_struct * p = new_param();
    p->type     = type;
    p->op_value = value;
  }

  void new_parameter()
  {
     switch (current_type) {
     case enum_parameter_group:
       new_param_group(yylval.group_type);
       break;
     case enum_parameter_integer:
       new_param_integer(yylval.integer_type);
       break;
     case enum_parameter_float:
       new_param_float(yylval.float_type);
       break;
     case enum_parameter_string: 
       new_param_string(yylval.string_type);
       break;
     case enum_parameter_logical:
       new_param_logical(yylval.logical_type);
       break;
     case enum_parameter_list:
       break;
     case enum_parameter_float_expr:
       new_param_expr(enum_parameter_float_expr,yylval.node_type);
       break;
     case enum_parameter_logical_expr:
       new_param_expr(enum_parameter_logical_expr,yylval.node_type);
       break;
    default:
       printf ("%s:%d Parse Error: unknown type %d\n",
	       __FILE__,__LINE__,current_type);
       exit(1);
       break;
     }
  }
%}


%union { 
  int logical_type;  
  int integer_type; 
  double float_type;  
  char * string_type; 
  char * group_type;
  struct node_expr * node_type;
  }

/* %token <string_type>  GROUP_NAME */
%token <string_type>  STRING
%token <string_type>  IDENTIFIER
%token <string_type> VARIABLE
%token <float_type>  FLOAT
%token <integer_type> INTEGER
%token <logical_type> LOGICAL

%type <integer_type> cie
%type <float_type>  cse
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
group_name parameter_group  {  }

parameter_group :
   '{' parameter_list '}'        { current_group[--current_group_level] = 0; }
 | '{' parameter_list ';' '}'    { current_group[--current_group_level] = 0; }

parameter_list : 
                      parameter_assignment  {  }
 | parameter_list ';' parameter_assignment  {  }
 |                    group  {  }
 | parameter_list ';' group  {  }

 
group_name :
  IDENTIFIER                    { current_group[current_group_level++] = $1; }

parameter_name :
  IDENTIFIER                    { current_parameter = $1;} 

parameter_assignment : 
  parameter_name '=' parameter_value { new_parameter(); }
 ;

parameter_value : 
 STRING { current_type = enum_parameter_string;       yylval.string_type = $1; }
 | cie  { current_type = enum_parameter_integer;      yylval.integer_type = $1;}
 | cse  { current_type = enum_parameter_float;       yylval.float_type = $1;}
 | cle  { current_type = enum_parameter_logical;      yylval.logical_type = $1; }
 | vse  { current_type = enum_parameter_float_expr;  yylval.node_type = $1; }
 | vle  { current_type = enum_parameter_logical_expr; yylval.node_type = $1; }
 | list { current_type = enum_parameter_list; }
 ;

list: LIST_BEGIN list_elements LIST_END {  }
    | LIST_BEGIN LIST_END {  }

LIST_BEGIN:
 '[' { 
   struct param_struct * p = new_param_sentinel();
   p->list_value = param_curr;
   new_param_list(p);
   param_curr = p;
 }
LIST_END:
 ']' { param_curr = param_curr->list_value; }


list_elements:
                      parameter_value    { new_parameter(); }
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
| FLOAT { $$ = $1;}
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
'(' vse ')'    { $$ = $2; }
 | vse '+' cse { $$ = new_node_operation ($1, enum_op_add,new_node_float($3)); }
 | cse '+' vse { $$ = new_node_operation (new_node_float($1), enum_op_add,$3); }
 | vse '+' vse { $$ = new_node_operation ($1, enum_op_add,$3); }
 | vse '-' cse { $$ = new_node_operation ($1, enum_op_sub,new_node_float($3)); }
 | cse '-' vse { $$ = new_node_operation (new_node_float($1), enum_op_sub,$3); }
 | vse '-' vse { $$ = new_node_operation ($1, enum_op_sub,$3); }
 | vse '*' cse { $$ = new_node_operation ($1, enum_op_mul,new_node_float($3)); }
 | cse '*' vse { $$ = new_node_operation (new_node_float($1), enum_op_mul,$3); }
 | vse '*' vse { $$ = new_node_operation ($1, enum_op_mul,$3); }
 | vse '/' cse { $$ = new_node_operation ($1, enum_op_div,new_node_float($3)); }
 | cse '/' vse { $$ = new_node_operation (new_node_float($1), enum_op_div,$3); }
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
 | VARIABLE { $$ = new_node_variable ($1);  }
 ;


vle: 
 '(' vle ')' { $$ = $2; }
 | vse LE cse  { $$ = new_node_operation ($1, enum_op_le,new_node_float($3)); }
 | cse LE vse  { $$ = new_node_operation (new_node_float($1), enum_op_le,$3); }
 | vse LE vse  { $$ = new_node_operation ($1, enum_op_le,$3); }
 | vse GE cse  { $$ = new_node_operation ($1, enum_op_ge,new_node_float($3)); }
 | cse GE vse  { $$ = new_node_operation (new_node_float($1), enum_op_ge,$3); }
 | vse GE vse  { $$ = new_node_operation ($1, enum_op_ge,$3); }
 | vse '<' cse { $$ = new_node_operation ($1, enum_op_lt,new_node_float($3)); }
 | cse '<' vse { $$ = new_node_operation (new_node_float($1), enum_op_lt,$3); }
 | vse '<' vse { $$ = new_node_operation ($1, enum_op_lt,$3); }
 | vse '>' cse { $$ = new_node_operation ($1, enum_op_gt,new_node_float($3)); }
 | cse '>' vse { $$ = new_node_operation (new_node_float($1), enum_op_gt,$3); }
 | vse '>' vse { $$ = new_node_operation ($1, enum_op_gt,$3); }
 | vse EQ cse  { $$ = new_node_operation ($1, enum_op_eq,new_node_float($3)); }
 | cse EQ vse  { $$ = new_node_operation (new_node_float($1), enum_op_eq,$3); }
 | vse EQ vse  { $$ = new_node_operation ($1, enum_op_eq,$3); }
 | vse NE cse  { $$ = new_node_operation ($1, enum_op_ne,new_node_float($3)); }
 | cse NE vse  { $$ = new_node_operation (new_node_float($1), enum_op_ne,$3); }
 | vse NE vse  { $$ = new_node_operation ($1, enum_op_ne,$3); }
 | vle OR cle  { $$ = new_node_operation ($1, enum_op_or,new_node_logical($3)); }
 | cle OR vle  { $$ = new_node_operation (new_node_logical($1), enum_op_or,$3); }
 | vle OR vle  { $$ = new_node_operation ($1, enum_op_or,$3); }
 | vle AND cle { $$ = new_node_operation ($1, enum_op_and,new_node_logical($3)); }
 | cle AND vle { $$ = new_node_operation (new_node_logical($1), enum_op_and,$3); }
 | vle AND vle { $$ = new_node_operation ($1, enum_op_and,$3); }
;

%%

int cello_new_file(const char * filename);

struct param_struct * 
cello_parameters_read(const char * filename, FILE * fp)
{
  clear_groups(current_group);

  /* initialize the linked list with an initial sentinel (sentinel) node */
  param_head = param_curr = new_param_sentinel();

  /*   yydebug=1; */

  cello_new_file (filename);

  yyrestart(fp);

  yyparse();
  yylex_destroy();

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

void print_expression (struct node_expr * node,
		       FILE * fp)
{
  if (node == NULL) {
    fprintf (fp,"NULL");
  } else {
    char left,right;
    switch (node->type) {
    case enum_node_integer:
      fprintf (fp,"%d",node->integer_value);
      break;
    case enum_node_float:
      /* '#' format character forces a decimal point */
      fprintf (fp,FLOAT_FORMAT,node->float_value);
      break;
    case enum_node_variable:
      fprintf (fp,"%c",node->var_value);
      break;
    case enum_node_function:
      fprintf (fp,"%s(",node->function_name);
      print_expression(node->left,fp);
      fprintf (fp,")");
      break;
    case enum_node_operation:
      left  = (node->left->type == enum_node_operation) ? '(' : ' ';
      right = (node->left->type == enum_node_operation) ? ')' : ' ';
      fprintf (fp,"%c",left);
      print_expression(node->left,fp);
      fprintf (fp,"%c",right);
      fprintf (fp," %s ",op_name[node->op_value]);
      left  = (node->right->type == enum_node_operation) ? '(' : ' ';
      right = (node->right->type == enum_node_operation) ? ')' : ' ';
      fprintf (fp,"%c",left);
      print_expression(node->right,fp);
      fprintf (fp,"%c",right);
      break;
    default:
      break;
    }
    fflush(fp);
  }

}

void sprintf_expression (struct node_expr * node,
			 char * buffer)
/* WARNING: buffer is assumed to be big enough to hold the expression */
{
  if (node == NULL) {
    sprintf (buffer,"NULL");
  } else {
    char left,right;
    switch (node->type) {
    case enum_node_integer:
      sprintf (buffer,"%d",node->integer_value);
      buffer += strlen(buffer);
      break;
    case enum_node_float:
      /* '#' format character forces a decimal point */
      sprintf (buffer,FLOAT_FORMAT,node->float_value);
      buffer += strlen(buffer);
      break;
    case enum_node_variable:
      sprintf (buffer,"%c",node->var_value);
      buffer += strlen(buffer);
      break;
    case enum_node_function:
      sprintf (buffer,"%s(",node->function_name);
      buffer += strlen(buffer);
      sprintf_expression(node->left,buffer+strlen(buffer));
      buffer += strlen(buffer);
      sprintf (buffer,")");
      buffer += strlen(buffer);
      break;
    case enum_node_operation:
      left  = (node->left->type == enum_node_operation) ? '(' : ' ';
      right = (node->left->type == enum_node_operation) ? ')' : ' ';
      sprintf (buffer,"%c",left);
      buffer += strlen(buffer);
      sprintf_expression(node->left,buffer+strlen(buffer));
      buffer += strlen(buffer);
      sprintf (buffer,"%c",right);
      buffer += strlen(buffer);
      sprintf (buffer," %s ",op_name[node->op_value]);
      buffer += strlen(buffer);
      left  = (node->right->type == enum_node_operation) ? '(' : ' ';
      right = (node->right->type == enum_node_operation) ? ')' : ' ';
      sprintf (buffer,"%c",left);
      buffer += strlen(buffer);
      sprintf_expression(node->right,buffer+strlen(buffer));
      buffer += strlen(buffer);
      sprintf (buffer,"%c",right);
      buffer += strlen(buffer);
      break;
    default:
      break;
    }
  }
}

void cello_parameters_print_list(struct param_struct * head, int level)
{
  struct param_struct * p = head->next;
  int i;

  while (p && p->type != enum_parameter_sentinel) {

    if (p->group != NULL) {
      indent(level);
      printf ("%s ", parameter_name[p->type]);
      for (i=0; p->group[i] != NULL && i < MAX_GROUP_DEPTH; i++) {
	printf ("%s:",p->group[i]);
      }
      printf ("%s = ", p->parameter);
    } else {
      /* list element */
      indent(level);
      printf ("%s %s = ", 
	      parameter_name[p->type], p->parameter);
    }
    switch (p->type) {
    case enum_parameter_float:  
      /* '#' format character forces a decimal point */
      printf (FLOAT_FORMAT,p->float_value);  
      break;
    case enum_parameter_integer: 
      printf ("%d\n",p->integer_value); 
      break;
    case enum_parameter_string:  
      printf ("%s\n",p->string_value); 
      break;
    case enum_parameter_group:  
      printf ("Uh oh: GROUP %s (should be deleted)\n",p->string_value);
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
      print_expression(p->op_value,stdout); printf ("\n");
      break;
    case enum_parameter_float_expr:
      indent(level);
      print_expression(p->op_value,stdout); printf ("\n");
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

