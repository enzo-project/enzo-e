%{
#include <stdio.h>
#include "parse.tab.h"
  int yydebug;
#define YYDEBUG 1
%}

%token GROUP
%token PARAMETER

/* %union {int logical_type; double scalar_type; char string_type[100];} */

%token LOGICAL
%token SCALAR
%token STRING

/* %token <logical_type> LOGICAL */
/* %token <scalar_type> SCALAR */
/* %token <string_type> STRING */
/* %type <scalar_type> scalar_expression */
/* %type <logical_type> logical */

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
 | file named_group { }
 ;

named_group: group_name group { }
;

group :  '{' parameter_list '}' 
 {  }
 ;

parameter_list : 
   parameter_assignment 
 | parameter_assignment  ';' parameter_list
{  }
 ;

parameter_assignment : parameter_name '=' parameter_value
  {  }
 ;

parameter_name : PARAMETER 
   {  }
 ;

parameter_value : STRING | scalar_expression | logical | list | group
{  }
;

list: '[' list_elements ']'


list_elements:
  parameter_value 
| parameter_value ',' list_elements
{ }
;

/* logical_expression: '(' logical ')' */
/*  { } */

logical: 
 '(' logical ')' { }
 | scalar_expression LE scalar_expression { }
 | scalar_expression GE scalar_expression { }
 | scalar_expression '<' scalar_expression { }
 | scalar_expression '>' scalar_expression { }
 | scalar_expression EQ scalar_expression { }
 | scalar_expression NE scalar_expression { }
 | logical OR logical {  }
 | logical AND logical {  }
 | LOGICAL { }
;

group_name : GROUP 
{  }
 ;


scalar_expression: 
   scalar_expression '+' scalar_expression {  }
 | scalar_expression '-' scalar_expression {  }
 | scalar_expression '*' scalar_expression {  }
 | scalar_expression '/' scalar_expression {  }
 | SCALAR { }
 | CONSTANT { }
 ;


%%
main(int argc, char **argv)
{
  yydebug=1;
  yyparse();
}

yyerror(char *s)
{
  fprintf(stderr, "error: %s\n", s);
}
