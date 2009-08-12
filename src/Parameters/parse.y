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

%token SIN

%%

file : /* nothing */
 | file group { printf ("file group EOF\n") }
 ;

group : group_name '{' parameter_list '}' { printf ("parameter list\n"); }
 ;

group_name : GROUP { printf ("group\n"); }
 ;

parameter_list : /* nothing */
 | parameter_list parameter_assignment ';'
 ;

parameter_assignment : parameter_name '=' parameter_value
  { printf ("parameter_assignment\n"); }
 ;

parameter_name : PARAMETER { printf ("PARAMETER\n"); }
 ;

parameter_value : STRING | scalar_expression

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
  fprintf(stderr, "error: %s on line %d\n", s, line_count);
}
