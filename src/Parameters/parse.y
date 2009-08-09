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

%%

scalar_expression: scalar_factor
 | scalar_expression '+' scalar_factor { $$ = $1 + $3; }
 | scalar_expression '-' scalar_factor { $$ = $1 - $3; }
 ;

scalar_factor: scalar_term
 | scalar_factor      '*' scalar_term   { $$ = $1 * $3; }
 | scalar_factor      '/' scalar_term   { $$ = $1 / $3; }
 ;

