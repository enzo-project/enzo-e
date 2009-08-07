D [0-9]
L [A-Za-z_]
LL [a-z]
LU [A-Z]
E [Ee][+-]?{D}+

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
%}

%%

false     { printf (" FALSE "); } 
true      { printf (" TRUE "); } 
{LU}{L}+  { printf (" GROUP "); }  /*  NAME_GROUP */
{LL}{L}+  { printf (" PARAMETER "); } /* PARAMETER_NAME */

"="        { printf ("="); }         /* = */

"*"        { printf ("*"); }
"/"        { printf ("/"); }
"+"        { printf ("+"); }
"-"        { printf ("-"); }

"<"        { printf ("<"); }         /* = */
">"        { printf (">"); }
"<="        { printf ("LE"); }
">="        { printf ("GE"); }
"=="        { printf ("EQ"); }
"!="        { printf ("NE"); }
"&&"        { printf ("AND"); }
"||"        { printf ("OR"); }
"("          { printf ("("); }
")"          { printf (")"); }

"["        { printf ("["); } /* begin list */
"]"        { printf ("]"); } /* end list */
"{"        { printf ("{"); } /* begin group */
"}"        { printf ("}"); } /* end group ( and end assignment ) */
";"        { printf (";"); } /* end assignment */
"#".*"\n"      { printf ("\n"); /* ignore comments */ }

[xyzt]             { printf (" CONSTANT "); }
{D}+{E}*           { printf (" SCALAR "); }
{D}*"."{D}+({E})?  { printf (" SCALAR "); }
{D}+"."{D}*({E})?  { printf (" SCALAR "); }

\"(\\.|[^\\"])*\"   { printf (" STRING "); }

[ ] { printf (" "); }
^[ \t]*"\n" { printf (""); }

[^ \v\f\t\n]  { /* ignore white space */ }

. { printf (""); }

%%

main()
{
    yylex() ;
}
