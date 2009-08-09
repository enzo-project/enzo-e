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
 int line_count = 0;

%}

/* double acos(double); */
/* double acosh(double); */
/* double asin(double); */
/* double asinh(double); */
/* double atan(double); */
/* double atanh(double); */
/* double cbrt(double); */
/* double ceil(double); */
/* double cos(double); */
/* double cosh(double); */
/* double erfc(double); */
/* double erf(double); */
/* double exp(double); */
/* double expm1(double); */
/* double fabs(double); */
/* double floor(double); */
/* double gamma(double); */
/* double j0(double); */
/* double j1(double); */
/* double lgamma(double); */
/* double log10(double); */
/* double log1p(double); */
/* double logb(double); */
/* double log(double); */
/* double sin(double); */
/* double sinh(double); */
/* double sqrt(double); */
/* double tan(double); */
/* double tanh(double); */
/* double y0(double); */
/* double y1(double); */
/* int    ilogb(double); */
/* int    isnan(double); */
/* double rint(double); */
/*  */
/*  */
/* double atan2(double, double); */
/* double fmod(double, double); */
/* double hypot(double, double); */
/* double jn(int, double); */
/* double ldexp(double, int); */
/* double nextafter(double, double); */
/* double pow(double, double); */
/* double remainder(double, double); */
/* double scalb(double, double); */
/* double yn(int, double); */

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
"#".*      { printf (""); /* ignore comments */ }

[xyzt]             { printf (" CONSTANT "); }
{D}+{E}*           { printf (" SCALAR "); }
{D}*"."{D}+({E})?  { printf (" SCALAR "); }
{D}+"."{D}*({E})?  { printf (" SCALAR "); }

\"(\\.|[^\\"])*\"   { printf (" STRING "); }

[ ] { printf (" "); }
^[ \t]* { printf (""); }

[^ \v\f\t\n]  { /* ignore white space */ }
\n { line_count++; printf ("\n"); }

. { printf (""); }

%%

main()
{
    yylex() ;
}
