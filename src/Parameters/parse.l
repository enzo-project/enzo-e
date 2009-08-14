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

#include "parse.tab.h"

 double scalar_value = 0;

%}

/* double foo (double) */

acos   { return ACOS; }
acosh  { return ACOSH; }
asin   { return ASIN; }
asinh  { return ASINH; }
atan   { return ATAN; }
atanh  { return ATANH; }
cbrt   { return CBRT; }
ceil   { return CEIL; }
cos    { return COS; }
cosh   { return COSH; }
erfc   { return ERFC; }
erf    { return ERF; }
exp    { return EXP; }
expm1  { return EXPM1; }
fabs   { return FABS; }
floor  { return FLOOR; }
gamma  { return GAMMA; }
j0     { return J0; }
j1     { return J1; }
lgamma { return LGAMMA; }
log10  { return LOG10; }
log1p  { return LOG1P; }
logb   { return LOGB; }
log    { return LOG; }
sin    { return SIN; }
sinh   { return SINH; }
sqrt   { return SQRT; }
tan    { return TAN; }
tanh   { return TANH; }
y0     { return Y0; }
y1     { return Y1; }
rint   { return RINT; }

/* int foo (double) */

/* int    ilogb(double); */
/* int    isnan(double); */

/*  */
/* double foo (double,double) */

atan2     { return ATAN2; }
fmod      { return FMOD; }
hypot     { return HYPOT; }
nextafter { return NEXTAFTER; }
pow       { return POW; }
remainder { return REMAINDER; }
scalb     { return SCALB; }

/* double jn(int, double); */
/* double ldexp(double, int); */
/* double yn(int, double); */

%%

false     { return FALSE; } 
true      { return TRUE; } 
{LU}{L}+  { return GROUP; }  /*  NAME_GROUP */
{LL}{L}+  { return PARAMETER; } /* PARAMETER_NAME */

"="        { return '='; }         /* = */

"*"        { return '*'; }
"/"        { return '/'; }
"+"        { return '+'; }
"-"        { return '-'; }

"<"        { return '<'; }         /* = */
">"        { return '>'; }
"<="        { return LE; }
">="        { return GE; }
"=="        { return EQ; }
"!="        { return NE; }
"&&"        { return AND; }
"||"        { return OR; }
"("          { return '('; }
")"          { return ')'; }

"["        { return '['; } /* begin list */
"]"        { return ']'; } /* end list */
"{"        { return '{'; } /* begin group */
"}"        { return '}'; } /* end group ( and end assignment ) */
";"        { return ';'; } /* end assignment */
"#".*      {  /* ignore comments */ }

[xyzt]             { return CONSTANT; }
{D}+{E}*           { scalar_value = atof (yytext); return SCALAR; }
{D}*"."{D}+({E})?  { scalar_value = atof (yytext); return SCALAR; }
{D}+"."{D}*({E})?  { scalar_value = atof (yytext); return SCALAR; }

\"(\\.|[^\\"])*\"   { return STRING; }

[ ] {  }
^[ \t]* {  }

[^ \v\f\t\n]  { /* ignore white space */ }
\n { yylineno++; }

. {  }

%%

