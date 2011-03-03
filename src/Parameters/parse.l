%option noyywrap

D [0-9]
LN [A-Za-z0-9_]
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

 void yyerror(char *s);

 int    lineno       = 1;


%}

/* double foo (double) */


/* int foo (double) */

/* int    ilogb(double); */
/* int    isnan(double); */

/* double foo (double,double) */

/* atan2     { return ATAN2; } */
/* fmod      { return FMOD; } */
/* hypot     { return HYPOT; } */
/* nextafter { return NEXTAFTER; } */
/* pow       { return POW; } */
/* remainder { return REMAINDER; } */
/* scalb     { return SCALB; } */

/* double jn(int, double); */
/* double ldexp(double, int); */
/* double yn(int, double); */

%%


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

false     { yylval.logical_type = 0; return LOGICAL; } 
true      { yylval.logical_type = 1; return LOGICAL; } 
{LU}{LN}+ { /* free (yylval.string_type); ILLEGAL */
            yylval.string_type = strdup (yytext); /* MEMORY LEAK */
            return GROUP_NAME; } /* Group names begin with upper-case */
[xyzt]    { yylval.string_type = strdup (yytext);
            return VARIABLE; }
{LL}{LN}* { yylval.string_type = strdup (yytext); /* MEMORY LEAK */
            return IDENTIFIER; } /* group, parameter, and variables */
                                 /* begin with lower-case */

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
";"        { return ';'; } /* parameter assignment separator */
","        { return ','; } /* list element separator */
"#".*      {  /* ignore comments */ }

{D}+{E}            { yylval.scalar_type = atof (yytext); return SCALAR; }
{D}*"."{D}+({E})?  { yylval.scalar_type = atof (yytext); return SCALAR; }
{D}+"."{D}*({E})?  { yylval.scalar_type = atof (yytext); return SCALAR; }
{D}+               { yylval.integer_type = atoi (yytext); return INTEGER; }

\"(\\.|[^\\"])*\"  { yylval.string_type = strdup(yytext+1); /* copy, skipping open quote */
                     yylval.string_type [yyleng-2] = '\0';   /* remove close quote */
                     return STRING; }

[ ] {  }
^[ \t]* {  }

[^ \v\f\t\n]  { /* ignore white space */ }
\n { lineno++; }

. {  }

%%

void yyerror(char *s)
{
  fprintf(stderr, "\n\nERROR: %s in line %d\n\n\n", s,lineno);
  exit(1);
}
