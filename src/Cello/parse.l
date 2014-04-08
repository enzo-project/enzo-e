/* See LICENSE_CELLO file for license and copyright information */

%option noyywrap

%x IFILE

D [0-9]
LN [A-Za-z0-9_]
LL [a-z]
LU [A-Z]
LLU [A-Za-z]
E [Ee][+-]?{D}+
S [-]?

%{

#include "parse.tab.h"

 int    lineno       = 1;


char linebuf[500];

  struct bufstack {
    /* Adopted from "flex & bison", John Levine, O'Reilly Pub,  pp.28-30 */
    struct bufstack * prev;
    YY_BUFFER_STATE bs;
    int lineno;
    const char *filename;
    FILE *fp;
  } * curbs = 0;

  const char * curfilename;
  int cello_new_file (const char *);
  int cello_pop_file ();
    

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
[xyzt]    { yylval.string_type = strdup (yytext);
            return VARIABLE; }

{LLU}{LN}* { yylval.string_type = strdup (yytext);
            return IDENTIFIER; }

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
"}"        { return '}'; } /* end group */
";"        { return ';'; } /* parameter assignment separator */
","        { return ','; } /* list element separator */

include[ \t]*[\"<] { BEGIN IFILE; }
"#".*      {  /* ignore comments */ }
<IFILE>[^>\"]+ {
  { 
    int c;
    while ((c = input()) && c != '\n');
  }
  cello_new_file (strdup(yytext));
  BEGIN INITIAL;
 }
<IFILE>.|\n {
  fprintf (stderr, "%s:%d bad include line\n",
	   curfilename, yylineno);
  BEGIN INITIAL;
 }
<<EOF>> {
  if (!cello_pop_file()) yyterminate(); 
  }



{S}{D}+{E}            { yylval.float_type = atof (yytext); return FLOAT; }
{S}{D}*"."{D}+({E})?  { yylval.float_type = atof (yytext); return FLOAT; }
{S}{D}+"."{D}*({E})?  { yylval.float_type = atof (yytext); return FLOAT; }
{S}{D}+               { yylval.integer_type = atoi (yytext); return INTEGER; }

\"(\\.|[^\\"])*\"  { yylval.string_type = strdup(yytext+1); /* copy, skipping open quote " */
                     yylval.string_type [yyleng-2] = '\0';   /* remove close quote */
                     return STRING; }

[ ] {  }
^[ \t]* {  }

[^ \v\f\t\n]  { /* ignore white space */ }
\n { lineno++; }

. {  }

\n.* {strncpy (linebuf,yytext+1,sizeof(linebuf));
yyless(1);
}

%%

void yyerror(char *s)
{
  fprintf(stderr, "\n\nERROR: %s at '%s' in line %d:\n\n   \"%s\"\n\n\n",
		  s, yytext,lineno,linebuf);
  exit(1);
}

int cello_new_file(const char * filename)
/* Adopted from "flex & bison", John Levine, O'Reilly Pub,  pp.28-30 */
{
  FILE * fp = fopen (filename, "r");
  struct bufstack * bs = malloc (sizeof(struct bufstack));

  if (fp == NULL) {
    printf ("%s:%d Parse Error: cannot open file %s\n",
	    __FILE__,__LINE__,filename);
    exit(1);
  }

  if (curbs != NULL) curbs->lineno = yylineno;
  bs->prev = curbs;
  bs->bs = yy_create_buffer(fp,YY_BUF_SIZE);
  bs->fp = fp;
  bs->filename = filename;
  yy_switch_to_buffer(bs->bs);
  curbs = bs;
  yylineno = 1;
  curfilename = filename;
  return 1;
}

int cello_pop_file ()
/* Adopted from "flex & bison", John Levine, O'Reilly Pub,  pp.28-30 */
{
  struct bufstack * bs = curbs;
  struct bufstack * prevbs;

  if (bs == NULL) return 0;
  fclose (bs->fp);
  yy_delete_buffer(bs->bs);
  prevbs = bs->prev;
  free(bs);

  /* BUG FIX: moved curbs assignment before conditional return */

  curbs = prevbs;

  if (prevbs == 0) return 0;

  yy_switch_to_buffer(prevbs->bs);
  yylineno = curbs->lineno;
  curfilename = curbs->filename;
  return 1;
}

