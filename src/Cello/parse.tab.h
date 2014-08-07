/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_BUILD_CELLO_PARSE_TAB_H_INCLUDED
# define YY_YY_BUILD_CELLO_PARSE_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    STRING = 258,
    IDENTIFIER = 259,
    VARIABLE = 260,
    FLOAT = 261,
    INTEGER = 262,
    LOGICAL = 263,
    LE = 264,
    GE = 265,
    NE = 266,
    EQ = 267,
    AND = 268,
    OR = 269,
    ACOS = 270,
    ACOSH = 271,
    ASIN = 272,
    ASINH = 273,
    ATAN = 274,
    ATANH = 275,
    CBRT = 276,
    CEIL = 277,
    COS = 278,
    COSH = 279,
    ERFC = 280,
    ERF = 281,
    EXP = 282,
    EXPM1 = 283,
    FABS = 284,
    FLOOR = 285,
    J0 = 286,
    J1 = 287,
    LGAMMA = 288,
    LOG10 = 289,
    LOG1P = 290,
    LOGB = 291,
    LOG = 292,
    SIN = 293,
    SINH = 294,
    SQRT = 295,
    TAN = 296,
    TANH = 297,
    Y0 = 298,
    Y1 = 299,
    RINT = 300
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 365 "build/Cello/parse.y" /* yacc.c:1909  */
 
  int logical_type;  
  int integer_type; 
  double float_type;  
  char * string_type; 
  char * group_type;
  struct node_expr * node_type;
  

#line 110 "build/Cello/parse.tab.h" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_BUILD_CELLO_PARSE_TAB_H_INCLUDED  */
