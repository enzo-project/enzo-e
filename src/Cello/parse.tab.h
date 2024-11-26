/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

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

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_YY_BUILD_CELLO_PARSE_TAB_H_INCLUDED
# define YY_YY_BUILD_CELLO_PARSE_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    STRING = 258,                  /* STRING  */
    IDENTIFIER = 259,              /* IDENTIFIER  */
    VARIABLE = 260,                /* VARIABLE  */
    FLOAT = 261,                   /* FLOAT  */
    INTEGER = 262,                 /* INTEGER  */
    LOGICAL = 263,                 /* LOGICAL  */
    LE = 264,                      /* LE  */
    GE = 265,                      /* GE  */
    NE = 266,                      /* NE  */
    EQ = 267,                      /* EQ  */
    AND = 268,                     /* AND  */
    OR = 269,                      /* OR  */
    ACOS = 270,                    /* ACOS  */
    ACOSH = 271,                   /* ACOSH  */
    APPEND = 272,                  /* APPEND  */
    ASIN = 273,                    /* ASIN  */
    ASINH = 274,                   /* ASINH  */
    ATAN = 275,                    /* ATAN  */
    ATANH = 276,                   /* ATANH  */
    CBRT = 277,                    /* CBRT  */
    CEIL = 278,                    /* CEIL  */
    COS = 279,                     /* COS  */
    COSH = 280,                    /* COSH  */
    ERFC = 281,                    /* ERFC  */
    ERF = 282,                     /* ERF  */
    EXP = 283,                     /* EXP  */
    EXPM1 = 284,                   /* EXPM1  */
    FABS = 285,                    /* FABS  */
    FLOOR = 286,                   /* FLOOR  */
    J0 = 287,                      /* J0  */
    J1 = 288,                      /* J1  */
    LGAMMA = 289,                  /* LGAMMA  */
    LOG10 = 290,                   /* LOG10  */
    LOG1P = 291,                   /* LOG1P  */
    LOGB = 292,                    /* LOGB  */
    LOG = 293,                     /* LOG  */
    PI = 294,                      /* PI  */
    SIN = 295,                     /* SIN  */
    SINH = 296,                    /* SINH  */
    SQRT = 297,                    /* SQRT  */
    TAN = 298,                     /* TAN  */
    TANH = 299,                    /* TANH  */
    Y0 = 300,                      /* Y0  */
    Y1 = 301,                      /* Y1  */
    RINT = 302                     /* RINT  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 381 "build/Cello/parse.y"
 
  int logical_type;  
  int integer_type; 
  double float_type;  
  char * string_type; 
  char * group_type;
  struct node_expr * node_type;
  

#line 121 "build/Cello/parse.tab.h"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;


int yyparse (void);


#endif /* !YY_YY_BUILD_CELLO_PARSE_TAB_H_INCLUDED  */
