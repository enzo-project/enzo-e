/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     GROUP_NAME = 258,
     STRING = 259,
     SCALAR = 260,
     INTEGER = 261,
     LOGICAL = 262,
     IDENTIFIER = 263,
     VARIABLE = 264,
     LE = 265,
     GE = 266,
     NE = 267,
     EQ = 268,
     AND = 269,
     OR = 270,
     ACOS = 271,
     ACOSH = 272,
     ASIN = 273,
     ASINH = 274,
     ATAN = 275,
     ATANH = 276,
     CBRT = 277,
     CEIL = 278,
     COS = 279,
     COSH = 280,
     ERFC = 281,
     ERF = 282,
     EXP = 283,
     EXPM1 = 284,
     FABS = 285,
     FLOOR = 286,
     J0 = 287,
     J1 = 288,
     LGAMMA = 289,
     LOG10 = 290,
     LOG1P = 291,
     LOGB = 292,
     LOG = 293,
     SIN = 294,
     SINH = 295,
     SQRT = 296,
     TAN = 297,
     TANH = 298,
     Y0 = 299,
     Y1 = 300,
     RINT = 301
   };
#endif
/* Tokens.  */
#define GROUP_NAME 258
#define STRING 259
#define SCALAR 260
#define INTEGER 261
#define LOGICAL 262
#define IDENTIFIER 263
#define VARIABLE 264
#define LE 265
#define GE 266
#define NE 267
#define EQ 268
#define AND 269
#define OR 270
#define ACOS 271
#define ACOSH 272
#define ASIN 273
#define ASINH 274
#define ATAN 275
#define ATANH 276
#define CBRT 277
#define CEIL 278
#define COS 279
#define COSH 280
#define ERFC 281
#define ERF 282
#define EXP 283
#define EXPM1 284
#define FABS 285
#define FLOOR 286
#define J0 287
#define J1 288
#define LGAMMA 289
#define LOG10 290
#define LOG1P 291
#define LOGB 292
#define LOG 293
#define SIN 294
#define SINH 295
#define SQRT 296
#define TAN 297
#define TANH 298
#define Y0 299
#define Y1 300
#define RINT 301




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 385 "src/Parameters/parse.y"
{ 
  int logical_type;  
  int integer_type; 
  double scalar_type;  
  char * string_type; 
  char * subgroup_type;
  struct node_expr * node_type;
  }
/* Line 1489 of yacc.c.  */
#line 150 "src/Parameters/parse.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

