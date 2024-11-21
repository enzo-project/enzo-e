/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 1 "build/Cello/parse.y"

/* See LICENSE_CELLO file for license and copyright information */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#ifndef __APPLE__
	#include <malloc.h>
#endif
	
#define YYDEBUG 1

  /* Quiet a few -Wall errors */

int yylex (void);
void yyrestart  (FILE * input_file );
void yyerror(char *s);
void yylex_destroy();

#include "parse.tab.h"

#include "parse.h"

const char * node_name[] = {
  "node_unknown",
  "node_operation",
  "node_float",
  "node_integer",
  "node_variable",
  "node_function"
  };


const char * op_name[] = {
    "+",
    "-",
    "*",
    "/",
    "^",
    "<=",
    "<",
    ">=",
    ">",
    "==",
    "!=",
    "&&",
    "||"};

  /* ANY CHANGES HERE MUST BE REFLECTED IN parse.h enum_parameter[] */
  const char * parameter_name[]  = {
    "unknown",
    "sentinel",
    "group",
    "integer",
    "float",
    "string",
    "identifier",
    "logical",
    "list",
    "float_expr",
    "logical_expr",
    "function" };

  /* Structure for storing a single parameter / value pair in a linked list */


  struct node_expr * new_node_operation
    (struct node_expr * left, 
     enum enum_op oper,
     struct node_expr * right)
  {
    
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_operation;
    node->op_value      = oper;
    node->left          = left;
    node->right         = right;
    node->function_name = NULL;
    return node;
  }

  struct node_expr * new_node_float (double value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_float;
    node->float_value  = value;
    node->left          = NULL;
    node->right         = NULL;
    node->function_name = NULL;
    return node;
  }
  struct node_expr * new_node_logical (int value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_integer;
    node->integer_value = value;
    node->left          = NULL;
    node->right         = NULL;
    node->function_name = NULL;
    return node;
  }
  struct node_expr * new_node_variable (char * value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_variable;
    node->var_value     = value[0];
    node->left          = NULL;
    node->right         = NULL;
    node->function_name = NULL;
    free (value);
    return node;
  }
  struct node_expr * new_node_function
    (double (*function)(double),
     char * function_name,
     struct node_expr * argument)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_function;
    node->fun_value     = function;
    node->left          = argument;
    node->right         = NULL;
    node->function_name = strdup(function_name);
    return node;
  }

  /* The head of the linked list of parameter / value pairs */

  struct param_struct * param_head = NULL; /* head of entire list */
  struct param_struct * param_curr = NULL; /* head of current list */

  /* The current groups and parameter type */


  char *              current_parameter = NULL;
  char *              current_group[MAX_GROUP_DEPTH];
  int                 current_group_level = 0;
  enum enum_parameter current_type      = enum_parameter_sentinel;

  void clear_groups (char * groups[]) {
    int i;
    for (i=0; i<MAX_GROUP_DEPTH; i++) {
      groups[i] = 0; 
    }
  }

  void copy_groups (char * group_dest[], char * group_src[]) {
    int i;
    for (i=0; i<MAX_GROUP_DEPTH; i++) {
      /* MEMORY LEAK */
      group_dest[i] = (group_src[i]) ? strdup(group_src[i]) : 0;
    }
  }

  /* Insert a parameter into the list */

  void insert_param(struct param_struct * head, struct param_struct * new)
  {
     new->next  = head->next;
     head->next = new;
  }

  /* Delete a parameter from the list given a pointer to the previous element */

  void delete_param(struct param_struct * previous)
  {
    struct param_struct * item = previous->next;
    previous->next = item->next;
    free (item);     
  }

  struct param_struct * reverse_param(struct param_struct * old_head)
  {
    /* Keep sentinel the same */

    struct param_struct * new_head = old_head;

    struct param_struct * p = old_head;
    struct param_struct * c = p->next;
    struct param_struct * n = c->next;

    do {
      /* If parameter is a list, recursively reverse it as well */
      if (c->type == enum_parameter_list) {
	c->list_value = reverse_param(c->list_value);
      }
      c->next = p;
      p       = c;
      c       = n;
      n = n->next;
    } while (p->type != enum_parameter_sentinel) ;

    new_head = p;
    return new_head;
  }

  /* Function for creating and inserting a new parameter / value pair */
  /* in the linked list */

  struct param_struct * new_param ()
  {
    /* Create the new node */

    /* MEMORY LEAK */
     struct param_struct * p = 
       (struct param_struct *) malloc (sizeof (struct param_struct));

   /* Fill in the non-type-specific values for the new node */

     /* MEMORY LEAK */
     
     copy_groups(p->group,current_group);

     p->parameter = (current_parameter) ? strdup(current_parameter) : 0;

     free (current_parameter);
     current_parameter = 0;

     current_type = enum_parameter_unknown;

     insert_param(param_curr,p);

   /* Clear variables for the next assignment */

     return p;
  }

  /* Find the parameter if it exists, for use in e.g. List append */


  struct param_struct * find_param ()
  {
    int i;
    struct param_struct * p = param_head;
    struct param_struct * q;

    int match;
    do {
      match = 1;
      for (i=0; i<MAX_GROUP_DEPTH && p->group[i] && current_group[i]; i++) {
	if (strcmp(p->group[i],current_group[i]) != 0) match=0;
      }
      if (p->parameter==NULL) 
	match=0;
      else {
	if (strcmp(p->parameter,current_parameter) != 0) match=0;
      }
      q=p;
      p = p -> next;
    } while (!match && p->type != enum_parameter_sentinel);

    return match ? q : NULL;
  }


  /* New integer parameter assignment */

  void new_param_integer (int value)
  {
    struct param_struct * p = new_param();
    p->type          = enum_parameter_integer;
    p->integer_value = value;
  }


  /* New floating-point parameter assignment */

  void new_param_float (double value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_float;
    p->float_value = value;
  }

  /* New logical parameter assignment */

  void new_param_logical (int value)
  {
    struct param_struct * p = new_param();
    p->type          = enum_parameter_logical;
    p->logical_value = value;
  }

  /* New string parameter assignment */
  void new_param_string (char * value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_string;
    p->string_value = value;
  }

  /* New subgroup  */
  void new_param_group (char * value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_group;
    p->string_value = value;
  }

  /* New empty parameter assignment: FIRST NODE IN LIST IS A SENTINEL  */
  struct param_struct * new_param_sentinel ()
  {
    /* MEMORY LEAK */
    struct param_struct * p = 
      (struct param_struct *) malloc (sizeof (struct param_struct));

    clear_groups(p->group);
    p->parameter = NULL;
    p->type      = enum_parameter_sentinel;
    p->next       = p;
    p->list_value = NULL;

    return p;
  }

  /* New list parameter assignment */

  void new_param_list (struct param_struct * curr)
  {
    struct param_struct * p = new_param();
    p->type       = enum_parameter_list;
    p->list_value = curr;
  }

  /* New string parameter assignment */

  void new_param_expr (enum enum_parameter type,
		       struct node_expr * value)
  {
    struct param_struct * p = new_param();
    p->type     = type;
    p->op_value = value;
  }

  void new_parameter()
  {
     switch (current_type) {
     case enum_parameter_group:
       new_param_group(yylval.group_type);
       break;
     case enum_parameter_integer:
       new_param_integer(yylval.integer_type);
       break;
     case enum_parameter_float:
       new_param_float(yylval.float_type);
       break;
     case enum_parameter_string: 
       new_param_string(yylval.string_type);
       break;
     case enum_parameter_logical:
       new_param_logical(yylval.logical_type);
       break;
     case enum_parameter_list:
       break;
     case enum_parameter_float_expr:
       new_param_expr(enum_parameter_float_expr,yylval.node_type);
       break;
     case enum_parameter_logical_expr:
       new_param_expr(enum_parameter_logical_expr,yylval.node_type);
       break;
    default:
       printf ("%s:%d Parse Error: unknown type %d\n",
	       __FILE__,__LINE__,current_type);
       exit(1);
       break;
     }
  }


#line 450 "build/Cello/parse.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

#include "parse.tab.h"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_STRING = 3,                     /* STRING  */
  YYSYMBOL_IDENTIFIER = 4,                 /* IDENTIFIER  */
  YYSYMBOL_VARIABLE = 5,                   /* VARIABLE  */
  YYSYMBOL_FLOAT = 6,                      /* FLOAT  */
  YYSYMBOL_INTEGER = 7,                    /* INTEGER  */
  YYSYMBOL_LOGICAL = 8,                    /* LOGICAL  */
  YYSYMBOL_LE = 9,                         /* LE  */
  YYSYMBOL_GE = 10,                        /* GE  */
  YYSYMBOL_NE = 11,                        /* NE  */
  YYSYMBOL_EQ = 12,                        /* EQ  */
  YYSYMBOL_AND = 13,                       /* AND  */
  YYSYMBOL_OR = 14,                        /* OR  */
  YYSYMBOL_15_ = 15,                       /* '<'  */
  YYSYMBOL_16_ = 16,                       /* '>'  */
  YYSYMBOL_17_ = 17,                       /* '+'  */
  YYSYMBOL_18_ = 18,                       /* '-'  */
  YYSYMBOL_19_ = 19,                       /* '*'  */
  YYSYMBOL_20_ = 20,                       /* '/'  */
  YYSYMBOL_21_ = 21,                       /* '^'  */
  YYSYMBOL_ACOS = 22,                      /* ACOS  */
  YYSYMBOL_ACOSH = 23,                     /* ACOSH  */
  YYSYMBOL_APPEND = 24,                    /* APPEND  */
  YYSYMBOL_ASIN = 25,                      /* ASIN  */
  YYSYMBOL_ASINH = 26,                     /* ASINH  */
  YYSYMBOL_ATAN = 27,                      /* ATAN  */
  YYSYMBOL_ATANH = 28,                     /* ATANH  */
  YYSYMBOL_CBRT = 29,                      /* CBRT  */
  YYSYMBOL_CEIL = 30,                      /* CEIL  */
  YYSYMBOL_COS = 31,                       /* COS  */
  YYSYMBOL_COSH = 32,                      /* COSH  */
  YYSYMBOL_ERFC = 33,                      /* ERFC  */
  YYSYMBOL_ERF = 34,                       /* ERF  */
  YYSYMBOL_EXP = 35,                       /* EXP  */
  YYSYMBOL_EXPM1 = 36,                     /* EXPM1  */
  YYSYMBOL_FABS = 37,                      /* FABS  */
  YYSYMBOL_FLOOR = 38,                     /* FLOOR  */
  YYSYMBOL_J0 = 39,                        /* J0  */
  YYSYMBOL_J1 = 40,                        /* J1  */
  YYSYMBOL_LGAMMA = 41,                    /* LGAMMA  */
  YYSYMBOL_LOG10 = 42,                     /* LOG10  */
  YYSYMBOL_LOG1P = 43,                     /* LOG1P  */
  YYSYMBOL_LOGB = 44,                      /* LOGB  */
  YYSYMBOL_LOG = 45,                       /* LOG  */
  YYSYMBOL_PI = 46,                        /* PI  */
  YYSYMBOL_SIN = 47,                       /* SIN  */
  YYSYMBOL_SINH = 48,                      /* SINH  */
  YYSYMBOL_SQRT = 49,                      /* SQRT  */
  YYSYMBOL_TAN = 50,                       /* TAN  */
  YYSYMBOL_TANH = 51,                      /* TANH  */
  YYSYMBOL_Y0 = 52,                        /* Y0  */
  YYSYMBOL_Y1 = 53,                        /* Y1  */
  YYSYMBOL_RINT = 54,                      /* RINT  */
  YYSYMBOL_55_ = 55,                       /* '{'  */
  YYSYMBOL_56_ = 56,                       /* '}'  */
  YYSYMBOL_57_ = 57,                       /* ';'  */
  YYSYMBOL_58_ = 58,                       /* '='  */
  YYSYMBOL_59_ = 59,                       /* '['  */
  YYSYMBOL_60_ = 60,                       /* ']'  */
  YYSYMBOL_61_ = 61,                       /* ','  */
  YYSYMBOL_62_ = 62,                       /* '('  */
  YYSYMBOL_63_ = 63,                       /* ')'  */
  YYSYMBOL_YYACCEPT = 64,                  /* $accept  */
  YYSYMBOL_file = 65,                      /* file  */
  YYSYMBOL_group = 66,                     /* group  */
  YYSYMBOL_parameter_group = 67,           /* parameter_group  */
  YYSYMBOL_parameter_list = 68,            /* parameter_list  */
  YYSYMBOL_group_name = 69,                /* group_name  */
  YYSYMBOL_parameter_name = 70,            /* parameter_name  */
  YYSYMBOL_parameter_assignment = 71,      /* parameter_assignment  */
  YYSYMBOL_parameter_append = 72,          /* parameter_append  */
  YYSYMBOL_parameter_value = 73,           /* parameter_value  */
  YYSYMBOL_list = 74,                      /* list  */
  YYSYMBOL_existing_list = 75,             /* existing_list  */
  YYSYMBOL_LIST_BEGIN = 76,                /* LIST_BEGIN  */
  YYSYMBOL_LIST_APPEND = 77,               /* LIST_APPEND  */
  YYSYMBOL_LIST_END = 78,                  /* LIST_END  */
  YYSYMBOL_list_elements = 79,             /* list_elements  */
  YYSYMBOL_80_1 = 80,                      /* $@1  */
  YYSYMBOL_cle = 81,                       /* cle  */
  YYSYMBOL_cse = 82,                       /* cse  */
  YYSYMBOL_cie = 83,                       /* cie  */
  YYSYMBOL_vse = 84,                       /* vse  */
  YYSYMBOL_vle = 85                        /* vle  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1108

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  64
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  22
/* YYNRULES -- Number of rules.  */
#define YYNRULES  164
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  347

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   302


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      62,    63,    19,    17,    61,    18,     2,    20,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    57,
      15,    58,    16,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    59,     2,    60,    21,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    55,     2,    56,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    53,    54
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   465,   465,   466,   470,   473,   474,   477,   478,   479,
     480,   481,   482,   483,   484,   488,   493,   496,   500,   504,
     505,   506,   507,   508,   509,   510,   516,   517,   519,   520,
     523,   531,   551,   558,   559,   559,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   579,   580,   581,   582,
     583,   584,   585,   586,   587,   588,   589,   590,   591,   592,
     593,   594,   595,   596,   597,   598,   599,   600,   602,   603,
     604,   605,   606,   607,   608,   609,   610,   611,   612,   613,
     614,   615,   616,   617,   618,   622,   623,   624,   625,   626,
     627,   628,   632,   633,   634,   635,   636,   637,   638,   639,
     640,   641,   642,   643,   644,   645,   646,   647,   648,   649,
     650,   651,   652,   653,   654,   655,   656,   657,   658,   659,
     660,   661,   662,   663,   665,   666,   667,   668,   669,   670,
     671,   672,   673,   674,   675,   676,   677,   678,   679,   680,
     685,   686,   687,   688,   689,   690,   691,   692,   693,   694,
     695,   696,   697,   698,   699,   700,   701,   702,   703,   704,
     705,   706,   707,   708,   709
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "STRING", "IDENTIFIER",
  "VARIABLE", "FLOAT", "INTEGER", "LOGICAL", "LE", "GE", "NE", "EQ", "AND",
  "OR", "'<'", "'>'", "'+'", "'-'", "'*'", "'/'", "'^'", "ACOS", "ACOSH",
  "APPEND", "ASIN", "ASINH", "ATAN", "ATANH", "CBRT", "CEIL", "COS",
  "COSH", "ERFC", "ERF", "EXP", "EXPM1", "FABS", "FLOOR", "J0", "J1",
  "LGAMMA", "LOG10", "LOG1P", "LOGB", "LOG", "PI", "SIN", "SINH", "SQRT",
  "TAN", "TANH", "Y0", "Y1", "RINT", "'{'", "'}'", "';'", "'='", "'['",
  "']'", "','", "'('", "')'", "$accept", "file", "group",
  "parameter_group", "parameter_list", "group_name", "parameter_name",
  "parameter_assignment", "parameter_append", "parameter_value", "list",
  "existing_list", "LIST_BEGIN", "LIST_APPEND", "LIST_END",
  "list_elements", "$@1", "cle", "cse", "cie", "vse", "vle", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-66)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-16)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -66,    84,   -66,   -66,   -66,   -38,     4,   -66,    31,    71,
      81,    63,   -66,   -66,   -66,   -66,    93,    87,   -66,   -66,
      97,   311,   -66,   -66,   -66,    96,   253,   -66,   -66,   -66,
     -66,   -66,   108,   110,   120,   121,   154,   157,   159,   161,
     175,   176,   177,   178,   179,   193,   205,   206,   207,   208,
     209,   210,   211,   212,   215,   -66,   247,   248,   249,   258,
     259,   260,   261,   262,   -66,   361,   101,   -66,   253,   -11,
      -5,   122,    10,   109,   -66,   -66,   -66,   -66,   -37,   461,
     461,   461,   461,   461,   461,   461,   461,   461,   461,   461,
     461,   461,   461,   461,   461,   461,   461,   461,   461,   461,
     461,   461,   461,   461,   461,   461,   461,   461,   461,   461,
      89,   115,    94,   459,    92,   -66,   -66,   -37,   411,   411,
     461,   461,   461,   461,   461,   461,   461,   461,   461,   461,
     461,     2,     2,     2,     2,     2,   461,   461,   461,   461,
     461,   461,   461,   461,   461,   461,   461,   411,   411,   311,
     -66,   461,   245,   357,   409,   499,   507,   512,   517,   522,
     527,   532,   537,   546,   584,   593,   598,   603,   608,   613,
     618,   623,   631,   670,   678,   683,   688,   693,   698,   703,
     708,   717,   755,   764,   769,   774,   779,   784,   789,   794,
     802,   841,   849,   854,   859,   864,   869,   874,   879,   888,
     926,   935,   940,   945,   950,   955,   960,   965,   973,  1012,
    1020,  1025,  1030,  1035,   -66,   -66,   -66,   -66,   -66,   -66,
     411,   -66,    -5,    10,   -66,   138,   164,   143,   148,   143,
     148,   143,   148,   143,   148,   143,   148,   143,   148,    88,
     155,    88,   155,   196,   201,   196,   201,   196,   201,     2,
     160,   160,   304,   304,   304,   143,   148,   143,   148,   143,
     148,   143,   148,   143,   148,   143,   148,    88,   155,    88,
     155,   196,   201,   196,   201,   196,   201,   -66,   -66,   138,
     164,   -66,  1040,  1045,   -66,   -66,   -66,   -66,   -66,   -66,
     -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,
     -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,
     -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,
     -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,
     -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,   -66,
     -66,   -66,   -66,   -66,   -66,   -66,   -66
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,    15,     3,     0,     0,     4,    16,     9,
       0,     0,     7,     8,    10,     5,     0,    13,    11,    12,
       0,     0,     6,    14,    31,     0,     0,    19,   139,    83,
      91,    45,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    84,     0,     0,     0,     0,
       0,     0,     0,     0,    30,     0,     0,    25,     0,    22,
      21,    20,    23,    24,    18,    32,    33,    29,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    17,    27,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      28,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    36,    46,    85,    92,   140,    26,
       0,    44,     0,     0,   163,    43,   160,    37,   142,    38,
     145,    42,   157,    41,   154,    39,   148,    40,   151,    47,
      94,    48,    97,    49,   100,    50,   103,    51,   106,     0,
      86,    87,    88,    89,    90,   141,   143,   144,   146,   156,
     158,   153,   155,   147,   149,   150,   152,    93,    95,    96,
      98,    99,   101,   102,   104,   105,   107,   162,   164,   159,
     161,    34,     0,     0,    52,   108,    53,   109,    54,   110,
      55,   111,    56,   112,    57,   113,    58,   114,    59,   115,
      60,   116,    61,   117,    62,   118,    63,   119,    64,   120,
      65,   121,    66,   122,    67,   123,    68,   124,    69,   125,
      70,   126,    71,   127,    72,   128,    73,   129,    74,   130,
      75,   131,    76,   132,    77,   133,    78,   134,    79,   135,
      80,   136,    81,   137,    82,   138,    35
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -66,   -66,   144,   -66,   -66,   -66,   -66,   317,   318,   -20,
     -66,   -66,   -66,   -66,    42,   263,   -66,    -2,   -47,   -65,
     106,     0
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,     1,     4,     7,    10,     5,    11,    12,    13,    76,
      67,    25,    68,    26,    77,    78,   346,    69,    70,    71,
      72,    73
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     112,    66,   118,   119,   120,   121,   122,   123,     8,    30,
     124,   125,   126,   127,   128,   129,   130,     6,   111,   136,
     137,   138,   139,    75,   149,   140,   141,   142,   143,   144,
     145,   146,   152,   154,   156,   158,   160,   162,   164,   166,
     168,   170,   172,   174,   176,   178,   180,   182,   184,   186,
     188,   190,   192,   194,   196,   198,   200,   202,   204,   206,
     208,   210,   212,   110,   249,   114,   250,   251,   252,   253,
     254,   222,   222,   227,   229,   231,   233,   235,   237,   239,
     241,   243,   245,   247,     2,     8,   -15,    20,     3,   255,
     257,   259,   261,   263,   265,   267,   269,   271,   273,   275,
     222,   222,   118,   119,   282,   147,   148,   128,   129,   130,
     116,   131,   132,   133,   134,   135,   221,   225,   224,   226,
     150,    21,   147,   148,   120,   121,   122,   123,    14,   281,
     124,   125,   126,   127,   128,   129,   130,    15,    16,   131,
     132,   133,   134,   135,    23,   277,   279,   278,   280,    22,
       9,   118,   214,    74,    17,   218,    24,   216,   115,   219,
     126,   127,   128,   129,   130,   142,   143,   144,   145,   146,
      79,   113,    80,   111,   144,   145,   146,   147,   215,   133,
     134,   135,    81,    82,   112,   153,   155,   157,   159,   161,
     163,   165,   167,   169,   171,   173,   175,   177,   179,   181,
     183,   185,   187,   189,   191,   193,   195,   197,   199,   201,
     203,   205,   207,   209,   211,   213,    83,   130,   110,    84,
     114,    85,   146,    86,   223,   223,   228,   230,   232,   234,
     236,   238,   240,   242,   244,   246,   248,    87,    88,    89,
      90,    91,   256,   258,   260,   262,   264,   266,   268,   270,
     272,   274,   276,   223,   223,    92,    27,   283,    28,    29,
      30,    31,   126,   127,   128,   129,   130,    93,    94,    95,
      96,    97,    98,    99,   100,    32,    33,   101,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,    60,    61,    62,    63,   284,   102,
     103,   104,    64,    75,    27,    65,    28,    29,    30,    31,
     105,   106,   107,   108,   109,   135,   113,    18,    19,     0,
       0,   117,     0,    32,    33,     0,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    28,    29,    30,    31,
      64,     0,     0,    65,   142,   143,   144,   145,   146,     0,
       0,     0,     0,    32,    33,     0,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    28,    29,     0,    31,
     285,     0,     0,    65,     0,     0,   126,   127,   128,   129,
     130,     0,     0,    32,    33,     0,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    28,    29,   136,   137,
     138,   139,   286,   220,   140,   141,   142,   143,   144,   145,
     146,     0,     0,    32,    33,     0,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,   142,   143,   144,   145,
     146,     0,   217,   151,   126,   127,   128,   129,   130,   142,
     143,   144,   145,   146,   126,   127,   128,   129,   130,   142,
     143,   144,   145,   146,   126,   127,   128,   129,   130,   142,
     143,   144,   145,   146,   126,   127,   128,   129,   130,     0,
       0,     0,   287,   142,   143,   144,   145,   146,     0,     0,
     288,     0,     0,     0,     0,   289,     0,     0,     0,     0,
     290,     0,     0,     0,     0,   291,     0,     0,     0,     0,
     292,     0,     0,     0,     0,   293,     0,     0,     0,     0,
     294,   126,   127,   128,   129,   130,     0,     0,     0,   295,
     142,   143,   144,   145,   146,   126,   127,   128,   129,   130,
     142,   143,   144,   145,   146,   126,   127,   128,   129,   130,
     142,   143,   144,   145,   146,   126,   127,   128,   129,   130,
     142,   143,   144,   145,   146,     0,     0,   296,   126,   127,
     128,   129,   130,     0,     0,     0,   297,     0,     0,     0,
       0,   298,     0,     0,     0,     0,   299,     0,     0,     0,
       0,   300,     0,     0,     0,     0,   301,     0,     0,     0,
       0,   302,     0,     0,     0,     0,   303,   142,   143,   144,
     145,   146,     0,     0,   304,   126,   127,   128,   129,   130,
     142,   143,   144,   145,   146,   126,   127,   128,   129,   130,
     142,   143,   144,   145,   146,   126,   127,   128,   129,   130,
     142,   143,   144,   145,   146,   126,   127,   128,   129,   130,
       0,     0,     0,   305,   142,   143,   144,   145,   146,     0,
       0,   306,     0,     0,     0,     0,   307,     0,     0,     0,
       0,   308,     0,     0,     0,     0,   309,     0,     0,     0,
       0,   310,     0,     0,     0,     0,   311,     0,     0,     0,
       0,   312,   126,   127,   128,   129,   130,     0,     0,     0,
     313,   142,   143,   144,   145,   146,   126,   127,   128,   129,
     130,   142,   143,   144,   145,   146,   126,   127,   128,   129,
     130,   142,   143,   144,   145,   146,   126,   127,   128,   129,
     130,   142,   143,   144,   145,   146,     0,     0,   314,   126,
     127,   128,   129,   130,     0,     0,     0,   315,     0,     0,
       0,     0,   316,     0,     0,     0,     0,   317,     0,     0,
       0,     0,   318,     0,     0,     0,     0,   319,     0,     0,
       0,     0,   320,     0,     0,     0,     0,   321,   142,   143,
     144,   145,   146,     0,     0,   322,   126,   127,   128,   129,
     130,   142,   143,   144,   145,   146,   126,   127,   128,   129,
     130,   142,   143,   144,   145,   146,   126,   127,   128,   129,
     130,   142,   143,   144,   145,   146,   126,   127,   128,   129,
     130,     0,     0,     0,   323,   142,   143,   144,   145,   146,
       0,     0,   324,     0,     0,     0,     0,   325,     0,     0,
       0,     0,   326,     0,     0,     0,     0,   327,     0,     0,
       0,     0,   328,     0,     0,     0,     0,   329,     0,     0,
       0,     0,   330,   126,   127,   128,   129,   130,     0,     0,
       0,   331,   142,   143,   144,   145,   146,   126,   127,   128,
     129,   130,   142,   143,   144,   145,   146,   126,   127,   128,
     129,   130,   142,   143,   144,   145,   146,   126,   127,   128,
     129,   130,   142,   143,   144,   145,   146,     0,     0,   332,
     126,   127,   128,   129,   130,     0,     0,     0,   333,     0,
       0,     0,     0,   334,     0,     0,     0,     0,   335,     0,
       0,     0,     0,   336,     0,     0,     0,     0,   337,     0,
       0,     0,     0,   338,     0,     0,     0,     0,   339,   142,
     143,   144,   145,   146,     0,     0,   340,   126,   127,   128,
     129,   130,   142,   143,   144,   145,   146,   126,   127,   128,
     129,   130,   142,   143,   144,   145,   146,   126,   127,   128,
     129,   130,   142,   143,   144,   145,   146,     0,     0,     0,
       0,     0,     0,     0,     0,   341,     0,     0,     0,     0,
       0,     0,     0,   342,     0,     0,     0,     0,   343,     0,
       0,     0,     0,   344,     0,     0,     0,     0,   345,     0,
       0,     0,     0,   215,     0,     0,     0,     0,   217
};

static const yytype_int16 yycheck[] =
{
      65,    21,    13,    14,     9,    10,    11,    12,     4,     7,
      15,    16,    17,    18,    19,    20,    21,    55,    65,     9,
      10,    11,    12,    60,    61,    15,    16,    17,    18,    19,
      20,    21,    79,    80,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   109,    65,    62,    65,   131,   132,   133,   134,
     135,   118,   119,   120,   121,   122,   123,   124,   125,   126,
     127,   128,   129,   130,     0,     4,    55,    24,     4,   136,
     137,   138,   139,   140,   141,   142,   143,   144,   145,   146,
     147,   148,    13,    14,   151,    13,    14,    19,    20,    21,
      68,    17,    18,    19,    20,    21,   118,   119,   118,   119,
      78,    58,    13,    14,     9,    10,    11,    12,    57,   149,
      15,    16,    17,    18,    19,    20,    21,    56,    57,    17,
      18,    19,    20,    21,    57,   147,   148,   147,   148,    56,
       6,    13,    63,    57,    10,    63,    59,    63,    57,   117,
      17,    18,    19,    20,    21,    17,    18,    19,    20,    21,
      62,    65,    62,   220,    19,    20,    21,    13,    63,    19,
      20,    21,    62,    62,   249,    79,    80,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,    62,    21,   220,    62,
     220,    62,    21,    62,   118,   119,   120,   121,   122,   123,
     124,   125,   126,   127,   128,   129,   130,    62,    62,    62,
      62,    62,   136,   137,   138,   139,   140,   141,   142,   143,
     144,   145,   146,   147,   148,    62,     3,   151,     5,     6,
       7,     8,    17,    18,    19,    20,    21,    62,    62,    62,
      62,    62,    62,    62,    62,    22,    23,    62,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    63,    62,
      62,    62,    59,    60,     3,    62,     5,     6,     7,     8,
      62,    62,    62,    62,    62,    21,   220,    10,    10,    -1,
      -1,    68,    -1,    22,    23,    -1,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,     5,     6,     7,     8,
      59,    -1,    -1,    62,    17,    18,    19,    20,    21,    -1,
      -1,    -1,    -1,    22,    23,    -1,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,     5,     6,    -1,     8,
      63,    -1,    -1,    62,    -1,    -1,    17,    18,    19,    20,
      21,    -1,    -1,    22,    23,    -1,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,     5,     6,     9,    10,
      11,    12,    63,    62,    15,    16,    17,    18,    19,    20,
      21,    -1,    -1,    22,    23,    -1,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    17,    18,    19,    20,
      21,    -1,    63,    62,    17,    18,    19,    20,    21,    17,
      18,    19,    20,    21,    17,    18,    19,    20,    21,    17,
      18,    19,    20,    21,    17,    18,    19,    20,    21,    17,
      18,    19,    20,    21,    17,    18,    19,    20,    21,    -1,
      -1,    -1,    63,    17,    18,    19,    20,    21,    -1,    -1,
      63,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,
      63,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,
      63,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,
      63,    17,    18,    19,    20,    21,    -1,    -1,    -1,    63,
      17,    18,    19,    20,    21,    17,    18,    19,    20,    21,
      17,    18,    19,    20,    21,    17,    18,    19,    20,    21,
      17,    18,    19,    20,    21,    17,    18,    19,    20,    21,
      17,    18,    19,    20,    21,    -1,    -1,    63,    17,    18,
      19,    20,    21,    -1,    -1,    -1,    63,    -1,    -1,    -1,
      -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,
      -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,
      -1,    63,    -1,    -1,    -1,    -1,    63,    17,    18,    19,
      20,    21,    -1,    -1,    63,    17,    18,    19,    20,    21,
      17,    18,    19,    20,    21,    17,    18,    19,    20,    21,
      17,    18,    19,    20,    21,    17,    18,    19,    20,    21,
      17,    18,    19,    20,    21,    17,    18,    19,    20,    21,
      -1,    -1,    -1,    63,    17,    18,    19,    20,    21,    -1,
      -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,
      -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,
      -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,
      -1,    63,    17,    18,    19,    20,    21,    -1,    -1,    -1,
      63,    17,    18,    19,    20,    21,    17,    18,    19,    20,
      21,    17,    18,    19,    20,    21,    17,    18,    19,    20,
      21,    17,    18,    19,    20,    21,    17,    18,    19,    20,
      21,    17,    18,    19,    20,    21,    -1,    -1,    63,    17,
      18,    19,    20,    21,    -1,    -1,    -1,    63,    -1,    -1,
      -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,
      -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,
      -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    17,    18,
      19,    20,    21,    -1,    -1,    63,    17,    18,    19,    20,
      21,    17,    18,    19,    20,    21,    17,    18,    19,    20,
      21,    17,    18,    19,    20,    21,    17,    18,    19,    20,
      21,    17,    18,    19,    20,    21,    17,    18,    19,    20,
      21,    -1,    -1,    -1,    63,    17,    18,    19,    20,    21,
      -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,
      -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,
      -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,    -1,
      -1,    -1,    63,    17,    18,    19,    20,    21,    -1,    -1,
      -1,    63,    17,    18,    19,    20,    21,    17,    18,    19,
      20,    21,    17,    18,    19,    20,    21,    17,    18,    19,
      20,    21,    17,    18,    19,    20,    21,    17,    18,    19,
      20,    21,    17,    18,    19,    20,    21,    -1,    -1,    63,
      17,    18,    19,    20,    21,    -1,    -1,    -1,    63,    -1,
      -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,
      -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,
      -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    17,
      18,    19,    20,    21,    -1,    -1,    63,    17,    18,    19,
      20,    21,    17,    18,    19,    20,    21,    17,    18,    19,
      20,    21,    17,    18,    19,    20,    21,    17,    18,    19,
      20,    21,    17,    18,    19,    20,    21,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,
      -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    63,    -1,
      -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    63
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,    65,     0,     4,    66,    69,    55,    67,     4,    66,
      68,    70,    71,    72,    57,    56,    57,    66,    71,    72,
      24,    58,    56,    57,    59,    75,    77,     3,     5,     6,
       7,     8,    22,    23,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,    52,    53,    54,    59,    62,    73,    74,    76,    81,
      82,    83,    84,    85,    57,    60,    73,    78,    79,    62,
      62,    62,    62,    62,    62,    62,    62,    62,    62,    62,
      62,    62,    62,    62,    62,    62,    62,    62,    62,    62,
      62,    62,    62,    62,    62,    62,    62,    62,    62,    62,
      81,    82,    83,    84,    85,    57,    78,    79,    13,    14,
       9,    10,    11,    12,    15,    16,    17,    18,    19,    20,
      21,    17,    18,    19,    20,    21,     9,    10,    11,    12,
      15,    16,    17,    18,    19,    20,    21,    13,    14,    61,
      78,    62,    82,    84,    82,    84,    82,    84,    82,    84,
      82,    84,    82,    84,    82,    84,    82,    84,    82,    84,
      82,    84,    82,    84,    82,    84,    82,    84,    82,    84,
      82,    84,    82,    84,    82,    84,    82,    84,    82,    84,
      82,    84,    82,    84,    82,    84,    82,    84,    82,    84,
      82,    84,    82,    84,    82,    84,    82,    84,    82,    84,
      82,    84,    82,    84,    63,    63,    63,    63,    63,    78,
      62,    81,    82,    84,    85,    81,    85,    82,    84,    82,
      84,    82,    84,    82,    84,    82,    84,    82,    84,    82,
      84,    82,    84,    82,    84,    82,    84,    82,    84,    62,
      83,    83,    83,    83,    83,    82,    84,    82,    84,    82,
      84,    82,    84,    82,    84,    82,    84,    82,    84,    82,
      84,    82,    84,    82,    84,    82,    84,    81,    85,    81,
      85,    73,    82,    84,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    80
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    64,    65,    65,    66,    67,    67,    68,    68,    68,
      68,    68,    68,    68,    68,    69,    70,    71,    72,    73,
      73,    73,    73,    73,    73,    73,    74,    74,    75,    75,
      76,    77,    78,    79,    80,    79,    81,    81,    81,    81,
      81,    81,    81,    81,    81,    81,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    83,    83,    83,    83,    83,
      83,    83,    84,    84,    84,    84,    84,    84,    84,    84,
      84,    84,    84,    84,    84,    84,    84,    84,    84,    84,
      84,    84,    84,    84,    84,    84,    84,    84,    84,    84,
      84,    84,    84,    84,    84,    84,    84,    84,    84,    84,
      84,    84,    84,    84,    84,    84,    84,    84,    84,    84,
      85,    85,    85,    85,    85,    85,    85,    85,    85,    85,
      85,    85,    85,    85,    85,    85,    85,    85,    85,    85,
      85,    85,    85,    85,    85
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     0,     2,     2,     3,     4,     1,     1,     1,
       2,     2,     2,     2,     3,     1,     1,     4,     4,     1,
       1,     1,     1,     1,     1,     1,     3,     2,     3,     2,
       1,     1,     1,     1,     0,     4,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     1,     3,     3,     3,     3,
       3,     3,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     1,     1,     3,     3,     3,     3,     3,
       3,     1,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YY_USE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 3: /* file: file group  */
#line 466 "build/Cello/parse.y"
              { }
#line 1897 "build/Cello/parse.tab.c"
    break;

  case 4: /* group: group_name parameter_group  */
#line 470 "build/Cello/parse.y"
                            {  }
#line 1903 "build/Cello/parse.tab.c"
    break;

  case 5: /* parameter_group: '{' parameter_list '}'  */
#line 473 "build/Cello/parse.y"
                                 { current_group[--current_group_level] = 0; }
#line 1909 "build/Cello/parse.tab.c"
    break;

  case 6: /* parameter_group: '{' parameter_list ';' '}'  */
#line 474 "build/Cello/parse.y"
                                 { current_group[--current_group_level] = 0; }
#line 1915 "build/Cello/parse.tab.c"
    break;

  case 7: /* parameter_list: parameter_assignment  */
#line 477 "build/Cello/parse.y"
                                         {  }
#line 1921 "build/Cello/parse.tab.c"
    break;

  case 8: /* parameter_list: parameter_append  */
#line 478 "build/Cello/parse.y"
                                     {  }
#line 1927 "build/Cello/parse.tab.c"
    break;

  case 9: /* parameter_list: group  */
#line 479 "build/Cello/parse.y"
                          {  }
#line 1933 "build/Cello/parse.tab.c"
    break;

  case 10: /* parameter_list: group ';'  */
#line 480 "build/Cello/parse.y"
                              {  }
#line 1939 "build/Cello/parse.tab.c"
    break;

  case 11: /* parameter_list: parameter_list parameter_assignment  */
#line 481 "build/Cello/parse.y"
                                         {  }
#line 1945 "build/Cello/parse.tab.c"
    break;

  case 12: /* parameter_list: parameter_list parameter_append  */
#line 482 "build/Cello/parse.y"
                                     {  }
#line 1951 "build/Cello/parse.tab.c"
    break;

  case 13: /* parameter_list: parameter_list group  */
#line 483 "build/Cello/parse.y"
                          {  }
#line 1957 "build/Cello/parse.tab.c"
    break;

  case 14: /* parameter_list: parameter_list group ';'  */
#line 484 "build/Cello/parse.y"
                              {  }
#line 1963 "build/Cello/parse.tab.c"
    break;

  case 15: /* group_name: IDENTIFIER  */
#line 488 "build/Cello/parse.y"
             {
  current_group[current_group_level++] = (yyvsp[0].string_type); 
}
#line 1971 "build/Cello/parse.tab.c"
    break;

  case 16: /* parameter_name: IDENTIFIER  */
#line 493 "build/Cello/parse.y"
             { current_parameter = (yyvsp[0].string_type);}
#line 1977 "build/Cello/parse.tab.c"
    break;

  case 17: /* parameter_assignment: parameter_name '=' parameter_value ';'  */
#line 496 "build/Cello/parse.y"
                                       { new_parameter(); }
#line 1983 "build/Cello/parse.tab.c"
    break;

  case 18: /* parameter_append: parameter_name APPEND existing_list ';'  */
#line 500 "build/Cello/parse.y"
                                        { new_parameter(); }
#line 1989 "build/Cello/parse.tab.c"
    break;

  case 19: /* parameter_value: STRING  */
#line 504 "build/Cello/parse.y"
       { current_type = enum_parameter_string;       yylval.string_type = (yyvsp[0].string_type); }
#line 1995 "build/Cello/parse.tab.c"
    break;

  case 20: /* parameter_value: cie  */
#line 505 "build/Cello/parse.y"
       { current_type = enum_parameter_integer;      yylval.integer_type = (yyvsp[0].integer_type);}
#line 2001 "build/Cello/parse.tab.c"
    break;

  case 21: /* parameter_value: cse  */
#line 506 "build/Cello/parse.y"
       { current_type = enum_parameter_float;        yylval.float_type = (yyvsp[0].float_type);}
#line 2007 "build/Cello/parse.tab.c"
    break;

  case 22: /* parameter_value: cle  */
#line 507 "build/Cello/parse.y"
       { current_type = enum_parameter_logical;      yylval.logical_type = (yyvsp[0].logical_type); }
#line 2013 "build/Cello/parse.tab.c"
    break;

  case 23: /* parameter_value: vse  */
#line 508 "build/Cello/parse.y"
       { current_type = enum_parameter_float_expr;   yylval.node_type = (yyvsp[0].node_type); }
#line 2019 "build/Cello/parse.tab.c"
    break;

  case 24: /* parameter_value: vle  */
#line 509 "build/Cello/parse.y"
       { current_type = enum_parameter_logical_expr; yylval.node_type = (yyvsp[0].node_type); }
#line 2025 "build/Cello/parse.tab.c"
    break;

  case 25: /* parameter_value: list  */
#line 510 "build/Cello/parse.y"
       { current_type = enum_parameter_list; }
#line 2031 "build/Cello/parse.tab.c"
    break;

  case 26: /* list: LIST_BEGIN list_elements LIST_END  */
#line 516 "build/Cello/parse.y"
                                        {  }
#line 2037 "build/Cello/parse.tab.c"
    break;

  case 27: /* list: LIST_BEGIN LIST_END  */
#line 517 "build/Cello/parse.y"
                          {  }
#line 2043 "build/Cello/parse.tab.c"
    break;

  case 28: /* existing_list: LIST_APPEND list_elements LIST_END  */
#line 519 "build/Cello/parse.y"
                                                  {  }
#line 2049 "build/Cello/parse.tab.c"
    break;

  case 29: /* existing_list: LIST_APPEND LIST_END  */
#line 520 "build/Cello/parse.y"
                           {  }
#line 2055 "build/Cello/parse.tab.c"
    break;

  case 30: /* LIST_BEGIN: '['  */
#line 523 "build/Cello/parse.y"
     { 
   struct param_struct * p = new_param_sentinel();
   p->list_value = param_curr; /* save param_curr */
   new_param_list(p);
   param_curr = p;
 }
#line 2066 "build/Cello/parse.tab.c"
    break;

  case 31: /* LIST_APPEND: '['  */
#line 531 "build/Cello/parse.y"
     { 
  
   struct param_struct * p = find_param();

   if (p) {
     /* s is sentinal in list corresponding to p in LIST_BEGIN */
     struct param_struct * s = p->list_value;
     s->list_value=param_curr;
     param_curr = s;
     current_type = enum_parameter_list;
   } else {
     /* new list */
     p = new_param_sentinel();
     p->list_value = param_curr; /* save param_curr */
     new_param_list(p);
     param_curr = p;
   }
 }
#line 2089 "build/Cello/parse.tab.c"
    break;

  case 32: /* LIST_END: ']'  */
#line 551 "build/Cello/parse.y"
    {
  current_type = enum_parameter_list;
  param_curr = param_curr->list_value; /* restore param_curr */ 
}
#line 2098 "build/Cello/parse.tab.c"
    break;

  case 33: /* list_elements: parameter_value  */
#line 558 "build/Cello/parse.y"
                                         { new_parameter(); }
#line 2104 "build/Cello/parse.tab.c"
    break;

  case 34: /* $@1: %empty  */
#line 559 "build/Cello/parse.y"
                                         { new_parameter(); }
#line 2110 "build/Cello/parse.tab.c"
    break;

  case 35: /* list_elements: list_elements ',' parameter_value $@1  */
#line 561 "build/Cello/parse.y"
{ }
#line 2116 "build/Cello/parse.tab.c"
    break;

  case 36: /* cle: '(' cle ')'  */
#line 566 "build/Cello/parse.y"
            { (yyval.logical_type) = (yyvsp[-1].logical_type); }
#line 2122 "build/Cello/parse.tab.c"
    break;

  case 37: /* cle: cse LE cse  */
#line 567 "build/Cello/parse.y"
               { (yyval.logical_type) = (yyvsp[-2].float_type) <= (yyvsp[0].float_type); }
#line 2128 "build/Cello/parse.tab.c"
    break;

  case 38: /* cle: cse GE cse  */
#line 568 "build/Cello/parse.y"
               { (yyval.logical_type) = (yyvsp[-2].float_type) >= (yyvsp[0].float_type); }
#line 2134 "build/Cello/parse.tab.c"
    break;

  case 39: /* cle: cse '<' cse  */
#line 569 "build/Cello/parse.y"
               { (yyval.logical_type) = (yyvsp[-2].float_type) <  (yyvsp[0].float_type); }
#line 2140 "build/Cello/parse.tab.c"
    break;

  case 40: /* cle: cse '>' cse  */
#line 570 "build/Cello/parse.y"
               { (yyval.logical_type) = (yyvsp[-2].float_type) >  (yyvsp[0].float_type); }
#line 2146 "build/Cello/parse.tab.c"
    break;

  case 41: /* cle: cse EQ cse  */
#line 571 "build/Cello/parse.y"
               { (yyval.logical_type) = (yyvsp[-2].float_type) == (yyvsp[0].float_type); }
#line 2152 "build/Cello/parse.tab.c"
    break;

  case 42: /* cle: cse NE cse  */
#line 572 "build/Cello/parse.y"
               { (yyval.logical_type) = (yyvsp[-2].float_type) != (yyvsp[0].float_type); }
#line 2158 "build/Cello/parse.tab.c"
    break;

  case 43: /* cle: cle OR cle  */
#line 573 "build/Cello/parse.y"
               { (yyval.logical_type) = (yyvsp[-2].logical_type) || (yyvsp[0].logical_type); }
#line 2164 "build/Cello/parse.tab.c"
    break;

  case 44: /* cle: cle AND cle  */
#line 574 "build/Cello/parse.y"
               { (yyval.logical_type) = (yyvsp[-2].logical_type) && (yyvsp[0].logical_type); }
#line 2170 "build/Cello/parse.tab.c"
    break;

  case 45: /* cle: LOGICAL  */
#line 575 "build/Cello/parse.y"
           { (yyval.logical_type) = (yyvsp[0].logical_type); }
#line 2176 "build/Cello/parse.tab.c"
    break;

  case 46: /* cse: '(' cse ')'  */
#line 579 "build/Cello/parse.y"
             { (yyval.float_type) = (yyvsp[-1].float_type); }
#line 2182 "build/Cello/parse.tab.c"
    break;

  case 47: /* cse: cse '+' cse  */
#line 580 "build/Cello/parse.y"
               { (yyval.float_type) = (yyvsp[-2].float_type) + (yyvsp[0].float_type);}
#line 2188 "build/Cello/parse.tab.c"
    break;

  case 48: /* cse: cse '-' cse  */
#line 581 "build/Cello/parse.y"
               { (yyval.float_type) = (yyvsp[-2].float_type) - (yyvsp[0].float_type);}
#line 2194 "build/Cello/parse.tab.c"
    break;

  case 49: /* cse: cse '*' cse  */
#line 582 "build/Cello/parse.y"
               { (yyval.float_type) = (yyvsp[-2].float_type) * (yyvsp[0].float_type);}
#line 2200 "build/Cello/parse.tab.c"
    break;

  case 50: /* cse: cse '/' cse  */
#line 583 "build/Cello/parse.y"
               { (yyval.float_type) = (yyvsp[-2].float_type) / (yyvsp[0].float_type);}
#line 2206 "build/Cello/parse.tab.c"
    break;

  case 51: /* cse: cse '^' cse  */
#line 584 "build/Cello/parse.y"
               { (yyval.float_type) = pow((double)(yyvsp[-2].float_type), (double)(yyvsp[0].float_type)); }
#line 2212 "build/Cello/parse.tab.c"
    break;

  case 52: /* cse: ACOS '(' cse ')'  */
#line 585 "build/Cello/parse.y"
                   { (yyval.float_type) = acos((yyvsp[-1].float_type)); }
#line 2218 "build/Cello/parse.tab.c"
    break;

  case 53: /* cse: ACOSH '(' cse ')'  */
#line 586 "build/Cello/parse.y"
                    { (yyval.float_type) = acosh((yyvsp[-1].float_type)); }
#line 2224 "build/Cello/parse.tab.c"
    break;

  case 54: /* cse: ASIN '(' cse ')'  */
#line 587 "build/Cello/parse.y"
                   { (yyval.float_type) = asin((yyvsp[-1].float_type)); }
#line 2230 "build/Cello/parse.tab.c"
    break;

  case 55: /* cse: ASINH '(' cse ')'  */
#line 588 "build/Cello/parse.y"
                    { (yyval.float_type) = asinh((yyvsp[-1].float_type)); }
#line 2236 "build/Cello/parse.tab.c"
    break;

  case 56: /* cse: ATAN '(' cse ')'  */
#line 589 "build/Cello/parse.y"
                   { (yyval.float_type) = atan((yyvsp[-1].float_type)); }
#line 2242 "build/Cello/parse.tab.c"
    break;

  case 57: /* cse: ATANH '(' cse ')'  */
#line 590 "build/Cello/parse.y"
                    { (yyval.float_type) = atanh((yyvsp[-1].float_type)); }
#line 2248 "build/Cello/parse.tab.c"
    break;

  case 58: /* cse: CBRT '(' cse ')'  */
#line 591 "build/Cello/parse.y"
                   { (yyval.float_type) = cbrt((yyvsp[-1].float_type)); }
#line 2254 "build/Cello/parse.tab.c"
    break;

  case 59: /* cse: CEIL '(' cse ')'  */
#line 592 "build/Cello/parse.y"
                   { (yyval.float_type) = ceil((yyvsp[-1].float_type)); }
#line 2260 "build/Cello/parse.tab.c"
    break;

  case 60: /* cse: COS '(' cse ')'  */
#line 593 "build/Cello/parse.y"
                  { (yyval.float_type) = cos((yyvsp[-1].float_type)); }
#line 2266 "build/Cello/parse.tab.c"
    break;

  case 61: /* cse: COSH '(' cse ')'  */
#line 594 "build/Cello/parse.y"
                   { (yyval.float_type) = cosh((yyvsp[-1].float_type)); }
#line 2272 "build/Cello/parse.tab.c"
    break;

  case 62: /* cse: ERFC '(' cse ')'  */
#line 595 "build/Cello/parse.y"
                   { (yyval.float_type) = erfc((yyvsp[-1].float_type)); }
#line 2278 "build/Cello/parse.tab.c"
    break;

  case 63: /* cse: ERF '(' cse ')'  */
#line 596 "build/Cello/parse.y"
                  { (yyval.float_type) = erf((yyvsp[-1].float_type)); }
#line 2284 "build/Cello/parse.tab.c"
    break;

  case 64: /* cse: EXP '(' cse ')'  */
#line 597 "build/Cello/parse.y"
                  { (yyval.float_type) = exp((yyvsp[-1].float_type)); }
#line 2290 "build/Cello/parse.tab.c"
    break;

  case 65: /* cse: EXPM1 '(' cse ')'  */
#line 598 "build/Cello/parse.y"
                    { (yyval.float_type) = expm1((yyvsp[-1].float_type)); }
#line 2296 "build/Cello/parse.tab.c"
    break;

  case 66: /* cse: FABS '(' cse ')'  */
#line 599 "build/Cello/parse.y"
                   { (yyval.float_type) = fabs((yyvsp[-1].float_type)); }
#line 2302 "build/Cello/parse.tab.c"
    break;

  case 67: /* cse: FLOOR '(' cse ')'  */
#line 600 "build/Cello/parse.y"
                    { (yyval.float_type) = floor((yyvsp[-1].float_type)); }
#line 2308 "build/Cello/parse.tab.c"
    break;

  case 68: /* cse: J0 '(' cse ')'  */
#line 602 "build/Cello/parse.y"
                 { (yyval.float_type) = j0((yyvsp[-1].float_type)); }
#line 2314 "build/Cello/parse.tab.c"
    break;

  case 69: /* cse: J1 '(' cse ')'  */
#line 603 "build/Cello/parse.y"
                 { (yyval.float_type) = j1((yyvsp[-1].float_type)); }
#line 2320 "build/Cello/parse.tab.c"
    break;

  case 70: /* cse: LGAMMA '(' cse ')'  */
#line 604 "build/Cello/parse.y"
                     { (yyval.float_type) = lgamma((yyvsp[-1].float_type)); }
#line 2326 "build/Cello/parse.tab.c"
    break;

  case 71: /* cse: LOG10 '(' cse ')'  */
#line 605 "build/Cello/parse.y"
                    { (yyval.float_type) = log10((yyvsp[-1].float_type)); }
#line 2332 "build/Cello/parse.tab.c"
    break;

  case 72: /* cse: LOG1P '(' cse ')'  */
#line 606 "build/Cello/parse.y"
                    { (yyval.float_type) = log1p((yyvsp[-1].float_type)); }
#line 2338 "build/Cello/parse.tab.c"
    break;

  case 73: /* cse: LOGB '(' cse ')'  */
#line 607 "build/Cello/parse.y"
                   { (yyval.float_type) = logb((yyvsp[-1].float_type)); }
#line 2344 "build/Cello/parse.tab.c"
    break;

  case 74: /* cse: LOG '(' cse ')'  */
#line 608 "build/Cello/parse.y"
                  { (yyval.float_type) = log((yyvsp[-1].float_type)); }
#line 2350 "build/Cello/parse.tab.c"
    break;

  case 75: /* cse: SIN '(' cse ')'  */
#line 609 "build/Cello/parse.y"
                  { (yyval.float_type) = sin((yyvsp[-1].float_type)); }
#line 2356 "build/Cello/parse.tab.c"
    break;

  case 76: /* cse: SINH '(' cse ')'  */
#line 610 "build/Cello/parse.y"
                   { (yyval.float_type) = sinh((yyvsp[-1].float_type)); }
#line 2362 "build/Cello/parse.tab.c"
    break;

  case 77: /* cse: SQRT '(' cse ')'  */
#line 611 "build/Cello/parse.y"
                   { (yyval.float_type) = sqrt((yyvsp[-1].float_type)); }
#line 2368 "build/Cello/parse.tab.c"
    break;

  case 78: /* cse: TAN '(' cse ')'  */
#line 612 "build/Cello/parse.y"
                  { (yyval.float_type) = tan((yyvsp[-1].float_type)); }
#line 2374 "build/Cello/parse.tab.c"
    break;

  case 79: /* cse: TANH '(' cse ')'  */
#line 613 "build/Cello/parse.y"
                   { (yyval.float_type) = tanh((yyvsp[-1].float_type)); }
#line 2380 "build/Cello/parse.tab.c"
    break;

  case 80: /* cse: Y0 '(' cse ')'  */
#line 614 "build/Cello/parse.y"
                 { (yyval.float_type) = y0((yyvsp[-1].float_type)); }
#line 2386 "build/Cello/parse.tab.c"
    break;

  case 81: /* cse: Y1 '(' cse ')'  */
#line 615 "build/Cello/parse.y"
                 { (yyval.float_type) = y1((yyvsp[-1].float_type)); }
#line 2392 "build/Cello/parse.tab.c"
    break;

  case 82: /* cse: RINT '(' cse ')'  */
#line 616 "build/Cello/parse.y"
                   { (yyval.float_type) = rint((yyvsp[-1].float_type)); }
#line 2398 "build/Cello/parse.tab.c"
    break;

  case 83: /* cse: FLOAT  */
#line 617 "build/Cello/parse.y"
        { (yyval.float_type) = (yyvsp[0].float_type);}
#line 2404 "build/Cello/parse.tab.c"
    break;

  case 84: /* cse: PI  */
#line 618 "build/Cello/parse.y"
     { (yyval.float_type) = M_PI ; }
#line 2410 "build/Cello/parse.tab.c"
    break;

  case 85: /* cie: '(' cie ')'  */
#line 622 "build/Cello/parse.y"
             { (yyval.integer_type) = (yyvsp[-1].integer_type); }
#line 2416 "build/Cello/parse.tab.c"
    break;

  case 86: /* cie: cie '+' cie  */
#line 623 "build/Cello/parse.y"
               { (yyval.integer_type) = (yyvsp[-2].integer_type) + (yyvsp[0].integer_type);}
#line 2422 "build/Cello/parse.tab.c"
    break;

  case 87: /* cie: cie '-' cie  */
#line 624 "build/Cello/parse.y"
               { (yyval.integer_type) = (yyvsp[-2].integer_type) - (yyvsp[0].integer_type);}
#line 2428 "build/Cello/parse.tab.c"
    break;

  case 88: /* cie: cie '*' cie  */
#line 625 "build/Cello/parse.y"
               { (yyval.integer_type) = (yyvsp[-2].integer_type) * (yyvsp[0].integer_type);}
#line 2434 "build/Cello/parse.tab.c"
    break;

  case 89: /* cie: cie '/' cie  */
#line 626 "build/Cello/parse.y"
               { (yyval.integer_type) = (yyvsp[-2].integer_type) / (yyvsp[0].integer_type);}
#line 2440 "build/Cello/parse.tab.c"
    break;

  case 90: /* cie: cie '^' cie  */
#line 627 "build/Cello/parse.y"
               { (yyval.integer_type) = pow((double)(yyvsp[-2].integer_type), (double)(yyvsp[0].integer_type));}
#line 2446 "build/Cello/parse.tab.c"
    break;

  case 91: /* cie: INTEGER  */
#line 628 "build/Cello/parse.y"
           { (yyval.integer_type) = (yyvsp[0].integer_type);}
#line 2452 "build/Cello/parse.tab.c"
    break;

  case 92: /* vse: '(' vse ')'  */
#line 632 "build/Cello/parse.y"
               { (yyval.node_type) = (yyvsp[-1].node_type); }
#line 2458 "build/Cello/parse.tab.c"
    break;

  case 93: /* vse: vse '+' cse  */
#line 633 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_add,new_node_float((yyvsp[0].float_type))); }
#line 2464 "build/Cello/parse.tab.c"
    break;

  case 94: /* vse: cse '+' vse  */
#line 634 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_add,(yyvsp[0].node_type)); }
#line 2470 "build/Cello/parse.tab.c"
    break;

  case 95: /* vse: vse '+' vse  */
#line 635 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_add,(yyvsp[0].node_type)); }
#line 2476 "build/Cello/parse.tab.c"
    break;

  case 96: /* vse: vse '-' cse  */
#line 636 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_sub,new_node_float((yyvsp[0].float_type))); }
#line 2482 "build/Cello/parse.tab.c"
    break;

  case 97: /* vse: cse '-' vse  */
#line 637 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_sub,(yyvsp[0].node_type)); }
#line 2488 "build/Cello/parse.tab.c"
    break;

  case 98: /* vse: vse '-' vse  */
#line 638 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_sub,(yyvsp[0].node_type)); }
#line 2494 "build/Cello/parse.tab.c"
    break;

  case 99: /* vse: vse '*' cse  */
#line 639 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_mul,new_node_float((yyvsp[0].float_type))); }
#line 2500 "build/Cello/parse.tab.c"
    break;

  case 100: /* vse: cse '*' vse  */
#line 640 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_mul,(yyvsp[0].node_type)); }
#line 2506 "build/Cello/parse.tab.c"
    break;

  case 101: /* vse: vse '*' vse  */
#line 641 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_mul,(yyvsp[0].node_type)); }
#line 2512 "build/Cello/parse.tab.c"
    break;

  case 102: /* vse: vse '/' cse  */
#line 642 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_div,new_node_float((yyvsp[0].float_type))); }
#line 2518 "build/Cello/parse.tab.c"
    break;

  case 103: /* vse: cse '/' vse  */
#line 643 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_div,(yyvsp[0].node_type)); }
#line 2524 "build/Cello/parse.tab.c"
    break;

  case 104: /* vse: vse '/' vse  */
#line 644 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_div,(yyvsp[0].node_type)); }
#line 2530 "build/Cello/parse.tab.c"
    break;

  case 105: /* vse: vse '^' cse  */
#line 645 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_pow, new_node_float((yyvsp[0].float_type))); }
#line 2536 "build/Cello/parse.tab.c"
    break;

  case 106: /* vse: cse '^' vse  */
#line 646 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_pow, (yyvsp[0].node_type)); }
#line 2542 "build/Cello/parse.tab.c"
    break;

  case 107: /* vse: vse '^' vse  */
#line 647 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_pow, (yyvsp[0].node_type)); }
#line 2548 "build/Cello/parse.tab.c"
    break;

  case 108: /* vse: ACOS '(' vse ')'  */
#line 648 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( acos, "acos", (yyvsp[-1].node_type)); }
#line 2554 "build/Cello/parse.tab.c"
    break;

  case 109: /* vse: ACOSH '(' vse ')'  */
#line 649 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( acosh, "acosh", (yyvsp[-1].node_type)); }
#line 2560 "build/Cello/parse.tab.c"
    break;

  case 110: /* vse: ASIN '(' vse ')'  */
#line 650 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( asin, "asin", (yyvsp[-1].node_type)); }
#line 2566 "build/Cello/parse.tab.c"
    break;

  case 111: /* vse: ASINH '(' vse ')'  */
#line 651 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( asinh, "asinh", (yyvsp[-1].node_type)); }
#line 2572 "build/Cello/parse.tab.c"
    break;

  case 112: /* vse: ATAN '(' vse ')'  */
#line 652 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( atan, "atan", (yyvsp[-1].node_type)); }
#line 2578 "build/Cello/parse.tab.c"
    break;

  case 113: /* vse: ATANH '(' vse ')'  */
#line 653 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( atanh, "atanh", (yyvsp[-1].node_type)); }
#line 2584 "build/Cello/parse.tab.c"
    break;

  case 114: /* vse: CBRT '(' vse ')'  */
#line 654 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( cbrt, "cbrt", (yyvsp[-1].node_type)); }
#line 2590 "build/Cello/parse.tab.c"
    break;

  case 115: /* vse: CEIL '(' vse ')'  */
#line 655 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( ceil, "ceil", (yyvsp[-1].node_type)); }
#line 2596 "build/Cello/parse.tab.c"
    break;

  case 116: /* vse: COS '(' vse ')'  */
#line 656 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( cos, "cos", (yyvsp[-1].node_type)); }
#line 2602 "build/Cello/parse.tab.c"
    break;

  case 117: /* vse: COSH '(' vse ')'  */
#line 657 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( cosh, "cosh", (yyvsp[-1].node_type)); }
#line 2608 "build/Cello/parse.tab.c"
    break;

  case 118: /* vse: ERFC '(' vse ')'  */
#line 658 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( erfc, "erfc", (yyvsp[-1].node_type)); }
#line 2614 "build/Cello/parse.tab.c"
    break;

  case 119: /* vse: ERF '(' vse ')'  */
#line 659 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( erf, "erf", (yyvsp[-1].node_type)); }
#line 2620 "build/Cello/parse.tab.c"
    break;

  case 120: /* vse: EXP '(' vse ')'  */
#line 660 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( exp, "exp", (yyvsp[-1].node_type)); }
#line 2626 "build/Cello/parse.tab.c"
    break;

  case 121: /* vse: EXPM1 '(' vse ')'  */
#line 661 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( expm1, "expm1", (yyvsp[-1].node_type)); }
#line 2632 "build/Cello/parse.tab.c"
    break;

  case 122: /* vse: FABS '(' vse ')'  */
#line 662 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( fabs, "fabs", (yyvsp[-1].node_type)); }
#line 2638 "build/Cello/parse.tab.c"
    break;

  case 123: /* vse: FLOOR '(' vse ')'  */
#line 663 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( floor, "floor", (yyvsp[-1].node_type)); }
#line 2644 "build/Cello/parse.tab.c"
    break;

  case 124: /* vse: J0 '(' vse ')'  */
#line 665 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( j0, "j0", (yyvsp[-1].node_type)); }
#line 2650 "build/Cello/parse.tab.c"
    break;

  case 125: /* vse: J1 '(' vse ')'  */
#line 666 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( j1, "j1", (yyvsp[-1].node_type)); }
#line 2656 "build/Cello/parse.tab.c"
    break;

  case 126: /* vse: LGAMMA '(' vse ')'  */
#line 667 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( lgamma, "lgamma", (yyvsp[-1].node_type)); }
#line 2662 "build/Cello/parse.tab.c"
    break;

  case 127: /* vse: LOG10 '(' vse ')'  */
#line 668 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( log10, "log10", (yyvsp[-1].node_type)); }
#line 2668 "build/Cello/parse.tab.c"
    break;

  case 128: /* vse: LOG1P '(' vse ')'  */
#line 669 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( log1p, "log1p", (yyvsp[-1].node_type)); }
#line 2674 "build/Cello/parse.tab.c"
    break;

  case 129: /* vse: LOGB '(' vse ')'  */
#line 670 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( logb, "logb", (yyvsp[-1].node_type)); }
#line 2680 "build/Cello/parse.tab.c"
    break;

  case 130: /* vse: LOG '(' vse ')'  */
#line 671 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( log, "log", (yyvsp[-1].node_type)); }
#line 2686 "build/Cello/parse.tab.c"
    break;

  case 131: /* vse: SIN '(' vse ')'  */
#line 672 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( sin, "sin", (yyvsp[-1].node_type)); }
#line 2692 "build/Cello/parse.tab.c"
    break;

  case 132: /* vse: SINH '(' vse ')'  */
#line 673 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( sinh, "sinh", (yyvsp[-1].node_type)); }
#line 2698 "build/Cello/parse.tab.c"
    break;

  case 133: /* vse: SQRT '(' vse ')'  */
#line 674 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( sqrt, "sqrt", (yyvsp[-1].node_type)); }
#line 2704 "build/Cello/parse.tab.c"
    break;

  case 134: /* vse: TAN '(' vse ')'  */
#line 675 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( tan, "tan", (yyvsp[-1].node_type)); }
#line 2710 "build/Cello/parse.tab.c"
    break;

  case 135: /* vse: TANH '(' vse ')'  */
#line 676 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( tanh, "tanh", (yyvsp[-1].node_type)); }
#line 2716 "build/Cello/parse.tab.c"
    break;

  case 136: /* vse: Y0 '(' vse ')'  */
#line 677 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( y0, "y0", (yyvsp[-1].node_type)); }
#line 2722 "build/Cello/parse.tab.c"
    break;

  case 137: /* vse: Y1 '(' vse ')'  */
#line 678 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( y1, "y1", (yyvsp[-1].node_type)); }
#line 2728 "build/Cello/parse.tab.c"
    break;

  case 138: /* vse: RINT '(' vse ')'  */
#line 679 "build/Cello/parse.y"
                      { (yyval.node_type) = new_node_function ( rint, "rint", (yyvsp[-1].node_type)); }
#line 2734 "build/Cello/parse.tab.c"
    break;

  case 139: /* vse: VARIABLE  */
#line 680 "build/Cello/parse.y"
            { (yyval.node_type) = new_node_variable ((yyvsp[0].string_type));  }
#line 2740 "build/Cello/parse.tab.c"
    break;

  case 140: /* vle: '(' vle ')'  */
#line 685 "build/Cello/parse.y"
             { (yyval.node_type) = (yyvsp[-1].node_type); }
#line 2746 "build/Cello/parse.tab.c"
    break;

  case 141: /* vle: vse LE cse  */
#line 686 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_le,new_node_float((yyvsp[0].float_type))); }
#line 2752 "build/Cello/parse.tab.c"
    break;

  case 142: /* vle: cse LE vse  */
#line 687 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_le,(yyvsp[0].node_type)); }
#line 2758 "build/Cello/parse.tab.c"
    break;

  case 143: /* vle: vse LE vse  */
#line 688 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_le,(yyvsp[0].node_type)); }
#line 2764 "build/Cello/parse.tab.c"
    break;

  case 144: /* vle: vse GE cse  */
#line 689 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_ge,new_node_float((yyvsp[0].float_type))); }
#line 2770 "build/Cello/parse.tab.c"
    break;

  case 145: /* vle: cse GE vse  */
#line 690 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_ge,(yyvsp[0].node_type)); }
#line 2776 "build/Cello/parse.tab.c"
    break;

  case 146: /* vle: vse GE vse  */
#line 691 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_ge,(yyvsp[0].node_type)); }
#line 2782 "build/Cello/parse.tab.c"
    break;

  case 147: /* vle: vse '<' cse  */
#line 692 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_lt,new_node_float((yyvsp[0].float_type))); }
#line 2788 "build/Cello/parse.tab.c"
    break;

  case 148: /* vle: cse '<' vse  */
#line 693 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_lt,(yyvsp[0].node_type)); }
#line 2794 "build/Cello/parse.tab.c"
    break;

  case 149: /* vle: vse '<' vse  */
#line 694 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_lt,(yyvsp[0].node_type)); }
#line 2800 "build/Cello/parse.tab.c"
    break;

  case 150: /* vle: vse '>' cse  */
#line 695 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_gt,new_node_float((yyvsp[0].float_type))); }
#line 2806 "build/Cello/parse.tab.c"
    break;

  case 151: /* vle: cse '>' vse  */
#line 696 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_gt,(yyvsp[0].node_type)); }
#line 2812 "build/Cello/parse.tab.c"
    break;

  case 152: /* vle: vse '>' vse  */
#line 697 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_gt,(yyvsp[0].node_type)); }
#line 2818 "build/Cello/parse.tab.c"
    break;

  case 153: /* vle: vse EQ cse  */
#line 698 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_eq,new_node_float((yyvsp[0].float_type))); }
#line 2824 "build/Cello/parse.tab.c"
    break;

  case 154: /* vle: cse EQ vse  */
#line 699 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_eq,(yyvsp[0].node_type)); }
#line 2830 "build/Cello/parse.tab.c"
    break;

  case 155: /* vle: vse EQ vse  */
#line 700 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_eq,(yyvsp[0].node_type)); }
#line 2836 "build/Cello/parse.tab.c"
    break;

  case 156: /* vle: vse NE cse  */
#line 701 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_ne,new_node_float((yyvsp[0].float_type))); }
#line 2842 "build/Cello/parse.tab.c"
    break;

  case 157: /* vle: cse NE vse  */
#line 702 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_ne,(yyvsp[0].node_type)); }
#line 2848 "build/Cello/parse.tab.c"
    break;

  case 158: /* vle: vse NE vse  */
#line 703 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_ne,(yyvsp[0].node_type)); }
#line 2854 "build/Cello/parse.tab.c"
    break;

  case 159: /* vle: vle OR cle  */
#line 704 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_or,new_node_logical((yyvsp[0].logical_type))); }
#line 2860 "build/Cello/parse.tab.c"
    break;

  case 160: /* vle: cle OR vle  */
#line 705 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_logical((yyvsp[-2].logical_type)), enum_op_or,(yyvsp[0].node_type)); }
#line 2866 "build/Cello/parse.tab.c"
    break;

  case 161: /* vle: vle OR vle  */
#line 706 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_or,(yyvsp[0].node_type)); }
#line 2872 "build/Cello/parse.tab.c"
    break;

  case 162: /* vle: vle AND cle  */
#line 707 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_and,new_node_logical((yyvsp[0].logical_type))); }
#line 2878 "build/Cello/parse.tab.c"
    break;

  case 163: /* vle: cle AND vle  */
#line 708 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation (new_node_logical((yyvsp[-2].logical_type)), enum_op_and,(yyvsp[0].node_type)); }
#line 2884 "build/Cello/parse.tab.c"
    break;

  case 164: /* vle: vle AND vle  */
#line 709 "build/Cello/parse.y"
               { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_and,(yyvsp[0].node_type)); }
#line 2890 "build/Cello/parse.tab.c"
    break;


#line 2894 "build/Cello/parse.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 713 "build/Cello/parse.y"


int cello_new_file(const char * filename);

struct param_struct * 
cello_parameters_read(const char * filename, FILE * fp)
{
  clear_groups(current_group);

  /* initialize the linked list with an initial sentinel (sentinel) node */
  param_head = param_curr = new_param_sentinel();

  /*   yydebug=1; */

  cello_new_file (filename);

  yyrestart(fp);

  yyparse();
  yylex_destroy();

  param_head = reverse_param(param_head);
  return param_head;
}

void indent (int level)
{
  int i;
  for (i=0; i<level; i++) {
    printf ("  "); 
  }
}

void print_expression (struct node_expr * node,
		       FILE * fp)
{
  if (node == NULL) {
    fprintf (fp,"NULL");
  } else {
    char left,right;
    switch (node->type) {
    case enum_node_integer:
      fprintf (fp,"%d",node->integer_value);
      break;
    case enum_node_float:
      /* '#' format character forces a decimal point */
      fprintf (fp,FLOAT_FORMAT,node->float_value);
      break;
    case enum_node_variable:
      fprintf (fp,"%c",node->var_value);
      break;
    case enum_node_function:
      fprintf (fp,"%s(",node->function_name);
      print_expression(node->left,fp);
      fprintf (fp,")");
      break;
    case enum_node_operation:
      left  = (node->left->type == enum_node_operation) ? '(' : ' ';
      right = (node->left->type == enum_node_operation) ? ')' : ' ';
      fprintf (fp,"%c",left);
      print_expression(node->left,fp);
      fprintf (fp,"%c",right);
      fprintf (fp," %s ",op_name[node->op_value]);
      left  = (node->right->type == enum_node_operation) ? '(' : ' ';
      right = (node->right->type == enum_node_operation) ? ')' : ' ';
      fprintf (fp,"%c",left);
      print_expression(node->right,fp);
      fprintf (fp,"%c",right);
      break;
    default:
      break;
    }
    fflush(fp);
  }

}

struct crude_str_buffer_{
  char* start;
  char* stop; // first address that isn't in the buffer
};

void print_to_buf_(struct crude_str_buffer_ * buffer, const char * fmt, ...) {
  if ((buffer->start + 1) >= buffer->stop) {
    return;
  }
  va_list args;
  va_start(args, fmt);
  size_t size = buffer->stop - buffer->start;
  int count = vsnprintf(buffer->start, size, fmt, args);
  va_end(args);

  if (count < 0) {
    buffer->start = buffer->stop;
  } else {
    // count returns the total number of counters that would be written, without
    // truncation, while not including the terminating null-byte
    buffer->start += count;
  }
}

void sprintf_expression_helper_(struct node_expr * node,
                                struct crude_str_buffer_ * buffer)
{
  if (node == NULL) {
    print_to_buf_(buffer,"NULL");
  } else {
    char left,right;
    switch (node->type) {
    case enum_node_integer:
      print_to_buf_ (buffer,"%d",node->integer_value);
      break;
    case enum_node_float:
      /* '#' format character forces a decimal point */
      print_to_buf_ (buffer,FLOAT_FORMAT,node->float_value);
      break;
    case enum_node_variable:
      print_to_buf_ (buffer,"%c",node->var_value);
      break;
    case enum_node_function:
      print_to_buf_ (buffer,"%s(",node->function_name);
      sprintf_expression_helper_(node->left,buffer);
      print_to_buf_ (buffer,")");
      break;
    case enum_node_operation:
      left  = (node->left->type == enum_node_operation) ? '(' : ' ';
      right = (node->left->type == enum_node_operation) ? ')' : ' ';
      print_to_buf_ (buffer,"%c",left);
      sprintf_expression_helper_(node->left,buffer);
      print_to_buf_ (buffer,"%c",right);
      print_to_buf_ (buffer," %s ",op_name[node->op_value]);
      left  = (node->right->type == enum_node_operation) ? '(' : ' ';
      right = (node->right->type == enum_node_operation) ? ')' : ' ';
      print_to_buf_ (buffer,"%c",left);
      sprintf_expression_helper_(node->right,buffer);
      print_to_buf_ (buffer,"%c",right);
      break;
    default:
      break;
    }
  }
}

void sprintf_expression (struct node_expr * node,
                         char * buffer, size_t buffer_size)
{
  struct crude_str_buffer_ buf = {buffer, buffer+buffer_size};
  sprintf_expression_helper_(node,&buf);
}


void cello_parameters_print_list(struct param_struct * head, int level);
void cello_print_parameter(struct param_struct * p, int level)
{
  int i;
  if (p->group != NULL) {
    indent(level);
    if (parameter_name[p->type])
    printf ("%s ", parameter_name[p->type]);
    for (i=0; p->group[i] != NULL && i < MAX_GROUP_DEPTH; i++) {
      printf ("%s:",p->group[i]);
    }
    printf ("[%s] = \n", p->parameter);
  } else {
    /* list element */
    indent(level);
    printf ("%s %s = ", 
	    parameter_name[p->type], p->parameter);
  }
  switch (p->type) {
  case enum_parameter_float:  
    /* '#' format character forces a decimal point */
    printf (FLOAT_FORMAT,p->float_value);  
    break;
  case enum_parameter_integer: 
    printf ("%d\n",p->integer_value); 
    break;
  case enum_parameter_string:  
    printf ("%s\n",p->string_value); 
    break;
  case enum_parameter_group:  
    printf ("Uh oh: GROUP %s (should be deleted)\n",p->string_value);
    break;
  case enum_parameter_logical:
    printf ("%s\n",p->logical_value ? "true" : "false");
    break;
  case enum_parameter_list:    
    indent(level);
    printf ("[\n"); 
    cello_parameters_print_list(p->list_value, level + 1);
    indent(level);
    printf ("]\n"); 
    break;
  case enum_parameter_logical_expr:
    indent(level);
    print_expression(p->op_value,stdout); printf ("\n");
    break;
  case enum_parameter_float_expr:
    indent(level);
    print_expression(p->op_value,stdout); printf ("\n");
    break;
  default: 
    indent(level);
    printf ("unknown type\n"); 
    break;
  }
}
void cello_parameters_print_list(struct param_struct * head, int level)
{
  struct param_struct * p = head->next;

  while (p && p->type != enum_parameter_sentinel) {

    cello_print_parameter(p,level);

    p = p->next;
  }
}

void cello_parameters_print()
{
  cello_parameters_print_list(param_head,0);
}

