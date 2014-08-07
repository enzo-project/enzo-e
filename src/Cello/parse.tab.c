/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "build/Cello/parse.y" /* yacc.c:339  */

/* See LICENSE_CELLO file for license and copyright information */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

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
  };


  void copy_groups (char * group_dest[], char * group_src[]) {
    int i;
    for (i=0; i<MAX_GROUP_DEPTH; i++) {
      /* MEMORY LEAK */
      group_dest[i] = (group_src[i]) ? strdup(group_src[i]) : 0;
    }
  };

  /* Function to update parameter's groups once the group is known */

/*   void update_group (char * group) */
/*     { */
/*       struct param_struct * p = param_curr; */
/*       while (p->next->type  != enum_parameter_sentinel &&  */
/* 	     p->next->group == NULL) { */
/* 	p->next->group = strdup(group); */
/*         p = p -> next; */
/*       } */
/*     } */

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

#line 429 "build/Cello/parse.tab.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parse.tab.h".  */
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
#line 365 "build/Cello/parse.y" /* yacc.c:355  */
 
  int logical_type;  
  int integer_type; 
  double float_type;  
  char * string_type; 
  char * group_type;
  struct node_expr * node_type;
  

#line 525 "build/Cello/parse.tab.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_BUILD_CELLO_PARSE_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 540 "build/Cello/parse.tab.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

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

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
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


#if ! defined yyoverflow || YYERROR_VERBOSE

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
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
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
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
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
#define YYLAST   969

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  61
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  19
/* YYNRULES -- Number of rules.  */
#define YYNRULES  150
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  325

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   300

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      59,    60,    19,    17,    58,    18,     2,    20,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    54,
      15,    55,    16,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    56,     2,    57,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    52,     2,    53,     2,     2,     2,     2,
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
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   443,   443,   444,   448,   451,   452,   455,   456,   457,
     458,   462,   465,   468,   472,   473,   474,   475,   476,   477,
     478,   481,   482,   485,   492,   496,   497,   497,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   517,   518,
     519,   520,   521,   522,   523,   524,   525,   526,   527,   528,
     529,   530,   531,   532,   533,   534,   535,   536,   537,   539,
     540,   541,   542,   543,   544,   545,   546,   547,   548,   549,
     550,   551,   552,   553,   554,   558,   559,   560,   561,   562,
     563,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,   578,   579,   580,   581,   582,   583,   584,   585,
     586,   587,   588,   589,   590,   591,   592,   593,   594,   595,
     597,   598,   599,   600,   601,   602,   603,   604,   605,   606,
     607,   608,   609,   610,   611,   612,   617,   618,   619,   620,
     621,   622,   623,   624,   625,   626,   627,   628,   629,   630,
     631,   632,   633,   634,   635,   636,   637,   638,   639,   640,
     641
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "STRING", "IDENTIFIER", "VARIABLE",
  "FLOAT", "INTEGER", "LOGICAL", "LE", "GE", "NE", "EQ", "AND", "OR",
  "'<'", "'>'", "'+'", "'-'", "'*'", "'/'", "ACOS", "ACOSH", "ASIN",
  "ASINH", "ATAN", "ATANH", "CBRT", "CEIL", "COS", "COSH", "ERFC", "ERF",
  "EXP", "EXPM1", "FABS", "FLOOR", "J0", "J1", "LGAMMA", "LOG10", "LOG1P",
  "LOGB", "LOG", "SIN", "SINH", "SQRT", "TAN", "TANH", "Y0", "Y1", "RINT",
  "'{'", "'}'", "';'", "'='", "'['", "']'", "','", "'('", "')'", "$accept",
  "file", "group", "parameter_group", "parameter_list", "group_name",
  "parameter_name", "parameter_assignment", "parameter_value", "list",
  "LIST_BEGIN", "LIST_END", "list_elements", "$@1", "cle", "cse", "cie",
  "vse", "vle", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,    60,    62,    43,    45,    42,
      47,   270,   271,   272,   273,   274,   275,   276,   277,   278,
     279,   280,   281,   282,   283,   284,   285,   286,   287,   288,
     289,   290,   291,   292,   293,   294,   295,   296,   297,   298,
     299,   300,   123,   125,    59,    61,    91,    93,    44,    40,
      41
};
# endif

#define YYPACT_NINF -41

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-41)))

#define YYTABLE_NINF -13

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -41,    20,   -41,   -41,   -41,   -40,    28,   -41,   -27,   -41,
     -32,   -19,   -41,   -41,    29,   309,   -41,   -41,   -41,   -41,
     -41,   -41,   -41,   -41,   -20,   -10,    26,    64,    66,    67,
      73,   113,   114,   115,   126,   128,   141,   152,   153,   162,
     163,   174,   175,   176,   177,   178,   179,   184,   187,   188,
     191,   192,   193,   194,   196,   -41,   356,   -41,   -41,   251,
      13,    -2,    85,   921,    17,   450,   450,   450,   450,   450,
     450,   450,   450,   450,   450,   450,   450,   450,   450,   450,
     450,   450,   450,   450,   450,   450,   450,   450,   450,   450,
     450,   450,   450,   450,   450,   450,   -12,   208,   -14,   448,
      24,   -41,   -41,   -41,   -23,   403,   403,   450,   450,   450,
     450,   450,   450,   450,   450,   450,   450,    22,    22,    22,
      22,   450,   450,   450,   450,   450,   450,   450,   450,   450,
     450,   403,   403,   450,   111,   150,   185,   189,   243,   306,
     352,   401,   485,   493,   497,   501,   505,   509,   513,   517,
     521,   529,   565,   573,   577,   581,   585,   589,   593,   597,
     601,   609,   645,   653,   657,   661,   665,   669,   673,   677,
     681,   689,   725,   733,   737,   741,   745,   749,   753,   757,
     761,   769,   805,   813,   817,   821,   825,   829,   833,   837,
     841,   849,   885,   893,   897,   901,   -41,   -41,   -41,   -41,
     -41,   309,   -41,   403,   -41,    -2,   921,   -41,    27,    32,
     102,   212,   102,   212,   102,   212,   102,   212,   102,   212,
     102,   212,    23,    71,    23,    71,   -41,   -41,   -41,   -41,
      22,    97,    97,   -41,   -41,   102,   212,   102,   212,   102,
     212,   102,   212,   102,   212,   102,   212,    23,    71,    23,
      71,   -41,   -41,   -41,   -41,   -41,   -41,    27,    32,   905,
     909,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,
     -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,
     -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,
     -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,
     -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,
     -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,   -41,
     -41,   -41,   -41,   -41,   -41
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,    11,     3,     0,     0,     4,    11,     9,
       0,     0,     7,     5,     0,     0,     6,    10,     8,    14,
     125,    74,    80,    37,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    23,     0,    13,    20,     0,
      17,    16,    15,    18,    19,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    24,    25,    22,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    28,    38,    75,    81,
     126,     0,    21,     0,    36,     0,     0,   149,    35,   146,
      29,   128,    30,   131,    34,   143,    33,   140,    31,   134,
      32,   137,    39,    83,    40,    86,    41,    89,    42,    92,
       0,    76,    77,    78,    79,   127,   129,   130,   132,   142,
     144,   139,   141,   133,   135,   136,   138,    82,    84,    85,
      87,    88,    90,    91,    93,   148,   150,   145,   147,     0,
       0,    43,    94,    44,    95,    45,    96,    46,    97,    47,
      98,    48,    99,    49,   100,    50,   101,    51,   102,    52,
     103,    53,   104,    54,   105,    55,   106,    56,   107,    57,
     108,    58,   109,    59,   110,    60,   111,    61,   112,    62,
     113,    63,   114,    64,   115,    65,   116,    66,   117,    67,
     118,    68,   119,    69,   120,    70,   121,    71,   122,    72,
     123,    73,   124,    26,    27
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -41,   -41,     5,   -41,   -41,   -41,   -41,    33,   -36,   -41,
     -41,   140,   -41,   -41,   108,   -15,   -31,    68,   110
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,     4,     7,    10,     5,    11,    12,    57,    58,
      59,   103,   104,   324,    60,   205,    62,   206,    64
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      61,   105,   106,   117,   118,   119,   120,   107,   108,   109,
     110,     9,     6,   111,   112,   113,   114,   115,   116,    17,
       2,    13,    14,   102,     3,    98,   105,   106,   -12,    22,
     131,   132,     8,     8,   101,   201,    15,   131,   132,    65,
     105,    97,   115,   116,    61,   131,   198,    18,   196,    66,
     134,   136,   138,   140,   142,   144,   146,   148,   150,   152,
     154,   156,   158,   160,   162,   164,   166,   168,   170,   172,
     174,   176,   178,   180,   182,   184,   186,   188,   190,   192,
     194,   230,    16,    63,   200,    67,   231,   232,   233,   234,
     129,   130,   210,   212,   214,   216,   218,   220,   222,   224,
     226,   228,   117,   118,   119,   120,   235,   237,   239,   241,
     243,   245,   247,   249,   251,   253,   119,   120,   259,   113,
     114,   115,   116,    68,    99,    69,    70,    63,   113,   114,
     115,   116,    71,   135,   137,   139,   141,   143,   145,   147,
     149,   151,   153,   155,   157,   159,   161,   163,   165,   167,
     169,   171,   173,   175,   177,   179,   181,   183,   185,   187,
     189,   191,   193,   195,    96,   323,   100,   127,   128,   129,
     130,   261,    72,    73,    74,   211,   213,   215,   217,   219,
     221,   223,   225,   227,   229,    75,    61,    76,    97,   236,
     238,   240,   242,   244,   246,   248,   250,   252,   254,    98,
      77,   260,   113,   114,   115,   116,   127,   128,   129,   130,
     262,    78,    79,   204,   208,   207,   209,   107,   108,   109,
     110,    80,    81,   111,   112,   113,   114,   115,   116,   127,
     128,   129,   130,    82,    83,    84,    85,    86,    87,   255,
     257,   256,   258,    88,   202,   263,    89,    90,     0,   264,
      91,    92,    93,    94,    19,    95,    20,    21,    22,    23,
     113,   114,   115,   116,     0,     0,     0,     0,   197,    63,
       0,    99,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    53,    54,   265,     0,     0,     0,    55,   101,     0,
      56,    96,    19,   100,    20,    21,    22,    23,     0,     0,
       0,     0,     0,   127,   128,   129,   130,     0,     0,     0,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    53,
      54,    20,    21,    22,    23,    55,   266,     0,    56,   113,
     114,   115,   116,     0,     0,     0,     0,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    20,    21,
       0,    23,   267,     0,     0,    56,     0,     0,   127,   128,
     129,   130,     0,     0,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    20,    21,   121,   122,   123,
     124,   268,   203,   125,   126,   127,   128,   129,   130,     0,
       0,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,   113,   114,   115,   116,     0,     0,   199,   133,
     127,   128,   129,   130,   113,   114,   115,   116,   127,   128,
     129,   130,   113,   114,   115,   116,   127,   128,   129,   130,
     113,   114,   115,   116,   127,   128,   129,   130,   113,   114,
     115,   116,     0,     0,     0,   269,   127,   128,   129,   130,
       0,     0,     0,   270,     0,     0,     0,   271,     0,     0,
       0,   272,     0,     0,     0,   273,     0,     0,     0,   274,
       0,     0,     0,   275,     0,     0,     0,   276,     0,     0,
       0,   277,   113,   114,   115,   116,     0,     0,     0,   278,
     127,   128,   129,   130,   113,   114,   115,   116,   127,   128,
     129,   130,   113,   114,   115,   116,   127,   128,   129,   130,
     113,   114,   115,   116,   127,   128,   129,   130,   113,   114,
     115,   116,     0,     0,     0,   279,   127,   128,   129,   130,
       0,     0,     0,   280,     0,     0,     0,   281,     0,     0,
       0,   282,     0,     0,     0,   283,     0,     0,     0,   284,
       0,     0,     0,   285,     0,     0,     0,   286,     0,     0,
       0,   287,   113,   114,   115,   116,     0,     0,     0,   288,
     127,   128,   129,   130,   113,   114,   115,   116,   127,   128,
     129,   130,   113,   114,   115,   116,   127,   128,   129,   130,
     113,   114,   115,   116,   127,   128,   129,   130,   113,   114,
     115,   116,     0,     0,     0,   289,   127,   128,   129,   130,
       0,     0,     0,   290,     0,     0,     0,   291,     0,     0,
       0,   292,     0,     0,     0,   293,     0,     0,     0,   294,
       0,     0,     0,   295,     0,     0,     0,   296,     0,     0,
       0,   297,   113,   114,   115,   116,     0,     0,     0,   298,
     127,   128,   129,   130,   113,   114,   115,   116,   127,   128,
     129,   130,   113,   114,   115,   116,   127,   128,   129,   130,
     113,   114,   115,   116,   127,   128,   129,   130,   113,   114,
     115,   116,     0,     0,     0,   299,   127,   128,   129,   130,
       0,     0,     0,   300,     0,     0,     0,   301,     0,     0,
       0,   302,     0,     0,     0,   303,     0,     0,     0,   304,
       0,     0,     0,   305,     0,     0,     0,   306,     0,     0,
       0,   307,   113,   114,   115,   116,     0,     0,     0,   308,
     127,   128,   129,   130,   113,   114,   115,   116,   127,   128,
     129,   130,   113,   114,   115,   116,   127,   128,   129,   130,
     113,   114,   115,   116,   127,   128,   129,   130,   113,   114,
     115,   116,     0,     0,     0,   309,   127,   128,   129,   130,
       0,     0,     0,   310,     0,     0,     0,   311,     0,     0,
       0,   312,     0,     0,     0,   313,     0,     0,     0,   314,
       0,     0,     0,   315,     0,     0,     0,   316,     0,     0,
       0,   317,   113,   114,   115,   116,     0,     0,     0,   318,
     127,   128,   129,   130,   113,   114,   115,   116,   127,   128,
     129,   130,   113,   114,   115,   116,   127,   128,   129,   130,
     121,   122,   123,   124,     0,     0,   125,   126,   127,   128,
     129,   130,     0,     0,     0,   319,     0,     0,     0,     0,
       0,     0,     0,   320,     0,     0,     0,   321,     0,     0,
       0,   322,     0,     0,     0,   197,     0,     0,     0,   199
};

static const yytype_int16 yycheck[] =
{
      15,    13,    14,    17,    18,    19,    20,     9,    10,    11,
      12,     6,    52,    15,    16,    17,    18,    19,    20,    14,
       0,    53,    54,    59,     4,    56,    13,    14,    55,     7,
      13,    14,     4,     4,    57,    58,    55,    13,    14,    59,
      13,    56,    19,    20,    59,    13,    60,    14,    60,    59,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    59,    53,    15,    60,    59,   117,   118,   119,   120,
      19,    20,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,    17,    18,    19,    20,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,    19,    20,   133,    17,
      18,    19,    20,    59,    56,    59,    59,    59,    17,    18,
      19,    20,    59,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    78,    79,    80,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      92,    93,    94,    95,    56,   201,    56,    17,    18,    19,
      20,    60,    59,    59,    59,   107,   108,   109,   110,   111,
     112,   113,   114,   115,   116,    59,   201,    59,   203,   121,
     122,   123,   124,   125,   126,   127,   128,   129,   130,   230,
      59,   133,    17,    18,    19,    20,    17,    18,    19,    20,
      60,    59,    59,   105,   106,   105,   106,     9,    10,    11,
      12,    59,    59,    15,    16,    17,    18,    19,    20,    17,
      18,    19,    20,    59,    59,    59,    59,    59,    59,   131,
     132,   131,   132,    59,   104,    60,    59,    59,    -1,    60,
      59,    59,    59,    59,     3,    59,     5,     6,     7,     8,
      17,    18,    19,    20,    -1,    -1,    -1,    -1,    60,   201,
      -1,   203,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    60,    -1,    -1,    -1,    56,    57,    -1,
      59,   203,     3,   203,     5,     6,     7,     8,    -1,    -1,
      -1,    -1,    -1,    17,    18,    19,    20,    -1,    -1,    -1,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,     5,     6,     7,     8,    56,    60,    -1,    59,    17,
      18,    19,    20,    -1,    -1,    -1,    -1,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,     5,     6,
      -1,     8,    60,    -1,    -1,    59,    -1,    -1,    17,    18,
      19,    20,    -1,    -1,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,     5,     6,     9,    10,    11,
      12,    60,    59,    15,    16,    17,    18,    19,    20,    -1,
      -1,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    17,    18,    19,    20,    -1,    -1,    60,    59,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    17,    18,    19,    20,    17,    18,    19,    20,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    -1,    -1,    -1,    60,    17,    18,    19,    20,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    -1,    -1,    -1,    60,    -1,    -1,    -1,    60,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    17,    18,    19,    20,    -1,    -1,    -1,    60,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    17,    18,    19,    20,    17,    18,    19,    20,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    -1,    -1,    -1,    60,    17,    18,    19,    20,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    -1,    -1,    -1,    60,    -1,    -1,    -1,    60,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    17,    18,    19,    20,    -1,    -1,    -1,    60,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    17,    18,    19,    20,    17,    18,    19,    20,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    -1,    -1,    -1,    60,    17,    18,    19,    20,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    -1,    -1,    -1,    60,    -1,    -1,    -1,    60,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    17,    18,    19,    20,    -1,    -1,    -1,    60,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    17,    18,    19,    20,    17,    18,    19,    20,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    -1,    -1,    -1,    60,    17,    18,    19,    20,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    -1,    -1,    -1,    60,    -1,    -1,    -1,    60,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    17,    18,    19,    20,    -1,    -1,    -1,    60,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    17,    18,    19,    20,    17,    18,    19,    20,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    -1,    -1,    -1,    60,    17,    18,    19,    20,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    -1,    -1,    -1,    60,    -1,    -1,    -1,    60,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    17,    18,    19,    20,    -1,    -1,    -1,    60,
      17,    18,    19,    20,    17,    18,    19,    20,    17,    18,
      19,    20,    17,    18,    19,    20,    17,    18,    19,    20,
       9,    10,    11,    12,    -1,    -1,    15,    16,    17,    18,
      19,    20,    -1,    -1,    -1,    60,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    60,    -1,    -1,    -1,    60,    -1,    -1,
      -1,    60,    -1,    -1,    -1,    60,    -1,    -1,    -1,    60
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    62,     0,     4,    63,    66,    52,    64,     4,    63,
      65,    67,    68,    53,    54,    55,    53,    63,    68,     3,
       5,     6,     7,     8,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    56,    59,    69,    70,    71,
      75,    76,    77,    78,    79,    59,    59,    59,    59,    59,
      59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
      59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
      59,    59,    59,    59,    59,    59,    75,    76,    77,    78,
      79,    57,    69,    72,    73,    13,    14,     9,    10,    11,
      12,    15,    16,    17,    18,    19,    20,    17,    18,    19,
      20,     9,    10,    11,    12,    15,    16,    17,    18,    19,
      20,    13,    14,    59,    76,    78,    76,    78,    76,    78,
      76,    78,    76,    78,    76,    78,    76,    78,    76,    78,
      76,    78,    76,    78,    76,    78,    76,    78,    76,    78,
      76,    78,    76,    78,    76,    78,    76,    78,    76,    78,
      76,    78,    76,    78,    76,    78,    76,    78,    76,    78,
      76,    78,    76,    78,    76,    78,    76,    78,    76,    78,
      76,    78,    76,    78,    76,    78,    60,    60,    60,    60,
      60,    58,    72,    59,    75,    76,    78,    79,    75,    79,
      76,    78,    76,    78,    76,    78,    76,    78,    76,    78,
      76,    78,    76,    78,    76,    78,    76,    78,    76,    78,
      59,    77,    77,    77,    77,    76,    78,    76,    78,    76,
      78,    76,    78,    76,    78,    76,    78,    76,    78,    76,
      78,    76,    78,    76,    78,    75,    79,    75,    79,    76,
      78,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      60,    60,    60,    69,    74
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    61,    62,    62,    63,    64,    64,    65,    65,    65,
      65,    66,    67,    68,    69,    69,    69,    69,    69,    69,
      69,    70,    70,    71,    72,    73,    74,    73,    75,    75,
      75,    75,    75,    75,    75,    75,    75,    75,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    77,    77,    77,    77,    77,
      77,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    79,
      79
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     2,     3,     4,     1,     3,     1,
       3,     1,     1,     3,     1,     1,     1,     1,     1,     1,
       1,     3,     2,     1,     1,     1,     0,     4,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     3,     3,
       3,     3,     3,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     1,     3,     3,     3,     3,     3,
       1,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
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

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



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

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
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
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
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
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
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


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
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
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

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

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
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

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

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
| yyreduce -- Do a reduction.  |
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
        case 3:
#line 444 "build/Cello/parse.y" /* yacc.c:1646  */
    { }
#line 1960 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 4:
#line 448 "build/Cello/parse.y" /* yacc.c:1646  */
    {  }
#line 1966 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 5:
#line 451 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_group[--current_group_level] = 0; }
#line 1972 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 6:
#line 452 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_group[--current_group_level] = 0; }
#line 1978 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 7:
#line 455 "build/Cello/parse.y" /* yacc.c:1646  */
    {  }
#line 1984 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 8:
#line 456 "build/Cello/parse.y" /* yacc.c:1646  */
    {  }
#line 1990 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 9:
#line 457 "build/Cello/parse.y" /* yacc.c:1646  */
    {  }
#line 1996 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 10:
#line 458 "build/Cello/parse.y" /* yacc.c:1646  */
    {  }
#line 2002 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 11:
#line 462 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_group[current_group_level++] = (yyvsp[0].string_type); }
#line 2008 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 12:
#line 465 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_parameter = (yyvsp[0].string_type);}
#line 2014 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 13:
#line 468 "build/Cello/parse.y" /* yacc.c:1646  */
    { new_parameter(); }
#line 2020 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 14:
#line 472 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_type = enum_parameter_string;       yylval.string_type = (yyvsp[0].string_type); }
#line 2026 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 15:
#line 473 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_type = enum_parameter_integer;      yylval.integer_type = (yyvsp[0].integer_type);}
#line 2032 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 16:
#line 474 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_type = enum_parameter_float;       yylval.float_type = (yyvsp[0].float_type);}
#line 2038 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 17:
#line 475 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_type = enum_parameter_logical;      yylval.logical_type = (yyvsp[0].logical_type); }
#line 2044 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 18:
#line 476 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_type = enum_parameter_float_expr;  yylval.node_type = (yyvsp[0].node_type); }
#line 2050 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 19:
#line 477 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_type = enum_parameter_logical_expr; yylval.node_type = (yyvsp[0].node_type); }
#line 2056 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 20:
#line 478 "build/Cello/parse.y" /* yacc.c:1646  */
    { current_type = enum_parameter_list; }
#line 2062 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 21:
#line 481 "build/Cello/parse.y" /* yacc.c:1646  */
    {  }
#line 2068 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 22:
#line 482 "build/Cello/parse.y" /* yacc.c:1646  */
    {  }
#line 2074 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 23:
#line 485 "build/Cello/parse.y" /* yacc.c:1646  */
    { 
   struct param_struct * p = new_param_sentinel();
   p->list_value = param_curr;
   new_param_list(p);
   param_curr = p;
 }
#line 2085 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 24:
#line 492 "build/Cello/parse.y" /* yacc.c:1646  */
    { param_curr = param_curr->list_value; }
#line 2091 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 25:
#line 496 "build/Cello/parse.y" /* yacc.c:1646  */
    { new_parameter(); }
#line 2097 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 497 "build/Cello/parse.y" /* yacc.c:1646  */
    { new_parameter(); }
#line 2103 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 27:
#line 499 "build/Cello/parse.y" /* yacc.c:1646  */
    { }
#line 2109 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 28:
#line 504 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-1].logical_type); }
#line 2115 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 29:
#line 505 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-2].float_type) <= (yyvsp[0].float_type); }
#line 2121 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 30:
#line 506 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-2].float_type) >= (yyvsp[0].float_type); }
#line 2127 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 31:
#line 507 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-2].float_type) <  (yyvsp[0].float_type); }
#line 2133 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 32:
#line 508 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-2].float_type) >  (yyvsp[0].float_type); }
#line 2139 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 33:
#line 509 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-2].float_type) == (yyvsp[0].float_type); }
#line 2145 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 34:
#line 510 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-2].float_type) != (yyvsp[0].float_type); }
#line 2151 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 35:
#line 511 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-2].logical_type) || (yyvsp[0].logical_type); }
#line 2157 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 36:
#line 512 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[-2].logical_type) && (yyvsp[0].logical_type); }
#line 2163 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 37:
#line 513 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.logical_type) = (yyvsp[0].logical_type); }
#line 2169 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 38:
#line 517 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = (yyvsp[-1].float_type); }
#line 2175 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 39:
#line 518 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = (yyvsp[-2].float_type) + (yyvsp[0].float_type);}
#line 2181 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 40:
#line 519 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = (yyvsp[-2].float_type) - (yyvsp[0].float_type);}
#line 2187 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 41:
#line 520 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = (yyvsp[-2].float_type) * (yyvsp[0].float_type);}
#line 2193 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 42:
#line 521 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = (yyvsp[-2].float_type) / (yyvsp[0].float_type);}
#line 2199 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 43:
#line 522 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = acos((yyvsp[-1].float_type)); }
#line 2205 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 44:
#line 523 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = acosh((yyvsp[-1].float_type)); }
#line 2211 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 45:
#line 524 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = asin((yyvsp[-1].float_type)); }
#line 2217 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 46:
#line 525 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = asinh((yyvsp[-1].float_type)); }
#line 2223 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 47:
#line 526 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = atan((yyvsp[-1].float_type)); }
#line 2229 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 48:
#line 527 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = atanh((yyvsp[-1].float_type)); }
#line 2235 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 49:
#line 528 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = cbrt((yyvsp[-1].float_type)); }
#line 2241 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 50:
#line 529 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = ceil((yyvsp[-1].float_type)); }
#line 2247 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 51:
#line 530 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = cos((yyvsp[-1].float_type)); }
#line 2253 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 52:
#line 531 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = cosh((yyvsp[-1].float_type)); }
#line 2259 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 53:
#line 532 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = erfc((yyvsp[-1].float_type)); }
#line 2265 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 54:
#line 533 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = erf((yyvsp[-1].float_type)); }
#line 2271 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 55:
#line 534 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = exp((yyvsp[-1].float_type)); }
#line 2277 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 56:
#line 535 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = expm1((yyvsp[-1].float_type)); }
#line 2283 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 57:
#line 536 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = fabs((yyvsp[-1].float_type)); }
#line 2289 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 58:
#line 537 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = floor((yyvsp[-1].float_type)); }
#line 2295 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 59:
#line 539 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = j0((yyvsp[-1].float_type)); }
#line 2301 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 60:
#line 540 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = j1((yyvsp[-1].float_type)); }
#line 2307 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 61:
#line 541 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = lgamma((yyvsp[-1].float_type)); }
#line 2313 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 62:
#line 542 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = log10((yyvsp[-1].float_type)); }
#line 2319 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 63:
#line 543 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = log1p((yyvsp[-1].float_type)); }
#line 2325 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 64:
#line 544 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = logb((yyvsp[-1].float_type)); }
#line 2331 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 65:
#line 545 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = log((yyvsp[-1].float_type)); }
#line 2337 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 66:
#line 546 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = sin((yyvsp[-1].float_type)); }
#line 2343 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 67:
#line 547 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = sinh((yyvsp[-1].float_type)); }
#line 2349 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 68:
#line 548 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = sqrt((yyvsp[-1].float_type)); }
#line 2355 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 69:
#line 549 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = tan((yyvsp[-1].float_type)); }
#line 2361 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 70:
#line 550 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = tanh((yyvsp[-1].float_type)); }
#line 2367 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 71:
#line 551 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = y0((yyvsp[-1].float_type)); }
#line 2373 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 72:
#line 552 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = y1((yyvsp[-1].float_type)); }
#line 2379 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 73:
#line 553 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = rint((yyvsp[-1].float_type)); }
#line 2385 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 74:
#line 554 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.float_type) = (yyvsp[0].float_type);}
#line 2391 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 75:
#line 558 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.integer_type) = (yyvsp[-1].integer_type); }
#line 2397 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 76:
#line 559 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.integer_type) = (yyvsp[-2].integer_type) + (yyvsp[0].integer_type);}
#line 2403 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 77:
#line 560 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.integer_type) = (yyvsp[-2].integer_type) - (yyvsp[0].integer_type);}
#line 2409 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 78:
#line 561 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.integer_type) = (yyvsp[-2].integer_type) * (yyvsp[0].integer_type);}
#line 2415 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 79:
#line 562 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.integer_type) = (yyvsp[-2].integer_type) / (yyvsp[0].integer_type);}
#line 2421 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 80:
#line 563 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.integer_type) = (yyvsp[0].integer_type);}
#line 2427 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 81:
#line 567 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = (yyvsp[-1].node_type); }
#line 2433 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 82:
#line 568 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_add,new_node_float((yyvsp[0].float_type))); }
#line 2439 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 83:
#line 569 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_add,(yyvsp[0].node_type)); }
#line 2445 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 84:
#line 570 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_add,(yyvsp[0].node_type)); }
#line 2451 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 85:
#line 571 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_sub,new_node_float((yyvsp[0].float_type))); }
#line 2457 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 86:
#line 572 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_sub,(yyvsp[0].node_type)); }
#line 2463 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 87:
#line 573 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_sub,(yyvsp[0].node_type)); }
#line 2469 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 88:
#line 574 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_mul,new_node_float((yyvsp[0].float_type))); }
#line 2475 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 89:
#line 575 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_mul,(yyvsp[0].node_type)); }
#line 2481 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 90:
#line 576 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_mul,(yyvsp[0].node_type)); }
#line 2487 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 91:
#line 577 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_div,new_node_float((yyvsp[0].float_type))); }
#line 2493 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 92:
#line 578 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_div,(yyvsp[0].node_type)); }
#line 2499 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 93:
#line 579 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_div,(yyvsp[0].node_type)); }
#line 2505 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 94:
#line 580 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( acos, "acos", (yyvsp[-1].node_type)); }
#line 2511 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 95:
#line 581 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( acosh, "acosh", (yyvsp[-1].node_type)); }
#line 2517 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 96:
#line 582 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( asin, "asin", (yyvsp[-1].node_type)); }
#line 2523 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 97:
#line 583 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( asinh, "asinh", (yyvsp[-1].node_type)); }
#line 2529 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 98:
#line 584 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( atan, "atan", (yyvsp[-1].node_type)); }
#line 2535 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 99:
#line 585 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( atanh, "atanh", (yyvsp[-1].node_type)); }
#line 2541 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 100:
#line 586 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( cbrt, "cbrt", (yyvsp[-1].node_type)); }
#line 2547 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 101:
#line 587 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( ceil, "ceil", (yyvsp[-1].node_type)); }
#line 2553 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 102:
#line 588 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( cos, "cos", (yyvsp[-1].node_type)); }
#line 2559 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 103:
#line 589 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( cosh, "cosh", (yyvsp[-1].node_type)); }
#line 2565 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 104:
#line 590 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( erfc, "erfc", (yyvsp[-1].node_type)); }
#line 2571 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 105:
#line 591 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( erf, "erf", (yyvsp[-1].node_type)); }
#line 2577 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 106:
#line 592 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( exp, "exp", (yyvsp[-1].node_type)); }
#line 2583 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 107:
#line 593 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( expm1, "expm1", (yyvsp[-1].node_type)); }
#line 2589 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 108:
#line 594 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( fabs, "fabs", (yyvsp[-1].node_type)); }
#line 2595 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 109:
#line 595 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( floor, "floor", (yyvsp[-1].node_type)); }
#line 2601 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 110:
#line 597 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( j0, "j0", (yyvsp[-1].node_type)); }
#line 2607 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 111:
#line 598 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( j1, "j1", (yyvsp[-1].node_type)); }
#line 2613 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 112:
#line 599 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( lgamma, "lgamma", (yyvsp[-1].node_type)); }
#line 2619 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 113:
#line 600 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( log10, "log10", (yyvsp[-1].node_type)); }
#line 2625 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 114:
#line 601 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( log1p, "log1p", (yyvsp[-1].node_type)); }
#line 2631 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 115:
#line 602 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( logb, "logb", (yyvsp[-1].node_type)); }
#line 2637 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 116:
#line 603 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( log, "log", (yyvsp[-1].node_type)); }
#line 2643 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 117:
#line 604 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( sin, "sin", (yyvsp[-1].node_type)); }
#line 2649 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 118:
#line 605 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( sinh, "sinh", (yyvsp[-1].node_type)); }
#line 2655 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 119:
#line 606 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( sqrt, "sqrt", (yyvsp[-1].node_type)); }
#line 2661 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 120:
#line 607 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( tan, "tan", (yyvsp[-1].node_type)); }
#line 2667 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 121:
#line 608 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( tanh, "tanh", (yyvsp[-1].node_type)); }
#line 2673 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 122:
#line 609 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( y0, "y0", (yyvsp[-1].node_type)); }
#line 2679 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 123:
#line 610 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( y1, "y1", (yyvsp[-1].node_type)); }
#line 2685 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 124:
#line 611 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_function ( rint, "rint", (yyvsp[-1].node_type)); }
#line 2691 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 125:
#line 612 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_variable ((yyvsp[0].string_type));  }
#line 2697 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 126:
#line 617 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = (yyvsp[-1].node_type); }
#line 2703 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 127:
#line 618 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_le,new_node_float((yyvsp[0].float_type))); }
#line 2709 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 128:
#line 619 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_le,(yyvsp[0].node_type)); }
#line 2715 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 129:
#line 620 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_le,(yyvsp[0].node_type)); }
#line 2721 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 130:
#line 621 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_ge,new_node_float((yyvsp[0].float_type))); }
#line 2727 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 131:
#line 622 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_ge,(yyvsp[0].node_type)); }
#line 2733 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 132:
#line 623 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_ge,(yyvsp[0].node_type)); }
#line 2739 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 133:
#line 624 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_lt,new_node_float((yyvsp[0].float_type))); }
#line 2745 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 134:
#line 625 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_lt,(yyvsp[0].node_type)); }
#line 2751 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 135:
#line 626 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_lt,(yyvsp[0].node_type)); }
#line 2757 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 136:
#line 627 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_gt,new_node_float((yyvsp[0].float_type))); }
#line 2763 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 137:
#line 628 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_gt,(yyvsp[0].node_type)); }
#line 2769 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 138:
#line 629 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_gt,(yyvsp[0].node_type)); }
#line 2775 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 139:
#line 630 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_eq,new_node_float((yyvsp[0].float_type))); }
#line 2781 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 140:
#line 631 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_eq,(yyvsp[0].node_type)); }
#line 2787 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 141:
#line 632 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_eq,(yyvsp[0].node_type)); }
#line 2793 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 142:
#line 633 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_ne,new_node_float((yyvsp[0].float_type))); }
#line 2799 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 143:
#line 634 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_float((yyvsp[-2].float_type)), enum_op_ne,(yyvsp[0].node_type)); }
#line 2805 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 144:
#line 635 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_ne,(yyvsp[0].node_type)); }
#line 2811 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 145:
#line 636 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_or,new_node_logical((yyvsp[0].logical_type))); }
#line 2817 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 146:
#line 637 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_logical((yyvsp[-2].logical_type)), enum_op_or,(yyvsp[0].node_type)); }
#line 2823 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 147:
#line 638 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_or,(yyvsp[0].node_type)); }
#line 2829 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 148:
#line 639 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_and,new_node_logical((yyvsp[0].logical_type))); }
#line 2835 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 149:
#line 640 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation (new_node_logical((yyvsp[-2].logical_type)), enum_op_and,(yyvsp[0].node_type)); }
#line 2841 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;

  case 150:
#line 641 "build/Cello/parse.y" /* yacc.c:1646  */
    { (yyval.node_type) = new_node_operation ((yyvsp[-2].node_type), enum_op_and,(yyvsp[0].node_type)); }
#line 2847 "build/Cello/parse.tab.c" /* yacc.c:1646  */
    break;


#line 2851 "build/Cello/parse.tab.c" /* yacc.c:1646  */
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
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
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

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

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

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
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
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
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
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 645 "build/Cello/parse.y" /* yacc.c:1906  */


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

void sprintf_expression (struct node_expr * node,
			 char * buffer)
/* WARNING: buffer is assumed to be big enough to hold the expression */
{
  if (node == NULL) {
    sprintf (buffer,"NULL");
  } else {
    char left,right;
    switch (node->type) {
    case enum_node_integer:
      sprintf (buffer,"%d",node->integer_value);
      buffer += strlen(buffer);
      break;
    case enum_node_float:
      /* '#' format character forces a decimal point */
      sprintf (buffer,FLOAT_FORMAT,node->float_value);
      buffer += strlen(buffer);
      break;
    case enum_node_variable:
      sprintf (buffer,"%c",node->var_value);
      buffer += strlen(buffer);
      break;
    case enum_node_function:
      sprintf (buffer,"%s(",node->function_name);
      buffer += strlen(buffer);
      sprintf_expression(node->left,buffer+strlen(buffer));
      buffer += strlen(buffer);
      sprintf (buffer,")");
      buffer += strlen(buffer);
      break;
    case enum_node_operation:
      left  = (node->left->type == enum_node_operation) ? '(' : ' ';
      right = (node->left->type == enum_node_operation) ? ')' : ' ';
      sprintf (buffer,"%c",left);
      buffer += strlen(buffer);
      sprintf_expression(node->left,buffer+strlen(buffer));
      buffer += strlen(buffer);
      sprintf (buffer,"%c",right);
      buffer += strlen(buffer);
      sprintf (buffer," %s ",op_name[node->op_value]);
      buffer += strlen(buffer);
      left  = (node->right->type == enum_node_operation) ? '(' : ' ';
      right = (node->right->type == enum_node_operation) ? ')' : ' ';
      sprintf (buffer,"%c",left);
      buffer += strlen(buffer);
      sprintf_expression(node->right,buffer+strlen(buffer));
      buffer += strlen(buffer);
      sprintf (buffer,"%c",right);
      buffer += strlen(buffer);
      break;
    default:
      break;
    }
  }
}

void cello_parameters_print_list(struct param_struct * head, int level)
{
  struct param_struct * p = head->next;
  int i;

  while (p && p->type != enum_parameter_sentinel) {

    if (p->group != NULL) {
      indent(level);
      printf ("%s ", parameter_name[p->type]);
      for (i=0; p->group[i] != NULL && i < MAX_GROUP_DEPTH; i++) {
	printf ("%s:",p->group[i]);
      }
      printf ("%s = ", p->parameter);
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
    p = p->next;
  }
}

void cello_parameters_print()
{
  cello_parameters_print_list(param_head,0);
}

