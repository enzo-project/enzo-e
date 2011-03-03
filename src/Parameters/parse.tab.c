
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
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
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 1 "src/Parameters/parse.y"

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

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>

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
  "node_scalar",
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
    "subgroup",
    "integer",
    "scalar",
    "string",
    "identifier",
    "logical",
    "list",
    "scalar_expr",
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

  struct node_expr * new_node_scalar (double value)
  {
    struct node_expr * node = malloc (sizeof (struct node_expr));

    node->type          = enum_node_scalar;
    node->scalar_value  = value;
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

  /* The current group, subgroup, and parameter type */

  char *              current_parameter = NULL;
  char *              current_group     = NULL;
  char *              current_subgroup  = NULL;
  enum enum_parameter current_type      = enum_parameter_sentinel;

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

  /* Function to update parameter's subgroups once the subgroup is known */

/*   void update_subgroup (char * subgroup) */
/*     { */
/*       struct param_struct * p = param_curr; */
/*       int inside_subgroup = 1; */
/*       while (p->next->type     != enum_parameter_sentinel &&  */
/* 	     p->next->subgroup == NULL) { */
/* 	if (p->next->type == enum_parameter_subgroup) { */
/* 	  inside_subgroup = 0; */
/*           delete_param(p); */
/*         } else if (inside_subgroup) { */
/*           p->next->subgroup = strdup(subgroup); */
/*           p = p -> next; */
/*         } */
/*       } */
/*     } */

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
     p->group     = (current_group)     ? strdup(current_group)     : 0;
     p->subgroup  = (current_subgroup)  ? strdup(current_subgroup)  : 0;
     p->parameter = (current_parameter) ? strdup(current_parameter) : 0;

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


  /* New scalar parameter assignment */

  void new_param_scalar (double value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_scalar;
    p->scalar_value = value;
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
  void new_param_subgroup (char * value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_subgroup;
    p->string_value = value;
  }

  /* New empty parameter assignment: FIRST NODE IN LIST IS A SENTINEL  */
  struct param_struct * new_param_sentinel ()
  {

    /* MEMORY LEAK */
    struct param_struct * p = 
      (struct param_struct *) malloc (sizeof (struct param_struct));

    p->group     = NULL;
    p->subgroup  = NULL;
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
     case enum_parameter_subgroup:
       new_param_subgroup(yylval.subgroup_type);
       break;
     case enum_parameter_integer:
       new_param_integer(yylval.integer_type);
       break;
     case enum_parameter_scalar:
       new_param_scalar(yylval.scalar_type);
       break;
     case enum_parameter_string: 
       new_param_string(yylval.string_type);
       break;
     case enum_parameter_logical:
       new_param_logical(yylval.logical_type);
       break;
     case enum_parameter_list:
       break;
     case enum_parameter_scalar_expr:
       new_param_expr(enum_parameter_scalar_expr,yylval.node_type);
       break;
     case enum_parameter_logical_expr:
       new_param_expr(enum_parameter_logical_expr,yylval.node_type);
       break;
    default:
       printf ("%s:%d Parse Error: unknown type %d\n",
	       __FILE__,__LINE__,current_type);
       break;
     }
  }

  char * strcat3 (const char * s1,const char * s2,const char * s3)
  {
    char * s = malloc (strlen(s1) + strlen(s2) + strlen(s3) + 1);

    strcpy(s,s1);
    strcpy(s+strlen(s1),s2);
    strcpy(s+strlen(s1)+strlen(s2),s3);
    return s;
  }

  char * ftoa (double f)
    { 
      char * a = malloc(25); 

      sprintf (a,"%24.16e",f);
      return a;
    }



/* Line 189 of yacc.c  */
#line 460 "src/Parameters/parse.tab.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     GROUP_NAME = 258,
     STRING = 259,
     IDENTIFIER = 260,
     VARIABLE = 261,
     SCALAR = 262,
     INTEGER = 263,
     LOGICAL = 264,
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



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 388 "src/Parameters/parse.y"
 
  int logical_type;  
  int integer_type; 
  double scalar_type;  
  char * string_type; 
  char * subgroup_type;
  struct node_expr * node_type;
  


/* Line 214 of yacc.c  */
#line 553 "src/Parameters/parse.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 565 "src/Parameters/parse.tab.c"

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
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
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
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
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
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
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
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
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

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   930

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  62
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  22
/* YYNRULES -- Number of rules.  */
#define YYNRULES  152
/* YYNRULES -- Number of states.  */
#define YYNSTATES  328

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   301

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      60,    61,    20,    18,    59,    19,     2,    21,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    55,
      16,    56,    17,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    57,     2,    58,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    53,     2,    54,     2,     2,     2,     2,
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
      15,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,    52
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,    10,    13,    16,    20,    25,
      27,    31,    33,    35,    37,    39,    41,    45,    47,    49,
      51,    53,    55,    57,    59,    63,    65,    67,    69,    70,
      75,    79,    83,    87,    91,    95,    99,   103,   107,   111,
     113,   117,   121,   125,   129,   133,   138,   143,   148,   153,
     158,   163,   168,   173,   178,   183,   188,   193,   198,   203,
     208,   213,   218,   223,   228,   233,   238,   243,   248,   253,
     258,   263,   268,   273,   278,   283,   288,   290,   294,   298,
     302,   306,   310,   312,   316,   320,   324,   328,   332,   336,
     340,   344,   348,   352,   356,   360,   364,   369,   374,   379,
     384,   389,   394,   399,   404,   409,   414,   419,   424,   429,
     434,   439,   444,   449,   454,   459,   464,   469,   474,   479,
     484,   489,   494,   499,   504,   509,   514,   519,   521,   525,
     529,   533,   537,   541,   545,   549,   553,   557,   561,   565,
     569,   573,   577,   581,   585,   589,   593,   597,   601,   605,
     609,   613,   617
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      63,     0,    -1,    -1,    63,    64,    -1,    69,    66,    -1,
      69,    65,    -1,    70,    66,    -1,    53,    67,    54,    -1,
      53,    67,    55,    54,    -1,    68,    -1,    67,    55,    68,
      -1,    72,    -1,    65,    -1,     3,    -1,     5,    -1,     5,
      -1,    71,    56,    73,    -1,     4,    -1,    81,    -1,    80,
      -1,    79,    -1,    82,    -1,    83,    -1,    74,    -1,    75,
      77,    76,    -1,    57,    -1,    58,    -1,    73,    -1,    -1,
      77,    59,    73,    78,    -1,    60,    79,    61,    -1,    80,
      10,    80,    -1,    80,    11,    80,    -1,    80,    16,    80,
      -1,    80,    17,    80,    -1,    80,    13,    80,    -1,    80,
      12,    80,    -1,    79,    15,    79,    -1,    79,    14,    79,
      -1,     9,    -1,    60,    80,    61,    -1,    80,    18,    80,
      -1,    80,    19,    80,    -1,    80,    20,    80,    -1,    80,
      21,    80,    -1,    22,    60,    80,    61,    -1,    23,    60,
      80,    61,    -1,    24,    60,    80,    61,    -1,    25,    60,
      80,    61,    -1,    26,    60,    80,    61,    -1,    27,    60,
      80,    61,    -1,    28,    60,    80,    61,    -1,    29,    60,
      80,    61,    -1,    30,    60,    80,    61,    -1,    31,    60,
      80,    61,    -1,    32,    60,    80,    61,    -1,    33,    60,
      80,    61,    -1,    34,    60,    80,    61,    -1,    35,    60,
      80,    61,    -1,    36,    60,    80,    61,    -1,    37,    60,
      80,    61,    -1,    38,    60,    80,    61,    -1,    39,    60,
      80,    61,    -1,    40,    60,    80,    61,    -1,    41,    60,
      80,    61,    -1,    42,    60,    80,    61,    -1,    43,    60,
      80,    61,    -1,    44,    60,    80,    61,    -1,    45,    60,
      80,    61,    -1,    46,    60,    80,    61,    -1,    47,    60,
      80,    61,    -1,    48,    60,    80,    61,    -1,    49,    60,
      80,    61,    -1,    50,    60,    80,    61,    -1,    51,    60,
      80,    61,    -1,    52,    60,    80,    61,    -1,     7,    -1,
      60,    81,    61,    -1,    81,    18,    81,    -1,    81,    19,
      81,    -1,    81,    20,    81,    -1,    81,    21,    81,    -1,
       8,    -1,    60,    82,    61,    -1,    82,    18,    80,    -1,
      80,    18,    82,    -1,    82,    18,    82,    -1,    82,    19,
      80,    -1,    80,    19,    82,    -1,    82,    19,    82,    -1,
      82,    20,    80,    -1,    80,    20,    82,    -1,    82,    20,
      82,    -1,    82,    21,    80,    -1,    80,    21,    82,    -1,
      82,    21,    82,    -1,    22,    60,    82,    61,    -1,    23,
      60,    82,    61,    -1,    24,    60,    82,    61,    -1,    25,
      60,    82,    61,    -1,    26,    60,    82,    61,    -1,    27,
      60,    82,    61,    -1,    28,    60,    82,    61,    -1,    29,
      60,    82,    61,    -1,    30,    60,    82,    61,    -1,    31,
      60,    82,    61,    -1,    32,    60,    82,    61,    -1,    33,
      60,    82,    61,    -1,    34,    60,    82,    61,    -1,    35,
      60,    82,    61,    -1,    36,    60,    82,    61,    -1,    37,
      60,    82,    61,    -1,    38,    60,    82,    61,    -1,    39,
      60,    82,    61,    -1,    40,    60,    82,    61,    -1,    41,
      60,    82,    61,    -1,    42,    60,    82,    61,    -1,    43,
      60,    82,    61,    -1,    44,    60,    82,    61,    -1,    45,
      60,    82,    61,    -1,    46,    60,    82,    61,    -1,    47,
      60,    82,    61,    -1,    48,    60,    82,    61,    -1,    49,
      60,    82,    61,    -1,    50,    60,    82,    61,    -1,    51,
      60,    82,    61,    -1,    52,    60,    82,    61,    -1,     6,
      -1,    60,    83,    61,    -1,    82,    10,    80,    -1,    80,
      10,    82,    -1,    82,    10,    82,    -1,    82,    11,    80,
      -1,    80,    11,    82,    -1,    82,    11,    82,    -1,    82,
      16,    80,    -1,    80,    16,    82,    -1,    82,    16,    82,
      -1,    82,    17,    80,    -1,    80,    17,    82,    -1,    82,
      17,    82,    -1,    82,    13,    80,    -1,    80,    13,    82,
      -1,    82,    13,    82,    -1,    82,    12,    80,    -1,    80,
      12,    82,    -1,    82,    12,    82,    -1,    83,    15,    79,
      -1,    79,    15,    83,    -1,    83,    15,    83,    -1,    83,
      14,    79,    -1,    79,    14,    83,    -1,    83,    14,    83,
      -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   467,   467,   468,   472,   474,   478,   481,   482,   485,
     486,   489,   490,   493,   497,   500,   503,   507,   508,   509,
     510,   511,   512,   513,   516,   519,   526,   530,   531,   531,
     538,   539,   540,   541,   542,   543,   544,   545,   546,   547,
     551,   552,   553,   554,   555,   556,   557,   558,   559,   560,
     561,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   573,   574,   575,   576,   577,   578,   579,   580,   581,
     582,   583,   584,   585,   586,   587,   588,   592,   593,   594,
     595,   596,   597,   601,   602,   603,   604,   605,   606,   607,
     608,   609,   610,   611,   612,   613,   614,   615,   616,   617,
     618,   619,   620,   621,   622,   623,   624,   625,   626,   627,
     628,   629,   631,   632,   633,   634,   635,   636,   637,   638,
     639,   640,   641,   642,   643,   644,   645,   646,   651,   652,
     653,   654,   655,   656,   657,   658,   659,   660,   661,   662,
     663,   664,   665,   666,   667,   668,   669,   670,   671,   672,
     673,   674,   675
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "GROUP_NAME", "STRING", "IDENTIFIER",
  "VARIABLE", "SCALAR", "INTEGER", "LOGICAL", "LE", "GE", "NE", "EQ",
  "AND", "OR", "'<'", "'>'", "'+'", "'-'", "'*'", "'/'", "ACOS", "ACOSH",
  "ASIN", "ASINH", "ATAN", "ATANH", "CBRT", "CEIL", "COS", "COSH", "ERFC",
  "ERF", "EXP", "EXPM1", "FABS", "FLOOR", "J0", "J1", "LGAMMA", "LOG10",
  "LOG1P", "LOGB", "LOG", "SIN", "SINH", "SQRT", "TAN", "TANH", "Y0", "Y1",
  "RINT", "'{'", "'}'", "';'", "'='", "'['", "']'", "','", "'('", "')'",
  "$accept", "file", "group", "named_parameter_group", "parameter_group",
  "parameter_list", "parameter_item", "group_name", "subgroup_name",
  "parameter_name", "parameter_assignment", "parameter_value", "list",
  "LIST_BEGIN", "LIST_END", "list_elements", "$@1", "cle", "cse", "cie",
  "vse", "vle", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,    60,    62,    43,    45,
      42,    47,   271,   272,   273,   274,   275,   276,   277,   278,
     279,   280,   281,   282,   283,   284,   285,   286,   287,   288,
     289,   290,   291,   292,   293,   294,   295,   296,   297,   298,
     299,   300,   301,   123,   125,    59,    61,    91,    93,    44,
      40,    41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    62,    63,    63,    64,    64,    65,    66,    66,    67,
      67,    68,    68,    69,    70,    71,    72,    73,    73,    73,
      73,    73,    73,    73,    74,    75,    76,    77,    78,    77,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    79,
      80,    80,    80,    80,    80,    80,    80,    80,    80,    80,
      80,    80,    80,    80,    80,    80,    80,    80,    80,    80,
      80,    80,    80,    80,    80,    80,    80,    80,    80,    80,
      80,    80,    80,    80,    80,    80,    80,    81,    81,    81,
      81,    81,    81,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    83,    83,
      83,    83,    83,    83,    83,    83,    83,    83,    83,    83,
      83,    83,    83,    83,    83,    83,    83,    83,    83,    83,
      83,    83,    83
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     2,     2,     2,     3,     4,     1,
       3,     1,     1,     1,     1,     1,     3,     1,     1,     1,
       1,     1,     1,     1,     3,     1,     1,     1,     0,     4,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       3,     3,     3,     3,     3,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     1,     3,     3,     3,
       3,     3,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,    13,     3,     0,    14,     0,     5,     4,
       0,    14,    12,     0,     9,     0,    11,     6,     7,     0,
       0,     8,    10,    17,   127,    76,    82,    39,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    25,
       0,    16,    23,     0,    20,    19,    18,    21,    22,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    30,    40,
      77,    83,   128,    26,     0,    24,     0,    38,     0,     0,
     151,    37,   148,    31,   130,    32,   133,    36,   145,    35,
     142,    33,   136,    34,   139,    41,    85,    42,    88,    43,
      91,    44,    94,     0,    78,    79,    80,    81,   129,   131,
     132,   134,   144,   146,   141,   143,   135,   137,   138,   140,
      84,    86,    87,    89,    90,    92,    93,    95,   150,   152,
     147,   149,     0,     0,    45,    96,    46,    97,    47,    98,
      48,    99,    49,   100,    50,   101,    51,   102,    52,   103,
      53,   104,    54,   105,    55,   106,    56,   107,    57,   108,
      58,   109,    59,   110,    60,   111,    61,   112,    62,   113,
      63,   114,    64,   115,    65,   116,    66,   117,    67,   118,
      68,   119,    69,   120,    70,   121,    71,   122,    72,   123,
      73,   124,    74,   125,    75,   126,    28,    29
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,     4,    12,     9,    13,    14,     5,    10,    15,
      16,    61,    62,    63,   205,   106,   327,    64,   208,    66,
     209,    68
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -44
static const yytype_int16 yypact[] =
{
     -44,    23,   -44,   -44,   -44,    29,   -44,    19,   -44,   -44,
     -26,   -25,   -44,   -43,   -44,   -23,   -44,   -44,   -44,    27,
     249,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -16,    65,
      66,   108,   109,   110,   123,   125,   137,   138,   144,   145,
     160,   161,   162,   167,   168,   169,   170,   171,   172,   173,
     174,   175,   176,   177,   182,   183,   184,   185,   186,   -44,
     316,   -44,   -44,   249,    15,    -3,     1,   881,    24,   410,
     410,   410,   410,   410,   410,   410,   410,   410,   410,   410,
     410,   410,   410,   410,   410,   410,   410,   410,   410,   410,
     410,   410,   410,   410,   410,   410,   410,   410,   410,   410,
     -13,   190,   -15,   408,    22,   -44,   -17,   363,   363,   410,
     410,   410,   410,   410,   410,   410,   410,   410,   410,    20,
      20,    20,    20,   410,   410,   410,   410,   410,   410,   410,
     410,   410,   410,   363,   363,   410,    67,   111,   198,   242,
     246,   310,   314,   361,   445,   453,   457,   461,   465,   469,
     473,   477,   481,   489,   525,   533,   537,   541,   545,   549,
     553,   557,   561,   569,   605,   613,   617,   621,   625,   629,
     633,   637,   641,   649,   685,   693,   697,   701,   705,   709,
     713,   717,   721,   729,   765,   773,   777,   781,   785,   789,
     793,   797,   801,   809,   845,   853,   857,   861,   -44,   -44,
     -44,   -44,   -44,   -44,   249,   -44,   363,   -44,    -3,   881,
     -44,    21,    31,    81,    98,    81,    98,    81,    98,    81,
      98,    81,    98,    81,    98,    93,   100,    93,   100,   -44,
     -44,   -44,   -44,    20,   102,   102,   -44,   -44,    81,    98,
      81,    98,    81,    98,    81,    98,    81,    98,    81,    98,
      93,   100,    93,   100,   -44,   -44,   -44,   -44,   -44,   -44,
      21,    31,   865,   869,   -44,   -44,   -44,   -44,   -44,   -44,
     -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,
     -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,
     -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,
     -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,
     -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44,
     -44,   -44,   -44,   -44,   -44,   -44,   -44,   -44
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -44,   -44,   -44,    42,   237,   -44,   152,   -44,   -44,   -44,
     -44,   -38,   -44,   -44,   -44,   -44,   -44,   105,   -20,   104,
      64,   107
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -16
static const yytype_int16 yytable[] =
{
      65,   107,   108,   119,   120,   121,   122,   109,   110,   111,
     112,    18,    19,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   122,     2,    11,   105,     3,     7,    26,   107,
     108,   -15,    11,    20,     6,   107,   133,   134,   133,   134,
     101,   203,   204,    65,    69,   133,   200,     8,   198,   136,
     138,   140,   142,   144,   146,   148,   150,   152,   154,   156,
     158,   160,   162,   164,   166,   168,   170,   172,   174,   176,
     178,   180,   182,   184,   186,   188,   190,   192,   194,   196,
     233,    21,     7,   202,    67,   115,   116,   117,   118,   213,
     215,   217,   219,   221,   223,   225,   227,   229,   231,   115,
     116,   117,   118,   238,   240,   242,   244,   246,   248,   250,
     252,   254,   256,   117,   118,   262,   129,   130,   131,   132,
     131,   132,   121,   122,   103,    70,    71,    67,   264,   129,
     130,   131,   132,   137,   139,   141,   143,   145,   147,   149,
     151,   153,   155,   157,   159,   161,   163,   165,   167,   169,
     171,   173,   175,   177,   179,   181,   183,   185,   187,   189,
     191,   193,   195,   197,   102,   100,   326,   104,    72,    73,
      74,    22,   265,   214,   216,   218,   220,   222,   224,   226,
     228,   230,   232,    75,    65,    76,   101,   239,   241,   243,
     245,   247,   249,   251,   253,   255,   257,    77,    78,   263,
     109,   110,   111,   112,    79,    80,   113,   114,   115,   116,
     117,   118,   207,   211,   210,   212,   115,   116,   117,   118,
      81,    82,    83,   234,   235,   236,   237,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,   258,   260,
     259,   261,    95,    96,    97,    98,    99,    17,     0,     0,
       0,   199,     0,    23,     0,    24,    25,    26,    27,   266,
     129,   130,   131,   132,   115,   116,   117,   118,    67,     0,
     103,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    56,
      57,    58,     0,   267,     0,     0,    59,   268,     0,    60,
       0,   100,     0,   104,     0,     0,     0,     0,     0,     0,
       0,     0,    24,    25,    26,    27,     0,     0,   129,   130,
     131,   132,   115,   116,   117,   118,     0,   102,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    56,    57,    58,    24,
      25,   269,    27,     0,     0,   270,    60,     0,     0,   129,
     130,   131,   132,     0,     0,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    24,    25,   123,   124,
     125,   126,   271,   206,   127,   128,   129,   130,   131,   132,
       0,     0,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      56,    57,    58,   115,   116,   117,   118,     0,     0,   201,
     135,   129,   130,   131,   132,   115,   116,   117,   118,   129,
     130,   131,   132,   115,   116,   117,   118,   129,   130,   131,
     132,   115,   116,   117,   118,   129,   130,   131,   132,   115,
     116,   117,   118,     0,     0,     0,   272,   129,   130,   131,
     132,     0,     0,     0,   273,     0,     0,     0,   274,     0,
       0,     0,   275,     0,     0,     0,   276,     0,     0,     0,
     277,     0,     0,     0,   278,     0,     0,     0,   279,     0,
       0,     0,   280,   115,   116,   117,   118,     0,     0,     0,
     281,   129,   130,   131,   132,   115,   116,   117,   118,   129,
     130,   131,   132,   115,   116,   117,   118,   129,   130,   131,
     132,   115,   116,   117,   118,   129,   130,   131,   132,   115,
     116,   117,   118,     0,     0,     0,   282,   129,   130,   131,
     132,     0,     0,     0,   283,     0,     0,     0,   284,     0,
       0,     0,   285,     0,     0,     0,   286,     0,     0,     0,
     287,     0,     0,     0,   288,     0,     0,     0,   289,     0,
       0,     0,   290,   115,   116,   117,   118,     0,     0,     0,
     291,   129,   130,   131,   132,   115,   116,   117,   118,   129,
     130,   131,   132,   115,   116,   117,   118,   129,   130,   131,
     132,   115,   116,   117,   118,   129,   130,   131,   132,   115,
     116,   117,   118,     0,     0,     0,   292,   129,   130,   131,
     132,     0,     0,     0,   293,     0,     0,     0,   294,     0,
       0,     0,   295,     0,     0,     0,   296,     0,     0,     0,
     297,     0,     0,     0,   298,     0,     0,     0,   299,     0,
       0,     0,   300,   115,   116,   117,   118,     0,     0,     0,
     301,   129,   130,   131,   132,   115,   116,   117,   118,   129,
     130,   131,   132,   115,   116,   117,   118,   129,   130,   131,
     132,   115,   116,   117,   118,   129,   130,   131,   132,   115,
     116,   117,   118,     0,     0,     0,   302,   129,   130,   131,
     132,     0,     0,     0,   303,     0,     0,     0,   304,     0,
       0,     0,   305,     0,     0,     0,   306,     0,     0,     0,
     307,     0,     0,     0,   308,     0,     0,     0,   309,     0,
       0,     0,   310,   115,   116,   117,   118,     0,     0,     0,
     311,   129,   130,   131,   132,   115,   116,   117,   118,   129,
     130,   131,   132,   115,   116,   117,   118,   129,   130,   131,
     132,   115,   116,   117,   118,   129,   130,   131,   132,   115,
     116,   117,   118,     0,     0,     0,   312,   129,   130,   131,
     132,     0,     0,     0,   313,     0,     0,     0,   314,     0,
       0,     0,   315,     0,     0,     0,   316,     0,     0,     0,
     317,     0,     0,     0,   318,     0,     0,     0,   319,     0,
       0,     0,   320,   115,   116,   117,   118,     0,     0,     0,
     321,   129,   130,   131,   132,   115,   116,   117,   118,   129,
     130,   131,   132,   115,   116,   117,   118,   129,   130,   131,
     132,   123,   124,   125,   126,     0,     0,   127,   128,   129,
     130,   131,   132,     0,     0,     0,   322,     0,     0,     0,
       0,     0,     0,     0,   323,     0,     0,     0,   324,     0,
       0,     0,   325,     0,     0,     0,   199,     0,     0,     0,
     201
};

static const yytype_int16 yycheck[] =
{
      20,    14,    15,    18,    19,    20,    21,    10,    11,    12,
      13,    54,    55,    16,    17,    18,    19,    20,    21,    18,
      19,    20,    21,     0,     5,    63,     3,    53,     8,    14,
      15,    56,     5,    56,     5,    14,    14,    15,    14,    15,
      60,    58,    59,    63,    60,    14,    61,     5,    61,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
      60,    54,    53,    61,    20,    18,    19,    20,    21,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,    18,
      19,    20,    21,   123,   124,   125,   126,   127,   128,   129,
     130,   131,   132,    20,    21,   135,    18,    19,    20,    21,
      20,    21,    20,    21,    60,    60,    60,    63,    61,    18,
      19,    20,    21,    69,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
      96,    97,    98,    99,    60,    60,   204,    60,    60,    60,
      60,    19,    61,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,    60,   204,    60,   206,   123,   124,   125,
     126,   127,   128,   129,   130,   131,   132,    60,    60,   135,
      10,    11,    12,    13,    60,    60,    16,    17,    18,    19,
      20,    21,   107,   108,   107,   108,    18,    19,    20,    21,
      60,    60,    60,   119,   120,   121,   122,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    60,    60,   133,   134,
     133,   134,    60,    60,    60,    60,    60,    10,    -1,    -1,
      -1,    61,    -1,     4,    -1,     6,     7,     8,     9,    61,
      18,    19,    20,    21,    18,    19,    20,    21,   204,    -1,
     206,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,    52,    -1,    61,    -1,    -1,    57,    61,    -1,    60,
      -1,   206,    -1,   206,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,     6,     7,     8,     9,    -1,    -1,    18,    19,
      20,    21,    18,    19,    20,    21,    -1,   233,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,     6,
       7,    61,     9,    -1,    -1,    61,    60,    -1,    -1,    18,
      19,    20,    21,    -1,    -1,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,     6,     7,    10,    11,
      12,    13,    61,    60,    16,    17,    18,    19,    20,    21,
      -1,    -1,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    18,    19,    20,    21,    -1,    -1,    61,
      60,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    18,    19,    20,    21,    18,    19,    20,
      21,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    -1,    -1,    -1,    61,    18,    19,    20,
      21,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    -1,    -1,    -1,    61,    -1,    -1,    -1,
      61,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    18,    19,    20,    21,    -1,    -1,    -1,
      61,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    18,    19,    20,    21,    18,    19,    20,
      21,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    -1,    -1,    -1,    61,    18,    19,    20,
      21,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    -1,    -1,    -1,    61,    -1,    -1,    -1,
      61,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    18,    19,    20,    21,    -1,    -1,    -1,
      61,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    18,    19,    20,    21,    18,    19,    20,
      21,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    -1,    -1,    -1,    61,    18,    19,    20,
      21,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    -1,    -1,    -1,    61,    -1,    -1,    -1,
      61,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    18,    19,    20,    21,    -1,    -1,    -1,
      61,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    18,    19,    20,    21,    18,    19,    20,
      21,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    -1,    -1,    -1,    61,    18,    19,    20,
      21,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    -1,    -1,    -1,    61,    -1,    -1,    -1,
      61,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    18,    19,    20,    21,    -1,    -1,    -1,
      61,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    18,    19,    20,    21,    18,    19,    20,
      21,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    -1,    -1,    -1,    61,    18,    19,    20,
      21,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    -1,    -1,    -1,    61,    -1,    -1,    -1,
      61,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    18,    19,    20,    21,    -1,    -1,    -1,
      61,    18,    19,    20,    21,    18,    19,    20,    21,    18,
      19,    20,    21,    18,    19,    20,    21,    18,    19,    20,
      21,    10,    11,    12,    13,    -1,    -1,    16,    17,    18,
      19,    20,    21,    -1,    -1,    -1,    61,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    61,    -1,    -1,    -1,    61,    -1,
      -1,    -1,    61,    -1,    -1,    -1,    61,    -1,    -1,    -1,
      61
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    63,     0,     3,    64,    69,     5,    53,    65,    66,
      70,     5,    65,    67,    68,    71,    72,    66,    54,    55,
      56,    54,    68,     4,     6,     7,     8,     9,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    57,
      60,    73,    74,    75,    79,    80,    81,    82,    83,    60,
      60,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    60,    60,    60,    60,
      79,    80,    81,    82,    83,    73,    77,    14,    15,    10,
      11,    12,    13,    16,    17,    18,    19,    20,    21,    18,
      19,    20,    21,    10,    11,    12,    13,    16,    17,    18,
      19,    20,    21,    14,    15,    60,    80,    82,    80,    82,
      80,    82,    80,    82,    80,    82,    80,    82,    80,    82,
      80,    82,    80,    82,    80,    82,    80,    82,    80,    82,
      80,    82,    80,    82,    80,    82,    80,    82,    80,    82,
      80,    82,    80,    82,    80,    82,    80,    82,    80,    82,
      80,    82,    80,    82,    80,    82,    80,    82,    80,    82,
      80,    82,    80,    82,    80,    82,    80,    82,    61,    61,
      61,    61,    61,    58,    59,    76,    60,    79,    80,    82,
      83,    79,    83,    80,    82,    80,    82,    80,    82,    80,
      82,    80,    82,    80,    82,    80,    82,    80,    82,    80,
      82,    80,    82,    60,    81,    81,    81,    81,    80,    82,
      80,    82,    80,    82,    80,    82,    80,    82,    80,    82,
      80,    82,    80,    82,    80,    82,    80,    82,    79,    83,
      79,    83,    80,    82,    61,    61,    61,    61,    61,    61,
      61,    61,    61,    61,    61,    61,    61,    61,    61,    61,
      61,    61,    61,    61,    61,    61,    61,    61,    61,    61,
      61,    61,    61,    61,    61,    61,    61,    61,    61,    61,
      61,    61,    61,    61,    61,    61,    61,    61,    61,    61,
      61,    61,    61,    61,    61,    61,    61,    61,    61,    61,
      61,    61,    61,    61,    61,    61,    73,    78
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

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
#ifndef	YYINITDEPTH
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
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
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
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
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

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
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
  int yytoken;
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

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

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
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
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
      if (yyn == 0 || yyn == YYTABLE_NINF)
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
  *++yyvsp = yylval;

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
     `$$ = $1'.

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

/* Line 1455 of yacc.c  */
#line 468 "src/Parameters/parse.y"
    { ;}
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 472 "src/Parameters/parse.y"
    { current_group = ""; 
                                              current_subgroup = "";  ;}
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 474 "src/Parameters/parse.y"
    { current_group = "";
                                              current_subgroup = "";  ;}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 478 "src/Parameters/parse.y"
    { current_subgroup = "";;}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 481 "src/Parameters/parse.y"
    { current_subgroup = ""; ;}
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 482 "src/Parameters/parse.y"
    { current_subgroup = ""; ;}
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 485 "src/Parameters/parse.y"
    {  ;}
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 486 "src/Parameters/parse.y"
    {  ;}
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 489 "src/Parameters/parse.y"
    {  ;}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 490 "src/Parameters/parse.y"
    {  ;}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 493 "src/Parameters/parse.y"
    { current_group = (yyvsp[(1) - (1)].string_type);
                                             current_subgroup = ""; ;}
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 497 "src/Parameters/parse.y"
    { current_subgroup = (yyvsp[(1) - (1)].string_type); ;}
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 500 "src/Parameters/parse.y"
    { current_parameter = (yyvsp[(1) - (1)].string_type);;}
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 503 "src/Parameters/parse.y"
    { new_parameter(); ;}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 507 "src/Parameters/parse.y"
    { current_type = enum_parameter_string;       yylval.string_type = (yyvsp[(1) - (1)].string_type); ;}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 508 "src/Parameters/parse.y"
    { current_type = enum_parameter_integer;      yylval.integer_type = (yyvsp[(1) - (1)].integer_type);;}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 509 "src/Parameters/parse.y"
    { current_type = enum_parameter_scalar;       yylval.scalar_type = (yyvsp[(1) - (1)].scalar_type);;}
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 510 "src/Parameters/parse.y"
    { current_type = enum_parameter_logical;      yylval.logical_type = (yyvsp[(1) - (1)].logical_type); ;}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 511 "src/Parameters/parse.y"
    { current_type = enum_parameter_scalar_expr;  yylval.node_type = (yyvsp[(1) - (1)].node_type); ;}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 512 "src/Parameters/parse.y"
    { current_type = enum_parameter_logical_expr; yylval.node_type = (yyvsp[(1) - (1)].node_type); ;}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 513 "src/Parameters/parse.y"
    { current_type = enum_parameter_list; ;}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 516 "src/Parameters/parse.y"
    {  ;}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 519 "src/Parameters/parse.y"
    { 
   struct param_struct * p = new_param_sentinel();
   p->list_value = param_curr;
   new_param_list(p);
   param_curr = p;
 ;}
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 526 "src/Parameters/parse.y"
    { param_curr = param_curr->list_value; ;}
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 530 "src/Parameters/parse.y"
    { new_parameter(); ;}
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 531 "src/Parameters/parse.y"
    { new_parameter(); ;}
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 533 "src/Parameters/parse.y"
    { ;}
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 538 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(2) - (3)].logical_type); ;}
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 539 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) <= (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 540 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) >= (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 541 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) <  (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 542 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) >  (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 543 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) == (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 544 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) != (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 545 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].logical_type) || (yyvsp[(3) - (3)].logical_type); ;}
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 546 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].logical_type) && (yyvsp[(3) - (3)].logical_type); ;}
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 547 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (1)].logical_type); ;}
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 551 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(2) - (3)].scalar_type); ;}
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 552 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (3)].scalar_type) + (yyvsp[(3) - (3)].scalar_type);;}
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 553 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (3)].scalar_type) - (yyvsp[(3) - (3)].scalar_type);;}
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 554 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (3)].scalar_type) * (yyvsp[(3) - (3)].scalar_type);;}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 555 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (3)].scalar_type) / (yyvsp[(3) - (3)].scalar_type);;}
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 556 "src/Parameters/parse.y"
    { (yyval.scalar_type) = acos((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 557 "src/Parameters/parse.y"
    { (yyval.scalar_type) = acosh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 558 "src/Parameters/parse.y"
    { (yyval.scalar_type) = asin((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 559 "src/Parameters/parse.y"
    { (yyval.scalar_type) = asinh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 560 "src/Parameters/parse.y"
    { (yyval.scalar_type) = atan((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 561 "src/Parameters/parse.y"
    { (yyval.scalar_type) = atanh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 562 "src/Parameters/parse.y"
    { (yyval.scalar_type) = cbrt((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 563 "src/Parameters/parse.y"
    { (yyval.scalar_type) = ceil((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 564 "src/Parameters/parse.y"
    { (yyval.scalar_type) = cos((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 565 "src/Parameters/parse.y"
    { (yyval.scalar_type) = cosh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 566 "src/Parameters/parse.y"
    { (yyval.scalar_type) = erfc((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 567 "src/Parameters/parse.y"
    { (yyval.scalar_type) = erf((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 568 "src/Parameters/parse.y"
    { (yyval.scalar_type) = exp((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 569 "src/Parameters/parse.y"
    { (yyval.scalar_type) = expm1((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 570 "src/Parameters/parse.y"
    { (yyval.scalar_type) = fabs((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 60:

/* Line 1455 of yacc.c  */
#line 571 "src/Parameters/parse.y"
    { (yyval.scalar_type) = floor((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 573 "src/Parameters/parse.y"
    { (yyval.scalar_type) = j0((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 574 "src/Parameters/parse.y"
    { (yyval.scalar_type) = j1((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 575 "src/Parameters/parse.y"
    { (yyval.scalar_type) = lgamma((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 576 "src/Parameters/parse.y"
    { (yyval.scalar_type) = log10((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 577 "src/Parameters/parse.y"
    { (yyval.scalar_type) = log1p((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 66:

/* Line 1455 of yacc.c  */
#line 578 "src/Parameters/parse.y"
    { (yyval.scalar_type) = logb((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 67:

/* Line 1455 of yacc.c  */
#line 579 "src/Parameters/parse.y"
    { (yyval.scalar_type) = log((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 68:

/* Line 1455 of yacc.c  */
#line 580 "src/Parameters/parse.y"
    { (yyval.scalar_type) = sin((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 69:

/* Line 1455 of yacc.c  */
#line 581 "src/Parameters/parse.y"
    { (yyval.scalar_type) = sinh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 70:

/* Line 1455 of yacc.c  */
#line 582 "src/Parameters/parse.y"
    { (yyval.scalar_type) = sqrt((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 71:

/* Line 1455 of yacc.c  */
#line 583 "src/Parameters/parse.y"
    { (yyval.scalar_type) = tan((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 72:

/* Line 1455 of yacc.c  */
#line 584 "src/Parameters/parse.y"
    { (yyval.scalar_type) = tanh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 73:

/* Line 1455 of yacc.c  */
#line 585 "src/Parameters/parse.y"
    { (yyval.scalar_type) = y0((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 74:

/* Line 1455 of yacc.c  */
#line 586 "src/Parameters/parse.y"
    { (yyval.scalar_type) = y1((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 587 "src/Parameters/parse.y"
    { (yyval.scalar_type) = rint((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 76:

/* Line 1455 of yacc.c  */
#line 588 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (1)].scalar_type);;}
    break;

  case 77:

/* Line 1455 of yacc.c  */
#line 592 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(2) - (3)].integer_type); ;}
    break;

  case 78:

/* Line 1455 of yacc.c  */
#line 593 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (3)].integer_type) + (yyvsp[(3) - (3)].integer_type);;}
    break;

  case 79:

/* Line 1455 of yacc.c  */
#line 594 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (3)].integer_type) - (yyvsp[(3) - (3)].integer_type);;}
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 595 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (3)].integer_type) * (yyvsp[(3) - (3)].integer_type);;}
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 596 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (3)].integer_type) / (yyvsp[(3) - (3)].integer_type);;}
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 597 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (1)].integer_type);;}
    break;

  case 83:

/* Line 1455 of yacc.c  */
#line 601 "src/Parameters/parse.y"
    { (yyval.node_type) = (yyvsp[(2) - (3)].node_type); ;}
    break;

  case 84:

/* Line 1455 of yacc.c  */
#line 602 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_add,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 603 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_add,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 604 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_add,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 605 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_sub,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 606 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_sub,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 89:

/* Line 1455 of yacc.c  */
#line 607 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_sub,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 90:

/* Line 1455 of yacc.c  */
#line 608 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_mul,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 91:

/* Line 1455 of yacc.c  */
#line 609 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_mul,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 92:

/* Line 1455 of yacc.c  */
#line 610 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_mul,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 93:

/* Line 1455 of yacc.c  */
#line 611 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_div,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 94:

/* Line 1455 of yacc.c  */
#line 612 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_div,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 95:

/* Line 1455 of yacc.c  */
#line 613 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_div,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 96:

/* Line 1455 of yacc.c  */
#line 614 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( acos, "acos", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 97:

/* Line 1455 of yacc.c  */
#line 615 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( acosh, "acosh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 616 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( asin, "asin", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 617 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( asinh, "asinh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 100:

/* Line 1455 of yacc.c  */
#line 618 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( atan, "atan", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 101:

/* Line 1455 of yacc.c  */
#line 619 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( atanh, "atanh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 620 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( cbrt, "cbrt", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 103:

/* Line 1455 of yacc.c  */
#line 621 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( ceil, "ceil", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 622 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( cos, "cos", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 105:

/* Line 1455 of yacc.c  */
#line 623 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( cosh, "cosh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 106:

/* Line 1455 of yacc.c  */
#line 624 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( erfc, "erfc", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 107:

/* Line 1455 of yacc.c  */
#line 625 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( erf, "erf", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 108:

/* Line 1455 of yacc.c  */
#line 626 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( exp, "exp", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 109:

/* Line 1455 of yacc.c  */
#line 627 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( expm1, "expm1", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 110:

/* Line 1455 of yacc.c  */
#line 628 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( fabs, "fabs", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 111:

/* Line 1455 of yacc.c  */
#line 629 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( floor, "floor", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 112:

/* Line 1455 of yacc.c  */
#line 631 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( j0, "j0", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 113:

/* Line 1455 of yacc.c  */
#line 632 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( j1, "j1", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 114:

/* Line 1455 of yacc.c  */
#line 633 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( lgamma, "lgamma", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 115:

/* Line 1455 of yacc.c  */
#line 634 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( log10, "log10", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 635 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( log1p, "log1p", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 117:

/* Line 1455 of yacc.c  */
#line 636 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( logb, "logb", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 118:

/* Line 1455 of yacc.c  */
#line 637 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( log, "log", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 638 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( sin, "sin", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 639 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( sinh, "sinh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 121:

/* Line 1455 of yacc.c  */
#line 640 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( sqrt, "sqrt", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 641 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( tan, "tan", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 642 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( tanh, "tanh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 124:

/* Line 1455 of yacc.c  */
#line 643 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( y0, "y0", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 125:

/* Line 1455 of yacc.c  */
#line 644 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( y1, "y1", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 126:

/* Line 1455 of yacc.c  */
#line 645 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( rint, "rint", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 127:

/* Line 1455 of yacc.c  */
#line 646 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_variable ((yyvsp[(1) - (1)].string_type));  ;}
    break;

  case 128:

/* Line 1455 of yacc.c  */
#line 651 "src/Parameters/parse.y"
    { ;}
    break;

  case 129:

/* Line 1455 of yacc.c  */
#line 652 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_le,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 653 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_le,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 654 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_le,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 132:

/* Line 1455 of yacc.c  */
#line 655 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_ge,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 656 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_ge,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 134:

/* Line 1455 of yacc.c  */
#line 657 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_ge,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 135:

/* Line 1455 of yacc.c  */
#line 658 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_lt,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 136:

/* Line 1455 of yacc.c  */
#line 659 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_lt,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 137:

/* Line 1455 of yacc.c  */
#line 660 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_lt,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 138:

/* Line 1455 of yacc.c  */
#line 661 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_gt,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 139:

/* Line 1455 of yacc.c  */
#line 662 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_gt,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 140:

/* Line 1455 of yacc.c  */
#line 663 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_gt,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 141:

/* Line 1455 of yacc.c  */
#line 664 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_eq,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 142:

/* Line 1455 of yacc.c  */
#line 665 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_eq,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 143:

/* Line 1455 of yacc.c  */
#line 666 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_eq,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 144:

/* Line 1455 of yacc.c  */
#line 667 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_ne,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 145:

/* Line 1455 of yacc.c  */
#line 668 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_ne,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 146:

/* Line 1455 of yacc.c  */
#line 669 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_ne,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 147:

/* Line 1455 of yacc.c  */
#line 670 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_or,new_node_logical((yyvsp[(3) - (3)].logical_type))); ;}
    break;

  case 148:

/* Line 1455 of yacc.c  */
#line 671 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_logical((yyvsp[(1) - (3)].logical_type)), enum_op_or,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 149:

/* Line 1455 of yacc.c  */
#line 672 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_or,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 150:

/* Line 1455 of yacc.c  */
#line 673 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_and,new_node_logical((yyvsp[(3) - (3)].logical_type))); ;}
    break;

  case 151:

/* Line 1455 of yacc.c  */
#line 674 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_logical((yyvsp[(1) - (3)].logical_type)), enum_op_and,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 152:

/* Line 1455 of yacc.c  */
#line 675 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_and,(yyvsp[(3) - (3)].node_type)); ;}
    break;



/* Line 1455 of yacc.c  */
#line 3229 "src/Parameters/parse.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
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

  /* Do not reclaim the symbols of the rule which action triggered
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
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
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

  *++yyvsp = yylval;


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

#if !defined(yyoverflow) || YYERROR_VERBOSE
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
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
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
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 680 "src/Parameters/parse.y"


struct param_struct * 
cello_parameters_read(FILE * fp)
{
  /* initialize the linked list with an initial sentinel (sentinel) node */
  param_head = param_curr = new_param_sentinel();

  /*   yydebug=1; */
  
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
    case enum_node_scalar:
      fprintf (fp,"%g",node->scalar_value);
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
    case enum_node_scalar:
      sprintf (buffer,"%g",node->scalar_value);
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

  while (p && p->type != enum_parameter_sentinel) {

    if (p->group != NULL) {
      indent(level);
      printf ("%s %s:%s:%s = ", 
	      parameter_name[p->type],p->group, p->subgroup, p->parameter);
    } else {
      /* list element */
      indent(level);
      printf ("%s %s = ", 
	      parameter_name[p->type], p->parameter);
    }
    switch (p->type) {
    case enum_parameter_scalar:  
      printf ("%g\n",p->scalar_value);  
      break;
    case enum_parameter_integer: 
      printf ("%d\n",p->integer_value); 
      break;
    case enum_parameter_string:  
      printf ("%s\n",p->string_value); 
      break;
    case enum_parameter_subgroup:  
      printf ("Uh oh: SUBGROUP %s (should be deleted)\n",p->string_value);
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
    case enum_parameter_scalar_expr:
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


