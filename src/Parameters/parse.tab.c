
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
    "group",
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
     
     copy_groups(p->group,current_group);

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
#line 479 "src/Parameters/parse.tab.c"

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
     STRING = 258,
     IDENTIFIER = 259,
     VARIABLE = 260,
     SCALAR = 261,
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



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 407 "src/Parameters/parse.y"
 
  int logical_type;  
  int integer_type; 
  double scalar_type;  
  char * string_type; 
  char * group_type;
  struct node_expr * node_type;
  


/* Line 214 of yacc.c  */
#line 571 "src/Parameters/parse.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 583 "src/Parameters/parse.tab.c"

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
#define YYLAST   969

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  61
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  19
/* YYNRULES -- Number of rules.  */
#define YYNRULES  150
/* YYNRULES -- Number of states.  */
#define YYNSTATES  325

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   300

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
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
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,    10,    14,    19,    21,    25,
      27,    31,    33,    35,    39,    41,    43,    45,    47,    49,
      51,    53,    57,    60,    62,    64,    66,    67,    72,    76,
      80,    84,    88,    92,    96,   100,   104,   108,   110,   114,
     118,   122,   126,   130,   135,   140,   145,   150,   155,   160,
     165,   170,   175,   180,   185,   190,   195,   200,   205,   210,
     215,   220,   225,   230,   235,   240,   245,   250,   255,   260,
     265,   270,   275,   280,   285,   287,   291,   295,   299,   303,
     307,   309,   313,   317,   321,   325,   329,   333,   337,   341,
     345,   349,   353,   357,   361,   366,   371,   376,   381,   386,
     391,   396,   401,   406,   411,   416,   421,   426,   431,   436,
     441,   446,   451,   456,   461,   466,   471,   476,   481,   486,
     491,   496,   501,   506,   511,   516,   518,   522,   526,   530,
     534,   538,   542,   546,   550,   554,   558,   562,   566,   570,
     574,   578,   582,   586,   590,   594,   598,   602,   606,   610,
     614
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      62,     0,    -1,    -1,    62,    63,    -1,    66,    64,    -1,
      52,    65,    53,    -1,    52,    65,    54,    53,    -1,    68,
      -1,    65,    54,    68,    -1,    63,    -1,    65,    54,    63,
      -1,     4,    -1,     4,    -1,    67,    55,    69,    -1,     3,
      -1,    77,    -1,    76,    -1,    75,    -1,    78,    -1,    79,
      -1,    70,    -1,    71,    73,    72,    -1,    71,    72,    -1,
      56,    -1,    57,    -1,    69,    -1,    -1,    73,    58,    69,
      74,    -1,    59,    75,    60,    -1,    76,     9,    76,    -1,
      76,    10,    76,    -1,    76,    15,    76,    -1,    76,    16,
      76,    -1,    76,    12,    76,    -1,    76,    11,    76,    -1,
      75,    14,    75,    -1,    75,    13,    75,    -1,     8,    -1,
      59,    76,    60,    -1,    76,    17,    76,    -1,    76,    18,
      76,    -1,    76,    19,    76,    -1,    76,    20,    76,    -1,
      21,    59,    76,    60,    -1,    22,    59,    76,    60,    -1,
      23,    59,    76,    60,    -1,    24,    59,    76,    60,    -1,
      25,    59,    76,    60,    -1,    26,    59,    76,    60,    -1,
      27,    59,    76,    60,    -1,    28,    59,    76,    60,    -1,
      29,    59,    76,    60,    -1,    30,    59,    76,    60,    -1,
      31,    59,    76,    60,    -1,    32,    59,    76,    60,    -1,
      33,    59,    76,    60,    -1,    34,    59,    76,    60,    -1,
      35,    59,    76,    60,    -1,    36,    59,    76,    60,    -1,
      37,    59,    76,    60,    -1,    38,    59,    76,    60,    -1,
      39,    59,    76,    60,    -1,    40,    59,    76,    60,    -1,
      41,    59,    76,    60,    -1,    42,    59,    76,    60,    -1,
      43,    59,    76,    60,    -1,    44,    59,    76,    60,    -1,
      45,    59,    76,    60,    -1,    46,    59,    76,    60,    -1,
      47,    59,    76,    60,    -1,    48,    59,    76,    60,    -1,
      49,    59,    76,    60,    -1,    50,    59,    76,    60,    -1,
      51,    59,    76,    60,    -1,     6,    -1,    59,    77,    60,
      -1,    77,    17,    77,    -1,    77,    18,    77,    -1,    77,
      19,    77,    -1,    77,    20,    77,    -1,     7,    -1,    59,
      78,    60,    -1,    78,    17,    76,    -1,    76,    17,    78,
      -1,    78,    17,    78,    -1,    78,    18,    76,    -1,    76,
      18,    78,    -1,    78,    18,    78,    -1,    78,    19,    76,
      -1,    76,    19,    78,    -1,    78,    19,    78,    -1,    78,
      20,    76,    -1,    76,    20,    78,    -1,    78,    20,    78,
      -1,    21,    59,    78,    60,    -1,    22,    59,    78,    60,
      -1,    23,    59,    78,    60,    -1,    24,    59,    78,    60,
      -1,    25,    59,    78,    60,    -1,    26,    59,    78,    60,
      -1,    27,    59,    78,    60,    -1,    28,    59,    78,    60,
      -1,    29,    59,    78,    60,    -1,    30,    59,    78,    60,
      -1,    31,    59,    78,    60,    -1,    32,    59,    78,    60,
      -1,    33,    59,    78,    60,    -1,    34,    59,    78,    60,
      -1,    35,    59,    78,    60,    -1,    36,    59,    78,    60,
      -1,    37,    59,    78,    60,    -1,    38,    59,    78,    60,
      -1,    39,    59,    78,    60,    -1,    40,    59,    78,    60,
      -1,    41,    59,    78,    60,    -1,    42,    59,    78,    60,
      -1,    43,    59,    78,    60,    -1,    44,    59,    78,    60,
      -1,    45,    59,    78,    60,    -1,    46,    59,    78,    60,
      -1,    47,    59,    78,    60,    -1,    48,    59,    78,    60,
      -1,    49,    59,    78,    60,    -1,    50,    59,    78,    60,
      -1,    51,    59,    78,    60,    -1,     5,    -1,    59,    79,
      60,    -1,    78,     9,    76,    -1,    76,     9,    78,    -1,
      78,     9,    78,    -1,    78,    10,    76,    -1,    76,    10,
      78,    -1,    78,    10,    78,    -1,    78,    15,    76,    -1,
      76,    15,    78,    -1,    78,    15,    78,    -1,    78,    16,
      76,    -1,    76,    16,    78,    -1,    78,    16,    78,    -1,
      78,    12,    76,    -1,    76,    12,    78,    -1,    78,    12,
      78,    -1,    78,    11,    76,    -1,    76,    11,    78,    -1,
      78,    11,    78,    -1,    79,    14,    75,    -1,    75,    14,
      79,    -1,    79,    14,    79,    -1,    79,    13,    75,    -1,
      75,    13,    79,    -1,    79,    13,    79,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   486,   486,   487,   491,   494,   495,   498,   499,   500,
     501,   505,   508,   511,   515,   516,   517,   518,   519,   520,
     521,   524,   525,   528,   535,   539,   540,   540,   547,   548,
     549,   550,   551,   552,   553,   554,   555,   556,   560,   561,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   582,
     583,   584,   585,   586,   587,   588,   589,   590,   591,   592,
     593,   594,   595,   596,   597,   601,   602,   603,   604,   605,
     606,   610,   611,   612,   613,   614,   615,   616,   617,   618,
     619,   620,   621,   622,   623,   624,   625,   626,   627,   628,
     629,   630,   631,   632,   633,   634,   635,   636,   637,   638,
     640,   641,   642,   643,   644,   645,   646,   647,   648,   649,
     650,   651,   652,   653,   654,   655,   660,   661,   662,   663,
     664,   665,   666,   667,   668,   669,   670,   671,   672,   673,
     674,   675,   676,   677,   678,   679,   680,   681,   682,   683,
     684
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "STRING", "IDENTIFIER", "VARIABLE",
  "SCALAR", "INTEGER", "LOGICAL", "LE", "GE", "NE", "EQ", "AND", "OR",
  "'<'", "'>'", "'+'", "'-'", "'*'", "'/'", "ACOS", "ACOSH", "ASIN",
  "ASINH", "ATAN", "ATANH", "CBRT", "CEIL", "COS", "COSH", "ERFC", "ERF",
  "EXP", "EXPM1", "FABS", "FLOOR", "J0", "J1", "LGAMMA", "LOG10", "LOG1P",
  "LOGB", "LOG", "SIN", "SINH", "SQRT", "TAN", "TANH", "Y0", "Y1", "RINT",
  "'{'", "'}'", "';'", "'='", "'['", "']'", "','", "'('", "')'", "$accept",
  "file", "group", "parameter_group", "parameter_list", "group_name",
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
     265,   266,   267,   268,   269,    60,    62,    43,    45,    42,
      47,   270,   271,   272,   273,   274,   275,   276,   277,   278,
     279,   280,   281,   282,   283,   284,   285,   286,   287,   288,
     289,   290,   291,   292,   293,   294,   295,   296,   297,   298,
     299,   300,   123,   125,    59,    61,    91,    93,    44,    40,
      41
};
# endif

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

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
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

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
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

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,     4,     7,    10,     5,    11,    12,    57,    58,
      59,   103,   104,   324,    60,   205,    62,   206,    64
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -41
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

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -41,   -41,     5,   -41,   -41,   -41,   -41,    33,   -36,   -41,
     -41,   140,   -41,   -41,   108,   -15,   -31,    68,   110
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -13
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
#line 487 "src/Parameters/parse.y"
    { ;}
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 491 "src/Parameters/parse.y"
    {  ;}
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 494 "src/Parameters/parse.y"
    { current_group[--current_group_level] = 0; ;}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 495 "src/Parameters/parse.y"
    { current_group[--current_group_level] = 0; ;}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 498 "src/Parameters/parse.y"
    {  ;}
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 499 "src/Parameters/parse.y"
    {  ;}
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 500 "src/Parameters/parse.y"
    {  ;}
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 501 "src/Parameters/parse.y"
    {  ;}
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 505 "src/Parameters/parse.y"
    { current_group[current_group_level++] = (yyvsp[(1) - (1)].string_type); ;}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 508 "src/Parameters/parse.y"
    { current_parameter = (yyvsp[(1) - (1)].string_type);;}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 511 "src/Parameters/parse.y"
    { new_parameter(); ;}
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 515 "src/Parameters/parse.y"
    { current_type = enum_parameter_string;       yylval.string_type = (yyvsp[(1) - (1)].string_type); ;}
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 516 "src/Parameters/parse.y"
    { current_type = enum_parameter_integer;      yylval.integer_type = (yyvsp[(1) - (1)].integer_type);;}
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 517 "src/Parameters/parse.y"
    { current_type = enum_parameter_scalar;       yylval.scalar_type = (yyvsp[(1) - (1)].scalar_type);;}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 518 "src/Parameters/parse.y"
    { current_type = enum_parameter_logical;      yylval.logical_type = (yyvsp[(1) - (1)].logical_type); ;}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 519 "src/Parameters/parse.y"
    { current_type = enum_parameter_scalar_expr;  yylval.node_type = (yyvsp[(1) - (1)].node_type); ;}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 520 "src/Parameters/parse.y"
    { current_type = enum_parameter_logical_expr; yylval.node_type = (yyvsp[(1) - (1)].node_type); ;}
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 521 "src/Parameters/parse.y"
    { current_type = enum_parameter_list; ;}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 524 "src/Parameters/parse.y"
    {  ;}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 525 "src/Parameters/parse.y"
    {  ;}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 528 "src/Parameters/parse.y"
    { 
   struct param_struct * p = new_param_sentinel();
   p->list_value = param_curr;
   new_param_list(p);
   param_curr = p;
 ;}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 535 "src/Parameters/parse.y"
    { param_curr = param_curr->list_value; ;}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 539 "src/Parameters/parse.y"
    { new_parameter(); ;}
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 540 "src/Parameters/parse.y"
    { new_parameter(); ;}
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 542 "src/Parameters/parse.y"
    { ;}
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 547 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(2) - (3)].logical_type); ;}
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 548 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) <= (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 549 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) >= (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 550 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) <  (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 551 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) >  (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 552 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) == (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 553 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].scalar_type) != (yyvsp[(3) - (3)].scalar_type); ;}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 554 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].logical_type) || (yyvsp[(3) - (3)].logical_type); ;}
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 555 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (3)].logical_type) && (yyvsp[(3) - (3)].logical_type); ;}
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 556 "src/Parameters/parse.y"
    { (yyval.logical_type) = (yyvsp[(1) - (1)].logical_type); ;}
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 560 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(2) - (3)].scalar_type); ;}
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 561 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (3)].scalar_type) + (yyvsp[(3) - (3)].scalar_type);;}
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 562 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (3)].scalar_type) - (yyvsp[(3) - (3)].scalar_type);;}
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 563 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (3)].scalar_type) * (yyvsp[(3) - (3)].scalar_type);;}
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 564 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (3)].scalar_type) / (yyvsp[(3) - (3)].scalar_type);;}
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 565 "src/Parameters/parse.y"
    { (yyval.scalar_type) = acos((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 566 "src/Parameters/parse.y"
    { (yyval.scalar_type) = acosh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 567 "src/Parameters/parse.y"
    { (yyval.scalar_type) = asin((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 568 "src/Parameters/parse.y"
    { (yyval.scalar_type) = asinh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 569 "src/Parameters/parse.y"
    { (yyval.scalar_type) = atan((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 570 "src/Parameters/parse.y"
    { (yyval.scalar_type) = atanh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 571 "src/Parameters/parse.y"
    { (yyval.scalar_type) = cbrt((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 572 "src/Parameters/parse.y"
    { (yyval.scalar_type) = ceil((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 573 "src/Parameters/parse.y"
    { (yyval.scalar_type) = cos((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 574 "src/Parameters/parse.y"
    { (yyval.scalar_type) = cosh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 575 "src/Parameters/parse.y"
    { (yyval.scalar_type) = erfc((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 576 "src/Parameters/parse.y"
    { (yyval.scalar_type) = erf((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 577 "src/Parameters/parse.y"
    { (yyval.scalar_type) = exp((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 578 "src/Parameters/parse.y"
    { (yyval.scalar_type) = expm1((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 579 "src/Parameters/parse.y"
    { (yyval.scalar_type) = fabs((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 580 "src/Parameters/parse.y"
    { (yyval.scalar_type) = floor((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 582 "src/Parameters/parse.y"
    { (yyval.scalar_type) = j0((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 60:

/* Line 1455 of yacc.c  */
#line 583 "src/Parameters/parse.y"
    { (yyval.scalar_type) = j1((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 584 "src/Parameters/parse.y"
    { (yyval.scalar_type) = lgamma((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 585 "src/Parameters/parse.y"
    { (yyval.scalar_type) = log10((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 586 "src/Parameters/parse.y"
    { (yyval.scalar_type) = log1p((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 587 "src/Parameters/parse.y"
    { (yyval.scalar_type) = logb((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 588 "src/Parameters/parse.y"
    { (yyval.scalar_type) = log((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 66:

/* Line 1455 of yacc.c  */
#line 589 "src/Parameters/parse.y"
    { (yyval.scalar_type) = sin((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 67:

/* Line 1455 of yacc.c  */
#line 590 "src/Parameters/parse.y"
    { (yyval.scalar_type) = sinh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 68:

/* Line 1455 of yacc.c  */
#line 591 "src/Parameters/parse.y"
    { (yyval.scalar_type) = sqrt((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 69:

/* Line 1455 of yacc.c  */
#line 592 "src/Parameters/parse.y"
    { (yyval.scalar_type) = tan((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 70:

/* Line 1455 of yacc.c  */
#line 593 "src/Parameters/parse.y"
    { (yyval.scalar_type) = tanh((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 71:

/* Line 1455 of yacc.c  */
#line 594 "src/Parameters/parse.y"
    { (yyval.scalar_type) = y0((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 72:

/* Line 1455 of yacc.c  */
#line 595 "src/Parameters/parse.y"
    { (yyval.scalar_type) = y1((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 73:

/* Line 1455 of yacc.c  */
#line 596 "src/Parameters/parse.y"
    { (yyval.scalar_type) = rint((yyvsp[(3) - (4)].scalar_type)); ;}
    break;

  case 74:

/* Line 1455 of yacc.c  */
#line 597 "src/Parameters/parse.y"
    { (yyval.scalar_type) = (yyvsp[(1) - (1)].scalar_type);;}
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 601 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(2) - (3)].integer_type); ;}
    break;

  case 76:

/* Line 1455 of yacc.c  */
#line 602 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (3)].integer_type) + (yyvsp[(3) - (3)].integer_type);;}
    break;

  case 77:

/* Line 1455 of yacc.c  */
#line 603 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (3)].integer_type) - (yyvsp[(3) - (3)].integer_type);;}
    break;

  case 78:

/* Line 1455 of yacc.c  */
#line 604 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (3)].integer_type) * (yyvsp[(3) - (3)].integer_type);;}
    break;

  case 79:

/* Line 1455 of yacc.c  */
#line 605 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (3)].integer_type) / (yyvsp[(3) - (3)].integer_type);;}
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 606 "src/Parameters/parse.y"
    { (yyval.integer_type) = (yyvsp[(1) - (1)].integer_type);;}
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 610 "src/Parameters/parse.y"
    { (yyval.node_type) = (yyvsp[(2) - (3)].node_type); ;}
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 611 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_add,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 83:

/* Line 1455 of yacc.c  */
#line 612 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_add,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 84:

/* Line 1455 of yacc.c  */
#line 613 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_add,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 614 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_sub,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 615 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_sub,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 616 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_sub,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 617 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_mul,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 89:

/* Line 1455 of yacc.c  */
#line 618 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_mul,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 90:

/* Line 1455 of yacc.c  */
#line 619 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_mul,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 91:

/* Line 1455 of yacc.c  */
#line 620 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_div,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 92:

/* Line 1455 of yacc.c  */
#line 621 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_div,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 93:

/* Line 1455 of yacc.c  */
#line 622 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_div,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 94:

/* Line 1455 of yacc.c  */
#line 623 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( acos, "acos", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 95:

/* Line 1455 of yacc.c  */
#line 624 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( acosh, "acosh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 96:

/* Line 1455 of yacc.c  */
#line 625 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( asin, "asin", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 97:

/* Line 1455 of yacc.c  */
#line 626 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( asinh, "asinh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 627 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( atan, "atan", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 628 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( atanh, "atanh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 100:

/* Line 1455 of yacc.c  */
#line 629 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( cbrt, "cbrt", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 101:

/* Line 1455 of yacc.c  */
#line 630 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( ceil, "ceil", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 631 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( cos, "cos", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 103:

/* Line 1455 of yacc.c  */
#line 632 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( cosh, "cosh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 633 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( erfc, "erfc", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 105:

/* Line 1455 of yacc.c  */
#line 634 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( erf, "erf", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 106:

/* Line 1455 of yacc.c  */
#line 635 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( exp, "exp", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 107:

/* Line 1455 of yacc.c  */
#line 636 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( expm1, "expm1", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 108:

/* Line 1455 of yacc.c  */
#line 637 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( fabs, "fabs", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 109:

/* Line 1455 of yacc.c  */
#line 638 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( floor, "floor", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 110:

/* Line 1455 of yacc.c  */
#line 640 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( j0, "j0", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 111:

/* Line 1455 of yacc.c  */
#line 641 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( j1, "j1", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 112:

/* Line 1455 of yacc.c  */
#line 642 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( lgamma, "lgamma", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 113:

/* Line 1455 of yacc.c  */
#line 643 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( log10, "log10", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 114:

/* Line 1455 of yacc.c  */
#line 644 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( log1p, "log1p", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 115:

/* Line 1455 of yacc.c  */
#line 645 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( logb, "logb", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 646 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( log, "log", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 117:

/* Line 1455 of yacc.c  */
#line 647 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( sin, "sin", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 118:

/* Line 1455 of yacc.c  */
#line 648 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( sinh, "sinh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 649 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( sqrt, "sqrt", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 650 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( tan, "tan", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 121:

/* Line 1455 of yacc.c  */
#line 651 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( tanh, "tanh", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 652 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( y0, "y0", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 653 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( y1, "y1", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 124:

/* Line 1455 of yacc.c  */
#line 654 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_function ( rint, "rint", (yyvsp[(3) - (4)].node_type)); ;}
    break;

  case 125:

/* Line 1455 of yacc.c  */
#line 655 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_variable ((yyvsp[(1) - (1)].string_type));  ;}
    break;

  case 126:

/* Line 1455 of yacc.c  */
#line 660 "src/Parameters/parse.y"
    { ;}
    break;

  case 127:

/* Line 1455 of yacc.c  */
#line 661 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_le,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 128:

/* Line 1455 of yacc.c  */
#line 662 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_le,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 129:

/* Line 1455 of yacc.c  */
#line 663 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_le,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 664 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_ge,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 665 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_ge,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 132:

/* Line 1455 of yacc.c  */
#line 666 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_ge,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 667 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_lt,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 134:

/* Line 1455 of yacc.c  */
#line 668 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_lt,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 135:

/* Line 1455 of yacc.c  */
#line 669 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_lt,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 136:

/* Line 1455 of yacc.c  */
#line 670 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_gt,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 137:

/* Line 1455 of yacc.c  */
#line 671 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_gt,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 138:

/* Line 1455 of yacc.c  */
#line 672 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_gt,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 139:

/* Line 1455 of yacc.c  */
#line 673 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_eq,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 140:

/* Line 1455 of yacc.c  */
#line 674 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_eq,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 141:

/* Line 1455 of yacc.c  */
#line 675 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_eq,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 142:

/* Line 1455 of yacc.c  */
#line 676 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_ne,new_node_scalar((yyvsp[(3) - (3)].scalar_type))); ;}
    break;

  case 143:

/* Line 1455 of yacc.c  */
#line 677 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_scalar((yyvsp[(1) - (3)].scalar_type)), enum_op_ne,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 144:

/* Line 1455 of yacc.c  */
#line 678 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_ne,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 145:

/* Line 1455 of yacc.c  */
#line 679 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_or,new_node_logical((yyvsp[(3) - (3)].logical_type))); ;}
    break;

  case 146:

/* Line 1455 of yacc.c  */
#line 680 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_logical((yyvsp[(1) - (3)].logical_type)), enum_op_or,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 147:

/* Line 1455 of yacc.c  */
#line 681 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_or,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 148:

/* Line 1455 of yacc.c  */
#line 682 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_and,new_node_logical((yyvsp[(3) - (3)].logical_type))); ;}
    break;

  case 149:

/* Line 1455 of yacc.c  */
#line 683 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation (new_node_logical((yyvsp[(1) - (3)].logical_type)), enum_op_and,(yyvsp[(3) - (3)].node_type)); ;}
    break;

  case 150:

/* Line 1455 of yacc.c  */
#line 684 "src/Parameters/parse.y"
    { (yyval.node_type) = new_node_operation ((yyvsp[(1) - (3)].node_type), enum_op_and,(yyvsp[(3) - (3)].node_type)); ;}
    break;



/* Line 1455 of yacc.c  */
#line 3232 "src/Parameters/parse.tab.c"
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
#line 689 "src/Parameters/parse.y"


struct param_struct * 
cello_parameters_read(FILE * fp)
{
  clear_groups(current_group);

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
    case enum_parameter_scalar:  
      printf ("%g\n",p->scalar_value);  
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


