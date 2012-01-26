/* ANY CHANGES HERE MUST BE REFLECTED IN parse.y parameter_name[] */

#define MAX_GROUP_DEPTH 10

enum enum_parameter {
  enum_parameter_unknown,
  enum_parameter_sentinel,
  enum_parameter_group,
  enum_parameter_integer,
  enum_parameter_float,
  enum_parameter_string,
  enum_parameter_identifier,
  enum_parameter_logical,
  enum_parameter_list,
  enum_parameter_float_expr,
  enum_parameter_logical_expr,
  enum_parameter_function
};


enum enum_node {
  enum_node_unknown,
  enum_node_operation,
  enum_node_float,
  enum_node_integer,
  enum_node_variable,
  enum_node_function
};

enum enum_op {
  enum_op_add,
  enum_op_sub,
  enum_op_mul,
  enum_op_div,
  enum_op_le,
  enum_op_lt,
  enum_op_ge,
  enum_op_gt,
  enum_op_eq,
  enum_op_ne,
  enum_op_and,
  enum_op_or
};

struct node_expr {
  enum enum_node type;
  union {
    enum enum_op op_value;       /* arthmetic / logical operation */
    double       float_value;   /* floating point number */
    int          integer_value;  /* integer / logical constant */
    char         var_value;      /* variable, e.g. x,y,z,t */
    double (*fun_value)(double); /* math.h function */
  };
  struct node_expr * left;
  struct node_expr * right;
  char * function_name;
};

struct param_struct {
  char * group[MAX_GROUP_DEPTH];
  char * parameter;
  enum enum_parameter type;
  union  {
    int                 logical_value; 
    int                 integer_value; 
    double              float_value; 
    char *              string_value;
    struct param_struct * list_value;
    struct node_expr    * op_value;    /* expression tree */
  };
  struct param_struct *   next;
};

