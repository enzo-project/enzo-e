  enum node_type {
    type_node_unknown,
    type_node_operation,
    type_node_scalar,
    type_node_integer,
    type_node_variable,
    type_node_function
  };

  enum op_type {
    type_op_add,
    type_op_sub,
    type_op_mul,
    type_op_div,
    type_op_le,
    type_op_lt,
    type_op_ge,
    type_op_gt,
    type_op_eq,
    type_op_ne,
    type_op_and,
    type_op_or
  };

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

