#ifndef PARAM_HPP
#define PARAM_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      param.hpp
 * @brief     Helper Param* classes for parameters
 * @author    James Bordner
 * @date      Sun Oct 11 14:55:25 PDT 2009
 * @bug       
 * @note      
 *
 * Classes Param* to represent and evaluate various parameter types and
 * expressions
 *
 * $Id$
 *
 *********************************************************************
 */


//======================================================================

class Param {

  friend class Parameters;

public:

  void evaluate_scalar  
  ( int n, double * result, double *x, double *y, double *z, double *t);

  void dealloc() { 
    switch (type_) {
    case type_string_: 
      dealloc_string_(); 
      break;
    case type_list_:   
      dealloc_param_(value_list_); 
      break;
    case type_logical_expr_:
    case type_scalar_expr_:
      dealloc_string_();
      break;
    }
  } 

private:

  enum type_param {
    type_unknown_,
    type_integer_,
    type_scalar_,
    type_logical_,
    type_string_,
    type_list_,
    type_scalar_expr_,
    type_logical_expr_
  };

  enum type_param type_;

  union {
    int                value_integer_; 
    double             value_scalar_; 
    bool               value_logical_; 
    char *             value_string_;
    struct param_type * value_list_;
    struct node_expr * value_expr_;
  };

  void set_integer_ (int value)    
  { 
    type_ = type_integer_; 
    value_integer_ = value; 
  };

  void set_scalar_  (double value) 
  { 
    type_ = type_scalar_; 
    value_scalar_ = value; 
  };

  void set_logical_ (int value)    
  { 
    type_ = type_logical_; 
    value_logical_ = (value != 0); 
  };

  void set_string_ (char * value) 
  { 
    type_ = type_string_; 
    value_string_ = value; 
  };

  void set_list_ (struct param_type * value) 
  { 
    type_ = type_list_; 
    value_list_ = value; 
  };

  void set_scalar_expr_ (struct node_expr * value)
  { 
    type_ = type_scalar_expr_;
    value_expr_ = value; 
  };

  void set_logical_expr_ (struct node_expr * value)
  { 
    type_ = type_logical_expr_;
    value_expr_ = value; 
  };

  void evaluate_scalar_
  ( struct node_expr * node, int n, 
    double * result, double *x, double *y, double *z, double *t);

//   bool evaluate_logical_
//   ( struct node_expr * node, int n, 
//     bool * result, double *x, double *y, double *z, double *t);

  void dealloc_string_() { free (value_string_); } 
  void dealloc_param_ (struct param_type * p);
  void dealloc_node_expr_ (struct node_expr * p);

};

//----------------------------------------------------------------------

#endif
