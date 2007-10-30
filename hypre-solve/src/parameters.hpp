//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Parameters class include file

/**
 * The Parameters class is used for storing key-value pairs of strings, where
 * the value may not necessarily be unique.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-08-30 
 *
 */

#define BUFFER_LENGTH 255

class Parameters 
{

  friend class ItParameters;

private:

  std::multimap<std::string, std::string> values_;

public:

  Parameters () throw () ;
  ~Parameters () throw () ;

  /// Read a file of key-value pairs.  May be more than one value per key.
  void read (std::string file) throw () ;

  /// Print all parameters to stdout.
  void print () throw () ;

  /// Associate the given value with the given key.
  void add_parameter (std::string key, std::string val) throw () ;

  /// Retrieve the ith value of the the given parameter.
  std::string ith_value  (std::string key, int i) const throw () ;

  /// Return the multiplicity of values for the given key.  May be 0.
  int num_values  (std::string key) const throw () ;

  /// Return the value for the given key.
  std::string value (std::string key) const throw ();

private:

  int readline_ (FILE *, char * buffer, int n) throw ();

};

/// ItParameters class

/**
 * 
 * An ItParameters object allows iterating through all key-value pairs in
 * a Parameters object.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-08-30
 *
 */

class ItParameters
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  Parameters * pparameters_;
  std::multimap<std::string, std::string>::iterator curr_;
  std::multimap<std::string, std::string>::iterator next_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItParameters (Parameters & parameters) throw ()
    : pparameters_(& parameters),
      curr_(),
      next_(parameters.values_.begin())
  { }

  ~ItParameters () throw () {};
  
  /// Iterate through all Grids in the Grid.
  int operator++ (int) { 
    
    curr_ = next_;
    int not_done = 1;
    if (next_ == pparameters_->values_.end()) {
      next_ = pparameters_->values_.begin();
      not_done = 0;
    }
    ++next_;
    return not_done;
  }

  std::string key () throw ()
  { return curr_->first; }
  std::string value () throw ()
  { return curr_->second; }

};

