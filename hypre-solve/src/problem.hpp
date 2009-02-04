//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Problem class include file

/**
 * 
 * @file      problem.hpp
 * @brief     Declaration of the Problem class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

#define BUFFER_LENGTH 255

class Problem {

  //----------------------------------------------------------------------


private:

  Parameters            parameters_; // Raw input parameters
  std::vector<Sphere *> spheres_;    // List of sphere masses
  std::vector<Point *>  points_;     // List of point masses
  Hierarchy             hierarchy_;  // AMR mesh hierarchy
  Domain                domain_;     // Problem domain

  //----------------------------------------------------------------------

public:

  Problem () throw ();

  ~Problem () throw ();

  Problem (const Problem & ) throw ();

  Problem & operator = (const Problem & ) throw ();

  //--------------------------------------------------

  void print () throw ();

  void write (FILE *fp = 0) throw ();

  /// Read in a problem file

  void read (std::string filename) throw ();

  /// Return the dimension

  int dimension ()          
  { 
    return hierarchy_.dimension(); 
  }

  /// Return the hierarchy
  
  Hierarchy & hierarchy ()  
  { 
    return hierarchy_; 
  }

  /// Return the domain
  
  Domain & domain ()        
  { 
    return domain_; 
  }

  /// Return a pointer to the ith Sphere
  
  Sphere & sphere (int i)   
  { 
    return * spheres_[i]; 
  }

  /// Return vector of spheres
  
  std::vector<Sphere *> spheres ()     
  { 
    return spheres_; 
  }

  /// Return the number of Spheres
  
  int num_spheres ()        
  { 
    return spheres_.size(); 
  }

  /// Return a pointer to the ith Point
  
  Point & point (int i)     
  { 
    return * points_[i]; 
  }

  /// Return vector of points
  
  std::vector<Point *> points ()     
  { 
    return points_; 
  }

  /// Return the number of Points
  
  int num_points ()         
  { 
    return points_.size(); 
  }

  /// Return a pointer to the ith Grid
  
  Grid & grid (int i)       
  { 
    return hierarchy_.grid(i); 
  }

  /// Return the number of Grids
  
  int num_grids ()          
  { 
    return hierarchy_.num_grids(); 
  }

  /// Return a pointer to the ith Level
  
  Level & level (int i)     
  { 
    return hierarchy_.level(i); 
  }

  /// Return the number of Levels
  
  int num_levels ()         
  { 
    return hierarchy_.num_levels(); 
  }

  /// Return the Parameters object
  
  Parameters & parameters ()     
  { 
    return parameters_; 
  }

};
