#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoReconstructor* EnzoReconstructor::construct_reconstructor(std::string name)
{
  // some repeated code from construct_riemann

  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(name.size(), ' ');
  std::transform(name.begin(), name.end(), formatted.begin(),
		 ::tolower);
  EnzoReconstructor* out;
  if (formatted == std::string("nn")){
    out = new EnzoReconstructorNN();
  } else if (formatted == std::string("plm")){
    out = new EnzoReconstructorPLM();
  } else {
    ASSERT("EnzoReconstructor",
	   "The only allowed solvers are NN & PLM", false);
    out = NULL; // Deals with compiler warning
  }
  return out;
}
