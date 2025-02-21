# See LICENSE_CELLO file for license and copyright information

# create Enzo::inference, which represents source files related to implementing
# inference.
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_inference
   EnzoLevelArray.cpp
   EnzoLevelArray.hpp
   EnzoMethodInference.cpp
   EnzoMethodInference.hpp
   EnzoMethodFBNetDeposit.cpp
   EnzoMethodFBNetDeposit.hpp
   EnzoObjectFeedbackSphere.hpp 
   EnzoStarNetUtils.cpp
   EnzoStarNetUtils.hpp
   Index3.cpp
   Index3.hpp
)
add_library(Enzo::inference ALIAS Enzo_inference)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_inference PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_inference PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_inference PRIVATE ${Cello_TARGET_LINK_OPTIONS})
