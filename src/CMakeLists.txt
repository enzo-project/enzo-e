


add_subdirectory(External)

#add_library(cello_library)
#add_dependencies(cello_library external_library)
add_subdirectory(Cello)


add_subdirectory(Enzo)
add_dependencies(enzo_e_exe cello_library)
add_dependencies(enzo_e_exe external_library)

