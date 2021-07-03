#Define external_library target, including all CPP source files and hpp,dir headers
add_library(external_library
    STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/defs.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/except.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/except.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/expr.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/expreval.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/expr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/func.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/funclist.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/funclist.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/node.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/node.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/parser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/parser.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/vallist.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/External/ExprEval/vallist.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/pngconf.h
    ${CMAKE_CURRENT_SOURCE_DIR}/External/pngwriter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/External/pngwriter.h
    )


#Define the charm interface CI files
set( charm_interfaces
    Cello/main.CI
    Cello/main_enzo.CI
    Cello/main_mesh.CI
    Cello/main_simulation.CI
    Cello/mesh.CI
    Cello/simulation.ci
    Enzo/enzo.CI
    )

# Use the preproccess_charm_interafaces function to create "*.decl.h" and
# "*.def.h" header files. These header files should have the "GENERATED"
# property for cmake, so they should be useable as dependencies for compiling
# source files
preprocess_charm_interfaces(charm_interface_headers ${charm_interfaces})

add_library(cello_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/cello.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/cello_Sync.cpp
	)

add_library(charm_library STATIC
    ${charm_interface_headers} #Include the charm interface headers as source dependencies
	${CMAKE_SOURCE_DIR}/src/Cello/charm_MappingArray.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/charm_MappingTree.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/charm_MsgCoarsen.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/charm_MsgRefine.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/charm_MsgRefresh.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/charm_mesh.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/charm_reductions.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/charm_simulation.cpp
	)

add_library(control_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/control_adapt.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/control_charm.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/control_compute.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/control_initialize.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/control_new_refresh.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/control_output.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/control_refresh.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/control_stopping.cpp
	)

add_library(disk_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/disk_FileHdf5.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/main.cpp
	)

add_library(error_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/error_Error.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/main.cpp
	)

add_library(data_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/data_Data.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_DataMsg.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_Face.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_FaceFluxes.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_FieldData.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_FieldDescr.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_FieldFace.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_FluxData.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_ParticleData.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_ParticleDescr.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/data_png_array.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/main_simulation.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/field_face_store.F
	)

add_library(data0_library INTERFACE
	)

add_library(io_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/io_ColormapRGB.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_Input.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_InputData.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_Io.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_IoBlock.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_IoFieldData.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_IoHierarchy.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_IoParticleData.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_Output.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_OutputCheckpoint.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_OutputData.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_OutputImage.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_Schedule.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_ScheduleInterval.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/io_ScheduleList.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/main_simulation.cpp
	)

add_library(main_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/main.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/main_enzo.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/main_mesh.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/main_simulation.cpp
	)

add_library(memory_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/memory_Memory.cpp
    #${CMAKE_SOURCE_DIR}/src/Cello/main.cpp
	)

add_library(mesh_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_Block.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_Factory.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_Hierarchy.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_Index.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_ItFace.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_ItNeighbor.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_ItNode.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_Node.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_NodeTrace.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_Refine.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_RefineDensity.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_RefineMask.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_RefineParticleCount.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_RefineShear.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_RefineSlope.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/mesh_Tree.cpp
    ${CMAKE_SOURCE_DIR}/src/Cello/main_mesh.cpp
	)

add_library(monitor_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/monitor_Monitor.cpp
    #${CMAKE_SOURCE_DIR}/src/Cello/main.cpp
	)

add_library(parallel_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/parallel_GroupProcess.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/parallel_GroupProcessCharm.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/parallel_GroupProcessMpi.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/parallel_GroupProcessSerial.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/parallel_ReduceMpi.cpp
	)

add_library(parameters_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/parameters.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/parameters_Config.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/parameters_Param.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/parameters_Parameters.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/parse.tab.c
	${CMAKE_SOURCE_DIR}/src/Cello/lex.yy.c
    #${CMAKE_SOURCE_DIR}/src/Cello/main.cpp
	)

add_library(performance_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/performance_Papi.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/performance_Performance.cpp
    #${CMAKE_SOURCE_DIR}/src/Cello/main.cpp
	)

add_library(problem_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/problem_BoundaryValue.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Initial.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_InitialFile.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_InitialTrace.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_InitialValue.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Mask.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_MaskExpr.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_MaskPng.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Method.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_MethodCloseFiles.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_MethodDebug.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_MethodFluxCorrect.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_MethodNull.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_MethodTrace.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Problem.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Prolong.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_ProlongInject.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_ProlongLinear.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Refresh.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Restrict.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_RestrictLinear.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_ScalarExpr.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Stopping.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/problem_Value.cpp
	)

add_library(compute_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/compute_Matrix.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/compute_Solver.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/compute_SolverNull.cpp
	)

# array_library is basically header-only
add_library(array_library INTERFACE
	)

add_library(simulation_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/simulation.cpp
	${CMAKE_SOURCE_DIR}/src/Cello/simulation_Simulation.cpp
	)

add_library(test_library STATIC
	${CMAKE_SOURCE_DIR}/src/Cello/test_Unit.cpp
	)

target_include_directories(array_library INTERFACE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(cello_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(charm_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(control_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(disk_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(error_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(data_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(io_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(main_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(memory_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(mesh_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(monitor_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(problem_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(compute_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(parallel_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(parameters_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(performance_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(simulation_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)
target_include_directories(test_library PRIVATE ${CMAKE_SOURCE_DIR}/src/Cello)




# specify the links between libraries (based on lines 60-80)
target_link_libraries(monitor_library
    PUBLIC ${Boost_LIBRARIES}
    INTERFACE parallel_library
)
target_link_libraries(error_library
    PUBLIC ${Boost_LIBRARIES}
    INTERFACE monitor_library
)
target_link_libraries(disk_library
    PUBLIC ${HDF5_LIBRARIES}
    INTERFACE error_library
    INTERFACE cello_library
)
# maybe add something related to PAPI in the future
target_link_libraries(performance_library
    PUBLIC ${Boost_LIBRARIES} #PAPI
    INTERFACE memory_library
    INTERFACE parameters_library
    INTERFACE error_library
)
target_link_libraries(array_library INTERFACE error_library)
target_link_libraries(parallel_library INTERFACE error_library)
target_link_libraries(data_library
    INTERFACE data_library
    INTERFACE io_library
    INTERFACE parallel_library
    INTERFACE cello_library
)
target_link_libraries(memory_library
    PUBLIC ${Boost_LIBRARIES}
    INTERFACE error_library
)
target_link_libraries(mesh_library
    PUBLIC ${PNG_LIBRARIES}
    PUBLIC ${HDF5_LIBRARIES}
    PUBLIC ${Boost_LIBRARIES}
    INTERFACE data_library
    INTERFACE io_library
    INTERFACE disk_library
    INTERFACE parallel_library
    INTERFACE external_library
    INTERFACE error_library
)
target_link_libraries(problem_library INTERFACE external_library)
target_link_libraries(parameters_library
    PUBLIC ${Boost_LIBRARIES}
    INTERFACE error_library
    INTERFACE cello_library
)
target_link_libraries(simulation_library INTERFACE problem_library)
target_link_libraries(test_library
    INTERFACE parallel_library
    INTERFACE cello_library
    PUBLIC charm_library
    INTERFACE data_library
    INTERFACE simulation_library
)

# this includes everything - may need to add something related to PAPI in the future
add_library(all_library INTERFACE)
target_link_libraries(all_library
    INTERFACE charm_library
    INTERFACE control_library
    INTERFACE simulation_library
    INTERFACE compute_library
    INTERFACE mesh_library
    INTERFACE data_library
    INTERFACE problem_library
    INTERFACE io_library
    INTERFACE disk_library
    INTERFACE memory_library
    INTERFACE parallel_library
    INTERFACE parameters_library
    INTERFACE error_library
    INTERFACE monitor_library
    INTERFACE performance_library
    INTERFACE test_library
    INTERFACE array_library
    INTERFACE cello_library
    INTERFACE external_library
    INTERFACE ${PNG_LIBRARIES}
    INTERFACE ${HDF5_LIBRARIES}
    INTERFACE ${Boost_LIBRARIES}
)
target_link_libraries(data_library INTERFACE mesh_library)
target_link_libraries(io_library INTERFACE mesh_library)


# Compile some unit tests
# this test doesn't work
#add_executable(test_FileHdf5 ${CMAKE_CURRENT_SOURCE_DIR}/Cello/test_FileHdf5.cpp)
#target_link_libraries(test_FileHdf5 disk_library test_library)

add_executable(test_Error ${CMAKE_CURRENT_SOURCE_DIR}/Cello/test_Error.cpp)
target_link_libraries(test_Error error_library test_library)

add_executable(test_memory ${CMAKE_CURRENT_SOURCE_DIR}/Cello/test_Memory.cpp)
target_link_libraries(test_memory memory_library test_library)

add_executable(test_monitor ${CMAKE_CURRENT_SOURCE_DIR}/Cello/test_Monitor.cpp)
target_link_libraries(test_monitor PUBLIC monitor_library test_library)
#target_link_libraries(test_monitor monitor_library all_library)

# these tests don't work:
add_executable(test_parameters ${CMAKE_CURRENT_SOURCE_DIR}/Cello/test_Parameters.cpp)
target_link_libraries(test_parameters parameters_library test_library)

add_executable(test_parse ${CMAKE_CURRENT_SOURCE_DIR}/Cello/test_Parse.cpp)
target_link_libraries(test_parse parameters_library test_library)


# this is a bit of a shortcut to help with the include directories:
include_directories(${CMAKE_SOURCE_DIR}/src/Enzo)
include_directories(${CMAKE_SOURCE_DIR}/src/Cello)
include_directories(${CMAKE_SOURCE_DIR}/src/External)
include_directories(${CMAKE_CURRENT_BINARY_DIR})



#Define enzo_e_exe executable, including all CPP,C,Fortran source files and hpp,dir headers
add_executable(enzo_e_exe
    ${charm_interface_headers} #Include the charm interface headers as source dependencies
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/calcdiss.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/calc_dt.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/calc_eigen.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/charm_enzo.hpp
   #${CMAKE_CURRENT_SOURCE_DIR}/Enzo/cic_deposit.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/cic_interp.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/cicinterp.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_defines.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo.hpp
    #${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo.dir
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo-e.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoBfieldMethod.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoBfieldMethodCT.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoBfieldMethodCT.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoBfieldMethod.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoBlock.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoBlock.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoBoundary.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoCenteredFieldRegistry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoCenteredFieldRegistry.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputeAcceleration.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputeAcceleration.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputeCicInterp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputeCicInterp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputeCoolingTime.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputeCoolingTime.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputePressure.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputePressure.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputeTemperature.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoComputeTemperature.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoConfig.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoConfig.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoControl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoEFltArrayMap.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoEFltArrayMap.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoEOSIdeal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoEOSIdeal.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoEquationOfState.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoFactory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoFactory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoFieldArrayFactory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialBCenter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialBCenter.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialBurkertBodenheimer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialBurkertBodenheimer.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialCloud.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialCloud.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialCollapse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialCollapse.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialCosmology.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialCosmology.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialGrackleTest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialGrackleTest.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialImplosion2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialImplosion2.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialInclinedWave.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialInclinedWave.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialIsolatedGalaxy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialIsolatedGalaxy.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialMusic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialMusic.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialPm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialPm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialPpmlTest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialPpmlTest.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialSedovArray2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialSedovArray2.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialSedovArray3.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialSedovArray3.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialSedovRandom.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialSedovRandom.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialShockTube.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialShockTube.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialSoup.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialSoup.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialTurbulence.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoInitialTurbulence.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoIntegrableUpdate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoIntegrableUpdate.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoLazyPassiveScalarFieldList.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoLazyPassiveScalarFieldList.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMatrixDiagonal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMatrixDiagonal.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMatrixIdentity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMatrixIdentity.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMatrixLaplace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMatrixLaplace.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodBackgroundAcceleration.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodBackgroundAcceleration.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodCheckGravity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodCheckGravity.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodComovingExpansion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodComovingExpansion.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodCosmology.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodCosmology.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodDistributedFeedback.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodDistributedFeedback.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodFeedback.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodFeedback.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodFire2Feedback.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodFire2Feedback.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodGrackle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodGrackle.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodGravity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodGravity.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodHeat.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodHeat.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodHydro.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodHydro.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodMHDVlct.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodMHDVlct.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodPmDeposit.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodPmDeposit.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodPmUpdate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodPmUpdate.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodPpm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodPpm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodPpml.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodPpml.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodStarMaker.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodStarMakerFIRE2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodStarMakerStochasticSF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodStarMakerStochasticSF.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodTurbulence.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoMethodTurbulence.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoPermutedCoordinates.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoPhysicsCosmology.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoPhysicsCosmology.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoProblem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoProblem.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoProlong.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoProlong.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoProlongMC1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoProlongMC1.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoProlongPoisson.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoProlongPoisson.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoReconstructor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoReconstructor.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoReconstructorNN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoReconstructorNN.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoReconstructorPLM.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRefineMass.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRefineMass.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRefineParticleMass.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRefineParticleMass.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRefineShock.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRefineShock.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRestrict.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRestrict.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRiemann.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRiemannHLLC.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRiemannHLLD.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRiemannHLL.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRiemann.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRiemannImpl.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRiemannLUT.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoRiemannUtils.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSimulation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSimulation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverBiCgStab.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverBiCgStab.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverCg.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverCg.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverDd.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverDd.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverDiagonal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverDiagonal.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverJacobi.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverJacobi.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverMg0.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSolverMg0.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSourceInternalEnergy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoSourceInternalEnergy.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoStopping.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoStopping.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_EnzoUnits.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_finalize.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_fortran.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/_enzo.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_IoEnzoBlock.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_IoEnzoBlock.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_reductions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_SetMinimumSupport.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_SolveHydroEquations.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_SolveMHDEquations.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/enzo_typedefs.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/euler.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/expand_terms.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/flux_hll.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/flux_twoshock.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/fortran.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/fortran_types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/grid_cic.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/interp1d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/interp2d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/interp3d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/interpolate.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/inteuler.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/int_lin3d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/intpos.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/intprim.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/intvar.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/ngpinterp.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/pgas2d_dual.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/pgas2d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/ppm_de.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_Conservative.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_HLLD.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_MAIN.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_monot.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_potokx.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_potoky.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_potokz.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_Primitive.F
   #${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_QDD6.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_sub.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/PPML_TimeStep.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/print_slice.F
   #${CMAKE_CURRENT_SOURCE_DIR}/Enzo/test_EnzoUnits.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/tscint1d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/tscint2d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/tscint3d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/turboinit2d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/turboinit.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/twoshock.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/xeuler_sweep.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/yeuler_sweep.F
    ${CMAKE_CURRENT_SOURCE_DIR}/Enzo/zeuler_sweep.F
    )


#Add the other directories in project for including
target_include_directories(enzo_e_exe
	PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)
target_include_directories(enzo_e_exe
	PRIVATE
    ${CMAKE_SOURCE_DIR}/src/External
)

target_include_directories(enzo_e_exe
	PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Cello
)

#Make sure mesh_library and external_library build before enzo_e_exe
add_dependencies(enzo_e_exe mesh_library)
add_dependencies(enzo_e_exe external_library)

#Link mesh_library and external_library to enzo_e_exe
target_link_libraries(enzo_e_exe mesh_library)
target_link_libraries(enzo_e_exe external_library)

#Link external libraries
#target_link_libraries(enzo_e_exe ${CHARM_LIBRARIES}) No need to link CHARM_LIBRARIES, linking with charmc takes care of it
target_link_libraries(enzo_e_exe ${ENZO_E_LIBRARIES})

#Link Boost
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    target_link_libraries(enzo_e_exe ${Boost_LIBRARIES})
endif()

# Link HDF5
if(HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  #target_link_libraries(mesh_library ${HDF5_LIBRARIES})
  target_link_libraries(enzo_e_exe ${HDF5_LIBRARIES})
endif()
