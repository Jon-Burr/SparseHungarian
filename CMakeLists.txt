cmake_minimum_required( VERSION 3.9 )
project( SparseHungarian VERSION 1.0.0 LANGUAGES CXX )

set( Eigen3_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
find_package( Eigen3 )

add_library( SparseHungarianLib SHARED
    src/SparseGroup.cxx src/Matching.cxx src/HungarianSolver.cxx
    )
target_include_directories( SparseHungarianLib
    PUBLIC
      $<INSTALL_INTERFACE:include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )
target_link_libraries( SparseHungarianLib
    PUBLIC
      Eigen3::Eigen
    )
target_compile_features( SparseHungarianLib
    PUBLIC cxx_alias_templates
    PRIVATE cxx_auto_type
    )
#install( TARGETS SparseHungarianLib
#         EXPORT SparseHungarianLibTargets
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib
#    )
#install( FILES 


# Build the simple reading script
# add_subdirectory( util )


#find_package( Eigen3 )


#add_library( SparseHungarianLib SHARED
#    src/SparseGroup.cxx src/Matching.cxx src/HungarianSolver.cxx
#    )

#target_include_directories( SparseHungarianLib
#    PUBLIC
#      $<INSTALL_INTERFACE:SparseHungarian>
#      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SparseHungarian>
#    )

#target_link_libraries( SparseHungarianLib
#    PUBLIC
#      Eigen3::Eigen
#    )
#target_compile_features( SparseHungarianLib
#    PUBLIC cxx_alias_templates
#   PRIVATE cxx_auto_type 
#    )

# Install the library
#install (TARGETS SparseHungarianLib
#    EXPORT SparseHungarian-export
#   LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/bin
#    )

#add_executable( testSparseHungarian util/test.cxx )
#target_link_libraries( testSparseHungarian SparseHungarianLib )
#target_compile_features( testSparseHungarian
#    PRIVATE cxx_auto_type
#    )
