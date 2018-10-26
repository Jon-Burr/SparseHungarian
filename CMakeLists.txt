cmake_minimum_required( VERSION 3.9 )
project( SparseHungarian VERSION 1.0.0 LANGUAGES CXX )

find_package( Eigen3 )
find_package( Boost REQUIRED program_options )

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

add_executable( MatchTestPoints util/MatchTestPoints.cxx )
target_link_libraries( MatchTestPoints SparseHungarianLib Boost::program_options)
target_compile_features( MatchTestPoints 
    PRIVATE cxx_auto_type )
