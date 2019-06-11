cmake_minimum_required( VERSION 2.8 )

if (CPLEX)
  project( DECIFER-CPLEX )
else()
  project( DECIFER-GUROBI )
endif()

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

option(CPLEX "Use CPLEX" ON)
option(PYTHON "Generate Python module" OFF)

if (CPLEX)
  ADD_DEFINITIONS(-DCPLEX -DIL_STD)
  if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR)
    set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ILOG;/opt/ibm/ILOG;$ENV{HOME}/Applications/IBM/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
    foreach( CplexInstallPrefix ${CPLEX_AUTODETECTION_PREFIXES} )
      file( GLOB_RECURSE CplexLibraryPaths ${CplexInstallPrefix}/*/cplex/lib/*/libcplex.a )
      if( CplexLibraryPaths )
        list( LENGTH CplexLibraryPaths NumPath )
        if( NumPath GREATER 1 )
          message( WARNING "Found multiple CPLEX installations in ${CplexInstallPrefix}, automatically selecting the first one, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually if this is not the right selection" )
        else()
          message( "Found CPLEX's installation in ${CPlexInstallPrefix}" )
        endif()
        list( GET CplexLibraryPaths 0 CplexLibraryPath )

        set( MyRegex "${CplexInstallPrefix}/(.*)/cplex/?/lib/(.*)/libcplex.a" )
        string( REGEX MATCH "${MyRegex}" ValidatedCplexLibPath ${CplexLibraryPath} )
        if( NOT ValidatedCplexLibPath )
          message( WARNING "Unable to recognize ${CplexLibraryPath} as a valid CPLEX path" )
        else()
          string( REGEX REPLACE "${MyRegex}" "\\1" CplexHome ${ValidatedCplexLibPath} )
          string( REGEX REPLACE "${MyRegex}" "\\2" CplexArch ${ValidatedCplexLibPath} )

          set( CPLEX_INC_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/include"            CACHE PATH "CPLEX's include dir" )
          set( CPLEX_LIB_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/lib/${CplexArch}"   CACHE PATH "CPLEX's library dir" )
          set( CONCERT_INC_DIR "${CplexInstallPrefix}/${CplexHome}/concert/include"          CACHE PATH "Concert's include dir" )
          set( CONCERT_LIB_DIR "${CplexInstallPrefix}/${CplexHome}/concert/lib/${CplexArch}" CACHE PATH "Concert's library dir" )

          break()
        endif()
      else()
        message( "Unable to find CPLEX and Concert installation in ${CplexInstallPrefix}" )
      endif()
    endforeach()

    if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
      message( "Unable to autodeted CPLEX's and Concert's installation directories, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually" )
    else()
      message( "Autodetected CPLEX's include dir   : ${CPLEX_INC_DIR}" )
      message( "Autodetected CPLEX's library dir   : ${CPLEX_LIB_DIR}" )
      message( "Autodetected Concert's include dir : ${CONCERT_INC_DIR}" )
      message( "Autodetected Concert's library dir : ${CONCERT_LIB_DIR}" )
    endif()
  endif()
else()
  find_package( GUROBI REQUIRED )
endif()

if (PYTHON)
  if (CPLEX)
    set( python_mod_src
      src/clusterilp.cpp
      src/clusterilpcplex.cpp
      src/statetree.cpp
      src/stategraph.cpp
      src/statetreessampler.cpp
      src/readmatrix.cpp
      src/basematrix.cpp
      src/utils.cpp
      src/solver.cpp
      src/python.cpp
      src/posteriorstatetree.cpp
      src/phylogenetictree.cpp
    )

    set( python_mod_hdr
      src/clusterilp.h
      src/clusterilpcplex.h
      src/statetree.h
      src/stategraph.h
      src/statetreessampler.h
      src/readmatrix.h
      src/basematrix.h
      src/solver.h
      src/utils.h
      src/posteriorstatetree.h
      src/phylogenetictree.h
    )
  else()
    set( python_mod_src
      src/clusterilp.cpp
      src/clusterilpgurobi.cpp
      src/statetree.cpp
      src/stategraph.cpp
      src/statetreessampler.cpp
      src/readmatrix.cpp
      src/basematrix.cpp
      src/utils.cpp
      src/solver.cpp
      src/python.cpp
      src/posteriorstatetree.cpp
      src/phylogenetictree.cpp
    )

    set( python_mod_hdr
      src/clusterilp.h
      src/clusterilpgurobi.h
      src/statetree.h
      src/stategraph.h
      src/statetreessampler.h
      src/readmatrix.h
      src/basematrix.h
      src/solver.h
      src/utils.h
      src/posteriorstatetree.h
      src/phylogenetictree.h
    )
  endif()
endif()

set( visualizepoststatetree_src
  src/visualizepoststatetreemain.cpp
  src/statetree.cpp
  src/posteriorstatetree.cpp
  src/utils.cpp
  src/readmatrix.cpp
  src/basematrix.cpp
)

set( visualizepoststatetree_hdr
  src/utils.h
  src/statetree.h
  src/posteriorstatetree.h
  src/readmatrix.h
  src/basematrix.h
)

set( visualizestatetree_src
  src/visualizestatetreemain.cpp
  src/statetree.cpp
  src/utils.cpp
)

set( visualizestatetree_hdr
  src/utils.h
  src/statetree.h
)

set( mergestatetrees_src
  src/mergestatetreesmain.cpp
  src/statetree.cpp
  src/stategraph.cpp
  src/utils.cpp
  src/readmatrix.cpp
  src/basematrix.cpp
)

set( mergestatetrees_hdr
  src/statetree.h
  src/stategraph.h
  src/utils.h
  src/readmatrix.h
  src/basematrix.h
)

set( generatestatetrees_src
  src/generatestatetreesmain.cpp
  src/statetree.cpp
  src/stategraph.cpp
  src/utils.cpp
  src/readmatrix.cpp
  src/basematrix.cpp
  src/solver.cpp
)

set( generatestatetrees_hdr
  src/statetree.h
  src/stategraph.h
  src/utils.h
  src/readmatrix.h
  src/basematrix.h
  src/solver.h
)

if (CPLEX)
  set( decifer_src
    src/incrementalsolver.cpp
    src/decifermain.cpp
    src/minclusterilp.cpp
    src/minclusterilpcplex.cpp
    src/hardclusterilp.cpp
    src/hardclusterilpcplex.cpp
    src/clusterilp.cpp
    src/clusterilpcplex.cpp
    src/softclusterilpcplexcallback.cpp
    src/softclusterilp.cpp
    src/softclusterilpcplex.cpp
    src/softclusterlpcplex.cpp
    src/softclusterlpcplexext.cpp
    src/statetree.cpp
    src/stategraph.cpp
    src/readmatrix.cpp
    src/basematrix.cpp
    src/utils.cpp
    src/em.cpp
    src/emcplex.cpp
    src/solver.cpp
    src/posteriorstatetree.cpp
    src/precluster.cpp
    src/hardpreclusterilpcplex.cpp
    src/softpreclusterilpcplex.cpp
  )

  set( decifer_hdr
    src/incrementalsolver.h
    src/utils.h
    src/statetree.h
    src/minclusterilp.h
    src/minclusterilpcplex.h
    src/hardclusterilp.h
    src/hardclusterilpcplex.h
    src/clusterilp.h
    src/clusterilpcplex.h
    src/softclusterilp.h
    src/softclusterilpcplex.h
    src/softclusterlpcplex.h
    src/softclusterlpcplexext.h
    src/stategraph.h
    src/readmatrix.h
    src/basematrix.h
    src/em.h
    src/emcplex.h
    src/dkm/dkm.hpp
    src/dkm/dkm_utils.hpp
    src/solver.h
    src/posteriorstatetree.h
    src/precluster.h
    src/hardpreclusterilpcplex.h
    src/softpreclusterilpcplex.h
  )

  set( precluster_src
    src/preclustermain.cpp
    src/hardclusterilp.cpp
    src/statetree.cpp
    src/stategraph.cpp
    src/readmatrix.cpp
    src/basematrix.cpp
    src/utils.cpp
    src/solver.cpp
    src/precluster.cpp
    src/hardpreclusterilpcplex.cpp
    src/softpreclusterilpcplex.cpp
    src/softclusterilp.cpp
  )

  set( precluster_hdr
    src/utils.h
    src/statetree.h
    src/hardclusterilp.h
    src/stategraph.h
    src/readmatrix.h
    src/basematrix.h
    src/solver.h
    src/precluster.h
    src/hardpreclusterilpcplex.h
    src/softpreclusterilpcplex.h
    src/softclusterilp.h
  )
else()
  set( decifer_src
    src/decifermain.cpp
    src/minclusterilp.cpp
    src/minclusterilpgurobi.cpp
    src/hardclusterilp.cpp
    src/hardclusterilpgurobi.cpp
    src/clusterilp.cpp
    src/clusterilpgurobi.cpp
    src/softclusterlp.cpp
    src/statetree.cpp
    src/stategraph.cpp
    src/readmatrix.cpp
    src/basematrix.cpp
    src/utils.cpp
    src/em.cpp
    src/emgurobi.cpp
    src/solver.cpp
    src/posteriorstatetree.cpp
  )

  set( decifer_hdr
    src/utils.h
    src/statetree.h
    src/minclusterilp.h
    src/minclusterilpgurobi.h
    src/hardclusterilp.h
    src/hardclusterilpgurobi.h
    src/clusterilp.h
    src/clusterilpgurobi.h
    src/softclusterlp.h
    src/stategraph.h
    src/readmatrix.h
    src/basematrix.h
    src/em.h
    src/emgurobi.h
    src/dkm/dkm.hpp
    src/dkm/dkm_utils.hpp
    src/solver.h
    src/posteriorstatetree.h
  )
endif()

set( generate_src
  src/generatemain.cpp
  src/statetree.cpp
  src/stategraph.cpp
  src/statetreessampler.cpp
  src/readmatrix.cpp
  src/basematrix.cpp
  src/phylogenetictree.cpp
  src/utils.cpp
  src/generativemodel.cpp
)

set( generate_hdr
  src/beta_distribution.hpp
  src/utils.h
  src/statetree.h
  src/stategraph.h
  src/statetreessampler.h
  src/phylogenetictree.h
  src/readmatrix.h
  src/basematrix.h
  src/generativemodel.h
)

set( generate_hdr
  src/beta_distribution.hpp
  src/utils.h
  src/statetree.h
  src/stategraph.h
  src/readmatrix.h
  src/basematrix.h
)

set( ccf_src
  src/ccfmain.cpp
  src/statetree.cpp
  src/stategraph.cpp
  src/readmatrix.cpp
  src/basematrix.cpp
  src/utils.cpp
)

set( ccf_hdr
  src/statetree.h
  src/stategraph.h
  src/readmatrix.h
  src/basematrix.h
  src/utils.h
)

set( summarize_src
  src/summarizemain.cpp
  src/posteriorstatetree.cpp
  src/solver.cpp
  src/statetree.cpp
  src/stategraph.cpp
  src/readmatrix.cpp
  src/basematrix.cpp
  src/utils.cpp
)

set( summarize_hdr
  src/solver.h
  src/posteriorstatetree.h
  src/statetree.h
  src/stategraph.h
  src/readmatrix.h
  src/basematrix.h
  src/utils.h
)

# first we can indicate the documentation build as an option and set it to OFF by default
option(BUILD_DOC "Build documentation" OFF)

if (BUILD_DOC)
  if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
    # check if Doxygen is installed
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        # set input and output files
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen build started")

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
      message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
  endif()
endif()

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set( Boost_USE_STATIC_LIBS ON )

if (PYTHON)
  find_package(PythonLibs)
  if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    find_package( Boost REQUIRED COMPONENTS thread system filesystem python2.7 )
  else()
    find_package( Boost REQUIRED COMPONENTS thread system filesystem python )
  endif()
else()
  find_package( Boost REQUIRED COMPONENTS thread system filesystem )
endif()

MESSAGE( STATUS "found boost:"
  "\nINCLUDE: ${Boost_INCLUDE_DIRS}"
  "\nLIB: ${Boost_LIBRARIES}"
)

include_directories( "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} "${LIBLEMON_ROOT}/lib" ${GUROBI_LIBRARIES} ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} )

if (PYTHON)
  FIND_PACKAGE(PythonInterp REQUIRED)
  FIND_PACKAGE(PythonLibsNew REQUIRED)
  # custom macro with most of the redundant code for making a python example module
  MACRO( addPythonExe _name _srccpp )
    ADD_EXECUTABLE(${_name} ${_srccpp} ${python_mod_src} ${python_mod_hdr} )
    # make the pyd library link against boost_numpy python and boost
    if (CPLEX)
      TARGET_INCLUDE_DIRECTORIES(${_name} PUBLIC "${PROJECT_SOURCE_DIR}/src" ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${CPLEX_INC_DIR} ${CONCERT_INC_DIR})
      TARGET_LINK_LIBRARIES(${_name} ${CommonLibs} ${Boost_LIBRARIES} ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} ilocplex cplex concert dl)
    else()
      TARGET_INCLUDE_DIRECTORIES(${_name} PUBLIC "${PROJECT_SOURCE_DIR}/src" ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR})
      TARGET_LINK_LIBRARIES(${_name} ${CommonLibs} boost_python ${GUROBI_LIBRARIES})
    endif()
  ENDMACRO()

  MACRO( addPythonMod _name _srccpp )
    PYTHON_ADD_MODULE(${_name} ${_srccpp} ${python_mod_src} ${python_mod_hdr} )
    # make the pyd library link against boost_numpy python and boost
    if (CPLEX)
      TARGET_INCLUDE_DIRECTORIES(${_name} PUBLIC "${PROJECT_SOURCE_DIR}/src" ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${CPLEX_INC_DIR} ${CONCERT_INC_DIR})
      TARGET_LINK_LIBRARIES(${_name} ${CommonLibs} ${Boost_LIBRARIES} ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} ilocplex cplex concert dl)
    else()
      TARGET_INCLUDE_DIRECTORIES(${_name} PUBLIC "${PROJECT_SOURCE_DIR}/src" ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR})
      TARGET_LINK_LIBRARIES(${_name} ${CommonLibs} boost_python ${GUROBI_LIBRARIES})
    endif()
  ENDMACRO()
endif()

set( CommonLibs
  m
  emon
  pthread
  ${Boost_LIBRARIES}
)

set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )

add_executable( visualizestatetree EXCLUDE_FROM_ALL ${visualizestatetree_src} ${visualizestatetree_hdr} )
target_link_libraries( visualizestatetree ${CommonLibs} )

add_executable( visualizepoststatetree ${visualizepoststatetree_src} ${visualizepoststatetree_hdr} )
target_link_libraries( visualizepoststatetree ${CommonLibs} )

add_executable( generate EXCLUDE_FROM_ALL ${generate_src} ${generate_hdr} )
target_link_libraries( generate ${CommonLibs} )

add_executable( mergestatetrees ${mergestatetrees_src} ${mergestatetrees_hdr} )
target_link_libraries( mergestatetrees ${CommonLibs} )

add_executable( generatestatetrees ${generatestatetrees_src} ${generatestatetrees_hdr} )
target_link_libraries( generatestatetrees ${CommonLibs} )

add_executable( ccf EXCLUDE_FROM_ALL ${ccf_src} ${ccf_hdr} )
target_link_libraries( ccf ${CommonLibs} )

add_executable( decifer ${decifer_src} ${decifer_hdr} )
if (CPLEX)
  target_include_directories( decifer PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} )
  target_link_libraries( decifer ${CommonLibs} ilocplex cplex concert dl )

  add_executable( precluster ${precluster_src} ${precluster_hdr} )
  target_include_directories( precluster PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} )
  target_link_libraries( precluster ${CommonLibs} ilocplex cplex concert dl )
else()
  target_include_directories( decifer PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
  target_link_libraries( decifer ${CommonLibs} ${GUROBI_LIBRARIES} )
endif()

add_executable( summarize ${summarize_src} ${summarize_hdr} )
target_include_directories( summarize PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} )
target_link_libraries( summarize ${CommonLibs} dl )

if (PYTHON)
  addPythonMod( decifermod src/python.cpp )
endif()
