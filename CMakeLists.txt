set(PROJECT_NAME "MuscleModelPython")
set(HEADER_DIRECTORY "muscle_model_python")

set(NAMESPACE_NAME "${PROJECT_NAME}")

set(LIBRARY_NAME "${PROJECT_NAME}")
set(EXECUTABLE_NAME "${PROJECT_NAME}Exec")
set(TEST_NAME "${PROJECT_NAME}Tests")

set(LIB_EXPORT_NAME "${LIBRARY_NAME}Targets")
set(LIB_CONFIG_NAME "${LIBRARY_NAME}Config")
set(LIB_VERSION_NAME "${LIB_CONFIG_NAME}Version")

cmake_minimum_required(VERSION 3.16)
project("${PROJECT_NAME}" VERSION 1.0.0)

include(GNUInstallDirs)

option(BUILD_EXECUTABLES "Build executables" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(Boost REQUIRED COMPONENTS system python)

set(LIB_SRC_FILES
	"src/Thelen2003MuscleGazebo.cpp"
	"src/muscle_model_python.cpp"
)

set(PYTHON_INSTALL_PATH "${CMAKE_INSTALL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/${LIBRARY_NAME}")

set(MUSCLE_MODEL_PYTHON_MODULE_DIR "${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PATH}")
set(MUSCLE_MODEL_DATA_DIR "data_files")
set(MUSCLE_MODEL_DATA_W_S_DIR "W_s")
set(MUSCLE_MODEL_DATA_W_E_DIR "W_e")
set(MUSCLE_MODEL_DATA_W_L_BAR_DIR "W_l_bar")
configure_file("include/muscle_model_plugin/muscle_model_python_config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/include/muscle_model_plugin/muscle_model_python_config.h" @ONLY)

##########################################
## Library
add_library(MuscleModelPython SHARED ${LIB_SRC_FILES})
add_library("${NAMESPACE_NAME}::MuscleModelPython" ALIAS MuscleModelPython)
target_compile_options(MuscleModelPython PUBLIC $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra>)

set_target_properties(MuscleModelPython PROPERTIES PREFIX "")

target_include_directories(MuscleModelPython
	PUBLIC
		"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"

	PRIVATE
)

target_link_libraries(MuscleModelPython
	PUBLIC
		Python3::Python

		Boost::boost
		Boost::python

	PRIVATE
)


##########################################
## Library
add_executable(TestMuscle "src/test_node_two_muscles.cpp")

target_link_libraries(TestMuscle ${NAMESPACE_NAME}::MuscleModelPython)


##########################################
## Install library files
set(INSTALL_CONFIGDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

message(!!!!!!!!!!!!!!!!!!!! ${PYTHON_INSTALL_PATH})

install(TARGETS ${LIBRARY_NAME}
	EXPORT ${LIB_EXPORT_NAME}
	LIBRARY DESTINATION "${PYTHON_INSTALL_PATH}"
	ARCHIVE DESTINATION "${PYTHON_INSTALL_PATH}"
	RUNTIME DESTINATION "${PYTHON_INSTALL_PATH}")

install(FILES "cmake/__init__.py" DESTINATION "${PYTHON_INSTALL_PATH}")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data_files" DESTINATION ${PYTHON_INSTALL_PATH})

# Install headers
#install(DIRECTORY ${HEADER_DIRECTORY} "${CMAKE_BINARY_DIR}/include/${HEADER_DIRECTORY}"
#	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
#	FILES_MATCHING
#		PATTERN "*.h"
#		PATTERN "*.hpp")

# Export target
#install(EXPORT ${LIB_EXPORT_NAME}
#	DESTINATION ${INSTALL_CONFIGDIR}
#	FILE "${LIB_EXPORT_NAME}.cmake"
#	NAMESPACE "${NAMESPACE_NAME}::")

# create version and config files
#include(CMakePackageConfigHelpers)
#write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${LIB_VERSION_NAME}.cmake"
#	VERSION ${PROJECT_VERSION}
#	COMPATIBILITY AnyNewerVersion)

#configure_package_config_file("${CMAKE_CURRENT_LIST_DIR}/cmake/ProjectConfig.cmake.in"
#	"${CMAKE_CURRENT_BINARY_DIR}/${LIB_CONFIG_NAME}.cmake"
#	INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

## Install all the helper files
#install(FILES	"${CMAKE_CURRENT_BINARY_DIR}/${LIB_CONFIG_NAME}.cmake"
#		"${CMAKE_CURRENT_BINARY_DIR}/${LIB_VERSION_NAME}.cmake"
#	DESTINATION ${INSTALL_CONFIGDIR})
