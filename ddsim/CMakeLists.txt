cmake_minimum_required(VERSION 3.10...3.16)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# When included as a subdirectory, the parent project already sets up the C++ standard.
# Keep separate project declaration only for standalone builds.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	project(ddsim
	        VERSION 1.0.2
	        DESCRIPTION "A front end to utilize the dd_package and qfr for simulation"
	        LANGUAGES CXX)

	set_property(GLOBAL PROPERTY USE_FOLDERS ON)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
	SET(COVERAGE OFF CACHE BOOL "Coverage")
	option(BUILD_SHARED_LIBS "Build shared libraries (.dll/.so) instead of static ones (.lib/.a)" OFF)

	set(default_build_type "Release")
	if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
		message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
		set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
		    STRING "Choose the type of build." FORCE)
		set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
	endif()
endif()

add_subdirectory(src)

# apps and test
option(BUILD_DDSIM_APPS "Build ddsim apps and tests" OFF)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR BUILD_DDSIM_APPS)
	option(BUILD_DDSIM_TESTS "Also build tests for DDSIM project" OFF)
	if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR BUILD_DDSIM_TESTS)
		enable_testing()
		include(GoogleTest)
		add_subdirectory(test)
	endif()

	set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Skip gbenchmark self tests")
	add_subdirectory(extern/benchmark ddsim_benchmark_build)
	add_subdirectory(apps)
endif()