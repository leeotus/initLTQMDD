# allow loader to find lib
if(APPLE)
	set(basePoint @loader_path)
else()
	set(basePoint $ORIGIN)
endif()
include(GNUInstallDirs)
file(RELATIVE_PATH relDir
     ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
     ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
     )
set(CMAKE_INSTALL_RPATH ${basePoint} ${basePoint}/${relDir})

# Explicit target name to avoid conflicts when included as subdirectory
set(DDSIM_TARGET ddsim)

# Source root of the ddsim project (parent of this src/ directory)
get_filename_component(DDSIM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

add_library(${DDSIM_TARGET}
		${DDSIM_SOURCE_DIR}/include/Simulator.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/Simulator.cpp
		${DDSIM_SOURCE_DIR}/include/QFRSimulator.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/QFRSimulator.cpp
		${DDSIM_SOURCE_DIR}/include/GroverSimulator.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/GroverSimulator.cpp
		${DDSIM_SOURCE_DIR}/include/ShorSimulator.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/ShorSimulator.cpp
        ${DDSIM_SOURCE_DIR}/include/ShorFastSimulator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ShorFastSimulator.cpp
		)
add_library(JKQ::${DDSIM_TARGET} ALIAS ${DDSIM_TARGET})

target_include_directories(${DDSIM_TARGET} PUBLIC
                           $<BUILD_INTERFACE:${DDSIM_SOURCE_DIR}/include>)
target_compile_features(${DDSIM_TARGET} PUBLIC cxx_std_14)
set_target_properties(${DDSIM_TARGET} PROPERTIES CXX_EXTENSIONS OFF)

# Link to the parent project's qfr library
target_link_libraries(${DDSIM_TARGET} PUBLIC qfr)

# add coverage compiler and linker flag if COVERAGE is set
if (COVERAGE)
	target_compile_options(${DDSIM_TARGET} PRIVATE --coverage)
	target_link_libraries(${DDSIM_TARGET} PRIVATE --coverage)
endif()

# set compiler flags depending on compiler
if(MSVC)
	target_compile_options(${DDSIM_TARGET} PRIVATE "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
	target_link_options(${DDSIM_TARGET} PRIVATE "$<$<CONFIG:Release>:/LTCG>")
else()
	target_compile_options(${DDSIM_TARGET} PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-DNDEBUG -g0 -O3 -mtune=native -march=native>)
endif()