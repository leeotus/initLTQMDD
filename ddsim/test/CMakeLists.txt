set(DDSIM_TARGET ddsim)

# fix cmake warnings for googletest when building shared libs
if (APPLE AND BUILD_SHARED_LIBS)
	set(CMAKE_MACOSX_RPATH ON)
endif ()

if (NOT TARGET gtest OR NOT TARGET gmock)
	# Prevent overriding the parent project's compiler/linker settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	add_subdirectory("${CMAKE_SOURCE_DIR}/extern/dd_package/extern/googletest"
			"ddsim/extern/googletest" EXCLUDE_FROM_ALL)
	mark_as_advanced(
			BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
			gmock_build_tests gtest_build_samples gtest_build_tests
			gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
	)
	set_target_properties(gtest PROPERTIES FOLDER extern)
	set_target_properties(gtest_main PROPERTIES FOLDER extern)
	set_target_properties(gmock PROPERTIES FOLDER extern)
	set_target_properties(gmock_main PROPERTIES FOLDER extern)
	target_compile_features(gtest PUBLIC cxx_std_14)
	target_compile_features(gtest_main PUBLIC cxx_std_14)
endif()

macro(package_add_test TESTNAME)
	add_executable(${TESTNAME} ${ARGN})
	target_link_libraries(${TESTNAME} PRIVATE ${DDSIM_TARGET} gtest_main)
	set_property(TARGET ${TESTNAME} PROPERTY CXX_EXTENSIONS OFF)
	if (COVERAGE)
		target_compile_options(${TESTNAME} PRIVATE --coverage)
		target_link_libraries(${TESTNAME} PRIVATE --coverage)
	endif()
	if(MSVC)
		target_compile_options(${TESTNAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
		target_link_options(${TESTNAME} PUBLIC "$<$<CONFIG:Release>:/LTCG>")
	else()
		target_compile_options(${TESTNAME} PUBLIC -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-DNDEBUG -g0 -O3 -mtune=native -march=native>)
	endif()

	gtest_discover_tests(${TESTNAME}
	                     WORKING_DIRECTORY ${PROJECT_DIR}
	                     PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
	                     )
	set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

package_add_test(${DDSIM_TARGET}_test
	${CMAKE_CURRENT_SOURCE_DIR}/test_qfr_sim.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test_shor_sim.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test_fast_shor_sim.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test_grover_sim.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test_stoch_noise_sim.cpp
)

add_custom_command(TARGET ${DDSIM_TARGET}_test
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${DDSIM_TARGET}_test>/${DDSIM_TARGET}_test ${CMAKE_BINARY_DIR}/${DDSIM_TARGET}_test
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/circuits $<TARGET_FILE_DIR:${DDSIM_TARGET}_test>/circuits
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${DDSIM_TARGET}_test>/circuits ${CMAKE_BINARY_DIR}/circuits
                   COMMENT "Copying circuits and creating symlinks for ${DDSIM_TARGET}_test"
                   VERBATIM)