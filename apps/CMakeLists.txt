set(APP_MAIN_NAME "ltqmdd")
set(LB_SIFTING_NAME "lb_sifting")
set(BENCHMARK_NAME "benchmark")
set(EQUIV_CHECK_NAME "equiv_check")
set(REAL2QASM_NAME "real2qasm")
add_executable(${PROJECT_NAME}_app ${CMAKE_CURRENT_SOURCE_DIR}/app.cpp)
add_executable(${APP_MAIN_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable(${LB_SIFTING_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/lb_main.cpp)
add_executable(${BENCHMARK_NAME} ${CMAKE_SOURCE_DIR}/test/benchmark.cpp)
add_executable(${EQUIV_CHECK_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/equiv_check.cpp)
add_executable(${REAL2QASM_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/real2qasm.cpp)

target_compile_features(${PROJECT_NAME}_app PUBLIC cxx_std_14)
target_compile_features(${APP_MAIN_NAME} PUBLIC cxx_std_14)
target_compile_features(${LB_SIFTING_NAME} PUBLIC cxx_std_14)
target_compile_features(${BENCHMARK_NAME} PUBLIC cxx_std_14)
target_compile_features(${EQUIV_CHECK_NAME} PUBLIC cxx_std_17)
target_compile_features(${REAL2QASM_NAME} PUBLIC cxx_std_17)

set_target_properties(${PROJECT_NAME}_app PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${APP_MAIN_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${LB_SIFTING_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${BENCHMARK_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${EQUIV_CHECK_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${REAL2QASM_NAME} PROPERTIES CXX_EXTENSIONS OFF)


target_link_libraries(${APP_MAIN_NAME} PUBLIC ${PROJECT_NAME})
if (COVERAGE)
        target_compile_options(${APP_MAIN_NAME} PUBLIC --coverage)
        target_link_libraries(${APP_MAIN_NAME} PUBLIC --coverage)
endif()
if(MSVC)
        target_compile_options(${APP_MAIN_NAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
        target_link_options(${APP_MAIN_NAME} PUBLIC "$<$<CONFIG:Release>:/LTCG>")
else()
        target_compile_options(${APP_MAIN_NAME} PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -g>)
endif()

add_custom_command(TARGET ${APP_MAIN_NAME}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${APP_MAIN_NAME}>/${APP_MAIN_NAME} ${CMAKE_BINARY_DIR}/${APP_MAIN_NAME}
                #    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/circuits/ $<TARGET_FILE_DIR:${APP_MAIN_NAME}>/circuits
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${APP_MAIN_NAME}>/circuits ${CMAKE_BINARY_DIR}/circuits
                   COMMENT "Copying circuits and creating symlinks for ${APP_MAIN_NAME}"
                   VERBATIM)

target_link_libraries(${LB_SIFTING_NAME} PUBLIC ${PROJECT_NAME})
if(MSVC)
        target_compile_options(${LB_SIFTING_NAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
        target_link_options(${LB_SIFTING_NAME} PUBLIC "$<$<CONFIG:Release>:/LTCG>")
else()
        target_compile_options(${LB_SIFTING_NAME} PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -g>)
endif()

add_custom_command(TARGET ${LB_SIFTING_NAME}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${LB_SIFTING_NAME}>/${LB_SIFTING_NAME} ${CMAKE_BINARY_DIR}/${LB_SIFTING_NAME}
                   COMMENT "Creating symlink for ${LB_SIFTING_NAME}"
                   VERBATIM)

target_link_libraries(${BENCHMARK_NAME} PUBLIC ${PROJECT_NAME})
if(MSVC)
        target_compile_options(${BENCHMARK_NAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
        target_link_options(${BENCHMARK_NAME} PUBLIC "$<$<CONFIG:Release>:/LTCG>")
else()
        target_compile_options(${BENCHMARK_NAME} PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -DNDEBUG>)
endif()

add_custom_command(TARGET ${BENCHMARK_NAME}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${BENCHMARK_NAME}>/${BENCHMARK_NAME} ${CMAKE_BINARY_DIR}/${BENCHMARK_NAME}
                   COMMENT "Creating symlink for ${BENCHMARK_NAME}"
                   VERBATIM)

target_link_libraries(${PROJECT_NAME}_app PRIVATE ${PROJECT_NAME})
if(MSVC)
	target_compile_options(${PROJECT_NAME}_app PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
	target_link_options(${PROJECT_NAME}_app PUBLIC "$<$<CONFIG:Release>:/LTCG>")
else()
	target_compile_options(${PROJECT_NAME}_app PUBLIC -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -DNDEBUG -g0>)
endif()

add_custom_command(TARGET ${PROJECT_NAME}_app
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_app ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_app
                   COMMENT "Creating symlink for ${PROJECT_NAME}_app"
                   VERBATIM)

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}_app
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Tight LB benchmark
set(TIGHT_LB_BENCH_NAME "benchmark_tight_lb")
add_executable(${TIGHT_LB_BENCH_NAME} ${CMAKE_SOURCE_DIR}/test/benchmark_tight_lb.cpp)
target_compile_features(${TIGHT_LB_BENCH_NAME} PUBLIC cxx_std_14)
set_target_properties(${TIGHT_LB_BENCH_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(${TIGHT_LB_BENCH_NAME} PUBLIC ${PROJECT_NAME})
if(MSVC)
        target_compile_options(${TIGHT_LB_BENCH_NAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
else()
        target_compile_options(${TIGHT_LB_BENCH_NAME} PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -DNDEBUG>)
endif()
add_custom_command(TARGET ${TIGHT_LB_BENCH_NAME}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${TIGHT_LB_BENCH_NAME}>/${TIGHT_LB_BENCH_NAME} ${CMAKE_BINARY_DIR}/${TIGHT_LB_BENCH_NAME}
                   COMMENT "Creating symlink for ${TIGHT_LB_BENCH_NAME}"
                   VERBATIM)
set_target_properties(${PROJECT_NAME}_app PROPERTIES EXPORT_NAME ${PROJECT_NAME}_app)

# IG Sifting benchmark
set(IG_BENCH_NAME "benchmark_ig_sifting")
add_executable(${IG_BENCH_NAME} ${CMAKE_SOURCE_DIR}/test/benchmark_ig_sifting.cpp)
target_compile_features(${IG_BENCH_NAME} PUBLIC cxx_std_14)
set_target_properties(${IG_BENCH_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(${IG_BENCH_NAME} PUBLIC ${PROJECT_NAME})
if(MSVC)
        target_compile_options(${IG_BENCH_NAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
else()
        target_compile_options(${IG_BENCH_NAME} PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -DNDEBUG>)
endif()
add_custom_command(TARGET ${IG_BENCH_NAME}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${IG_BENCH_NAME}>/${IG_BENCH_NAME} ${CMAKE_BINARY_DIR}/${IG_BENCH_NAME}
                   COMMENT "Creating symlink for ${IG_BENCH_NAME}"
                   VERBATIM)

# Equivalence Checker
target_link_libraries(${EQUIV_CHECK_NAME} PUBLIC ${PROJECT_NAME})
if(MSVC)
        target_compile_options(${EQUIV_CHECK_NAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
        target_link_options(${EQUIV_CHECK_NAME} PUBLIC "$<$<CONFIG:Release>:/LTCG>")
else()
        target_compile_options(${EQUIV_CHECK_NAME} PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native>)
endif()
add_custom_command(TARGET ${EQUIV_CHECK_NAME}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${EQUIV_CHECK_NAME}>/${EQUIV_CHECK_NAME} ${CMAKE_BINARY_DIR}/${EQUIV_CHECK_NAME}
                   COMMENT "Creating symlink for ${EQUIV_CHECK_NAME}"
                   VERBATIM)

# Real to QASM converter
target_link_libraries(${REAL2QASM_NAME} PUBLIC ${PROJECT_NAME})
if(NOT MSVC)
        target_compile_options(${REAL2QASM_NAME} PRIVATE -Wall $<$<CONFIG:RELEASE>:-O3>)
endif()
add_custom_command(TARGET ${REAL2QASM_NAME}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${REAL2QASM_NAME}>/${REAL2QASM_NAME} ${CMAKE_BINARY_DIR}/${REAL2QASM_NAME}
                   COMMENT "Creating symlink for ${REAL2QASM_NAME}"
                   VERBATIM)
