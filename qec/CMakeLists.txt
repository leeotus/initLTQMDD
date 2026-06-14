# This CMakeLists.txt is designed to be included from the project root
# via add_subdirectory(qec). It assumes the root CMakeLists.txt has already
# configured extern/dd_package and extern/json.

# QEC library
set(QEC_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QECCode.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/SteaneCode.cpp
)

add_library(qec STATIC ${QEC_SOURCES})

target_include_directories(qec PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/operations
    ${CMAKE_SOURCE_DIR}/include/parsers
    ${CMAKE_SOURCE_DIR}/include/algorithms
    ${CMAKE_SOURCE_DIR}/extern/dd_package/include
    ${CMAKE_SOURCE_DIR}/extern/json/single_include
)

target_link_libraries(qec PUBLIC JKQ::DDpackage)

# QEC Benchmark executable
add_executable(qec_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/apps/qec_benchmark.cpp
)

# DDSIM sources (simulator library)
set(DDSIM_SOURCES
    ${CMAKE_SOURCE_DIR}/ddsim/src/Simulator.cpp
    ${CMAKE_SOURCE_DIR}/ddsim/src/QFRSimulator.cpp
    ${CMAKE_SOURCE_DIR}/ddsim/src/GroverSimulator.cpp
    ${CMAKE_SOURCE_DIR}/ddsim/src/ShorSimulator.cpp
    ${CMAKE_SOURCE_DIR}/ddsim/src/ShorFastSimulator.cpp
)

# QFR core sources (include parsers to resolve importReal/importQASM/etc.)
set(QFR_CORE_SOURCES
    ${CMAKE_SOURCE_DIR}/src/QuantumComputation.cpp
    ${CMAKE_SOURCE_DIR}/src/CircuitOptimizer.cpp
    ${CMAKE_SOURCE_DIR}/src/operations/StandardOperation.cpp
    ${CMAKE_SOURCE_DIR}/src/operations/NonUnitaryOperation.cpp
    ${CMAKE_SOURCE_DIR}/src/operations/Operation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/parser_stubs.cpp
)

target_include_directories(qec_benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/operations
    ${CMAKE_SOURCE_DIR}/include/parsers
    ${CMAKE_SOURCE_DIR}/include/algorithms
    ${CMAKE_SOURCE_DIR}/ddsim/include
    ${CMAKE_SOURCE_DIR}/extern/dd_package/include
    ${CMAKE_SOURCE_DIR}/extern/json/single_include
)

target_link_libraries(qec_benchmark PRIVATE
    qec
    JKQ::DDpackage
)

# Compile DDSIM and QFR sources directly into the benchmark
target_sources(qec_benchmark PRIVATE ${DDSIM_SOURCES} ${QFR_CORE_SOURCES})

# Link pthread (required by DDSIM)
find_package(Threads REQUIRED)
target_link_libraries(qec_benchmark PRIVATE Threads::Threads)