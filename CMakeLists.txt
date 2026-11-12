cmake_minimum_required(VERSION 3.20)
project(predtmr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PREDTMR_BUILD_WASM "Build the WebAssembly bindings (phase B)" OFF)
option(PREDTMR_BUILD_TESTS "Build the unit tests" ON)

# The pure kernel is shared by every target (native CLI, tests, WASM module).
add_library(predtmr_core STATIC
    src/core/indicators.cpp
    src/core/predictor.cpp
    src/parser/parser.cpp
    src/parser/sequence.cpp
)
target_include_directories(predtmr_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(NOT EMSCRIPTEN)
    target_compile_options(predtmr_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(EMSCRIPTEN)
    # Emscripten toolchain: only the WASM module is produced. Native executables
    # (CLI, tests) are irrelevant in this configuration.
    if(PREDTMR_BUILD_WASM)
        add_executable(predtmr_wasm src/wasm/bindings.cpp)
        target_link_libraries(predtmr_wasm PRIVATE predtmr_core)
        target_link_options(predtmr_wasm PRIVATE
            --bind
            -O3
            -sMODULARIZE=1
            -sEXPORT_NAME=createPredTMR
            -sALLOW_MEMORY_GROWTH=1
            -sENVIRONMENT=web
            -sEXPORT_ES6=0
            -sFILESYSTEM=0
        )
        # Emit predtmr.js (JS glue) + predtmr.wasm directly into the web portal.
        set_target_properties(predtmr_wasm PROPERTIES
            OUTPUT_NAME predtmr
            SUFFIX ".js"
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/web
        )
    else()
        message(STATUS "Emscripten toolchain detected but PREDTMR_BUILD_WASM=OFF; no module built.")
    endif()
else()
    add_executable(predtmr src/cli/main.cpp)
    target_link_libraries(predtmr PRIVATE predtmr_core)
    target_compile_options(predtmr PRIVATE -Wall -Wextra -Wpedantic)

    if(PREDTMR_BUILD_TESTS)
        enable_testing()
        add_executable(predtmr_tests
            tests/main.cpp
            tests/test_core.cpp
            tests/test_parser.cpp
        )
        target_link_libraries(predtmr_tests PRIVATE predtmr_core)
        target_include_directories(predtmr_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
        target_compile_options(predtmr_tests PRIVATE -Wall -Wextra -Wpedantic)
        add_test(NAME predtmr_tests COMMAND predtmr_tests)
    endif()
endif()
