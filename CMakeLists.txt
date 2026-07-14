cmake_minimum_required(VERSION 3.20)
project(MOM CXX)

# C++ standard --------------------------------------------------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Source/execution charset --------------------------------------------------
# The MOM sources are UTF-8 encoded (they contain e.g. "C₂H₂", "✓", "→" in
# comments/strings and non-ASCII intermediate identifiers).  MSVC otherwise
# decodes source files using the host system codepage (e.g. 936/GBK), which
# corrupts these multibyte sequences and produces spurious C2065/C3688 errors.
# /utf-8 forces both the source- and execution-charset to UTF-8.
add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/utf-8>)

# Build options --------------------------------------------------------------
option(MOM_BUILD_SHARED "Build libMOM as a shared library (.so/.dylib)" OFF)
option(MOM_BUILD_TESTS "Build the standalone test executables" ON)
option(MOM_INSTALL "Generate install rules" ON)
option(MOM_USE_DICTIONARY "Enable dictionary support" OFF)

# Eigen ----------------------------------------------------------------------
if(DEFINED EIGEN3_DIR)
    message(STATUS "Using Eigen3 from: ${EIGEN3_DIR}")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_DIR})
else()
    find_package(Eigen3 3.4 QUIET NO_MODULE)
    if(NOT Eigen3_FOUND)
        message(STATUS "Eigen3 not found via find_package, fetching from gitlab")
        include(FetchContent)
        FetchContent_Declare(
            Eigen3
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 3.4.0
            GIT_SHALLOW TRUE
        )
        FetchContent_MakeAvailable(Eigen3)
    endif()
endif()

# std::expected (C++23) is used only in ParseConfig<> template bodies, which
# live at the bottom of each variant's .tpp and are compiled only when
# MOM_USE_DICTIONARY=ON.  The core library targets C++20; the dictionary
# build upgrades to C++23 below.

# Helper: apply common settings to a MOM target ------------------------------
function(mom_configure_target tgt access)
    target_include_directories(${tgt} ${access}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>
        $<INSTALL_INTERFACE:include>)
    target_link_libraries(${tgt} ${access} Eigen3::Eigen)
    target_compile_features(${tgt} ${access} cxx_std_20)
endfunction()

# MOM_headers: header-only INTERFACE target ---------------------------------
# Use this when you need to instantiate MOM templates with custom Thermo types.
add_library(MOM_headers INTERFACE)
mom_configure_target(MOM_headers INTERFACE)

# MOM: compiled library (static or shared) ----------------------------------
#
# Core sources — compiled unconditionally; no dependency on OpenSMOKE++.
set(_MOM_SOURCES
    src/ThreeEquations/ThreeEquations_BasicThermoData.cpp
    src/HMOM/HMOM_BasicThermoData.cpp
    src/BrookesMoss/BrookesMoss_BasicThermoData.cpp
    src/MetalOxide/MetalOxide_BasicThermoData.cpp
    src/Utilities/OutputFileColumns.cpp
    src/MOM/ManifestReader.cpp
    src/MOM/MechanismAdapter.cpp
    src/MOM/MomContext.cpp
)

# Dictionary support ---------------------------------------------------------
#
# Grammar sources define DefineRules() for each variant's OpenSMOKE++ grammar.
# They include Grammar headers which pull in OpenSMOKE++ headers, so they must
# only be compiled when MOM_USE_DICTIONARY=ON and DICTIONARY_DIR is available.
if(MOM_USE_DICTIONARY)
    if(NOT DEFINED DICTIONARY_DIR)
        message(FATAL_ERROR
            "MOM_USE_DICTIONARY=ON requires DICTIONARY_DIR to be set.\n"
            "Example: -DDICTIONARY_DIR=/path/to/dictionary")
    endif()

    list(APPEND _MOM_SOURCES
        src/ThreeEquations/ThreeEquations_Grammar.cpp
        src/HMOM/HMOM_Grammar.cpp
        src/BrookesMoss/BrookesMoss_Grammar.cpp
        src/MetalOxide/MetalOxide_Grammar.cpp
    )

    message(STATUS "MOM: Dictionary support enabled")
    message(STATUS "     Dictionary path : ${DICTIONARY_DIR}")
endif()

# Choose library type
set(_MOM_LIB_TYPE STATIC)
if(MOM_BUILD_SHARED)
    set(_MOM_LIB_TYPE SHARED)
endif()

add_library(MOM ${_MOM_LIB_TYPE} ${_MOM_SOURCES})
mom_configure_target(MOM PUBLIC)

# This definition tells the headers that a pre-compiled version exists
target_compile_definitions(MOM INTERFACE MOM_COMPILED_LIBRARY)

# Dictionary linking and compile definition
if(MOM_USE_DICTIONARY)
    target_include_directories(MOM PUBLIC ${DICTIONARY_DIR}/include)
    target_link_libraries(MOM PUBLIC "${DICTIONARY_DIR}/build/libdictionary_parser.a")

    # Expose MOM_USE_DICTIONARY so variant headers activate ParseConfig<> and Grammar includes.
    target_compile_definitions(MOM PUBLIC MOM_USE_DICTIONARY)

    # ParseConfig<> uses std::expected — upgrade to C++23 for dictionary builds.
    target_compile_features(MOM PUBLIC cxx_std_23)
endif()

# Compiler warnings (only for our own sources)
target_compile_options(MOM PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)

# Tests ----------------------------------------------------------------------
if(MOM_BUILD_TESTS)
    enable_testing()

    add_executable(test_loop_benchmark tests/test_loop_benchmark.cpp)
    target_link_libraries(test_loop_benchmark PRIVATE MOM)
    target_compile_options(test_loop_benchmark PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME LoopBenchmark COMMAND test_loop_benchmark)

    add_executable(test_verify_all_variants tests/test_verify_all_variants.cpp)
    target_link_libraries(test_verify_all_variants PRIVATE MOM)
    target_compile_options(test_verify_all_variants PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME VerifyAllVariants COMMAND test_verify_all_variants)

    add_executable(reporter_usage_example tests/reporter_usage_example.cpp)
    target_link_libraries(reporter_usage_example PRIVATE MOM)
    target_compile_options(reporter_usage_example PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME ReporterUsageExample COMMAND reporter_usage_example)

    add_executable(test_oh_fixture tests/test_oh_fixture.cpp)
    target_link_libraries(test_oh_fixture PRIVATE MOM)
    target_compile_options(test_oh_fixture PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME OHFixture COMMAND test_oh_fixture)

    # G1Closure + G1Negative: Python-based closure verification on the JSON
    # produced by OHFixture. Python3 is a HARD requirement -- if not found,
    # configuration fails rather than silently dropping to 7 tests.
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    add_test(
        NAME G1Closure
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/g1_closure_analysis.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/oh_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(G1Closure PROPERTIES DEPENDS OHFixture)

    # G1Negative: proves G1Closure fails on wrong stoichiometry (red/green proof).
    add_test(
        NAME G1Negative
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/g1_negative_test.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/oh_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(G1Negative PROPERTIES DEPENDS OHFixture)

    # --- M1 Condensation fixture (condensation-only vertical slice) --------
    add_executable(test_condensation_fixture tests/test_condensation_fixture.cpp)
    target_link_libraries(test_condensation_fixture PRIVATE MOM)
    target_compile_options(test_condensation_fixture PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME CondFixture COMMAND test_condensation_fixture)

    # CondClosure + CondNegative: Python-based closure verification on the JSON
    # produced by CondFixture.
    add_test(
        NAME CondClosure
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/cond_closure_analysis.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/condensation_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(CondClosure PROPERTIES DEPENDS CondFixture)

    # CondNegative: proves CondClosure fails when omega_gas[C2H2] is zeroed
    # (simulating the pre-fix bug where PAH gas consumption was gated by
    # nucleation_model_).
    add_test(
        NAME CondNegative
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/cond_negative_test.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/condensation_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(CondNegative PROPERTIES DEPENDS CondFixture)

    # --- M1 Nucleation fixture (nucleation-only vertical slice) ------------
    add_executable(test_nucleation_fixture tests/test_nucleation_fixture.cpp)
    target_link_libraries(test_nucleation_fixture PRIVATE MOM)
    target_compile_options(test_nucleation_fixture PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME NucFixture COMMAND test_nucleation_fixture)

    add_test(
        NAME NucClosure
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/nuc_closure_analysis.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/nucleation_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(NucClosure PROPERTIES DEPENDS NucFixture)

    add_test(
        NAME NucNegative
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/nuc_negative_test.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/nucleation_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(NucNegative PROPERTIES DEPENDS NucFixture)

    # --- M1 Surface growth fixture (HACA, tanh threshold states) ----------
    add_executable(test_surface_growth_fixture tests/test_surface_growth_fixture.cpp)
    target_link_libraries(test_surface_growth_fixture PRIVATE MOM)
    target_compile_options(test_surface_growth_fixture PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME SGFixture COMMAND test_surface_growth_fixture)

    # SGClosure + SGNegative: Python-based closure verification.
    add_test(
        NAME SGClosure
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/sg_closure_analysis.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/surface_growth_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(SGClosure PROPERTIES DEPENDS SGFixture)

    add_test(
        NAME SGNegative
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/sg_negative_test.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/surface_growth_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(SGNegative PROPERTIES DEPENDS SGFixture)

    # --- M1 O2 Oxidation fixture ------------------------------------------------
    add_executable(test_o2_oxidation_fixture tests/test_o2_oxidation_fixture.cpp)
    target_link_libraries(test_o2_oxidation_fixture PRIVATE MOM)
    target_compile_options(test_o2_oxidation_fixture PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME O2OxFixture COMMAND test_o2_oxidation_fixture)

    # O2OxClosure + O2OxNegative: Python-based closure verification.
    add_test(
        NAME O2OxClosure
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/o2ox_closure_analysis.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/o2_oxidation_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(O2OxClosure PROPERTIES DEPENDS O2OxFixture)

    add_test(
        NAME O2OxNegative
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/o2ox_negative_test.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/o2_oxidation_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(O2OxNegative PROPERTIES DEPENDS O2OxFixture)

    # --- M1 OH-oxidation fixture (extends G1 with tanh threshold states) ---
    add_executable(test_oh_oxidation_fixture tests/test_oh_oxidation_fixture.cpp)
    target_link_libraries(test_oh_oxidation_fixture PRIVATE MOM)
    target_compile_options(test_oh_oxidation_fixture PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME OHOxFixture COMMAND test_oh_oxidation_fixture)

    add_test(
        NAME OHOxClosure
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/ohox_closure_analysis.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/oh_oxidation_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(OHOxClosure PROPERTIES DEPENDS OHOxFixture)

    add_test(
        NAME OHOxNegative
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/ohox_negative_test.py
                --input ${CMAKE_CURRENT_BINARY_DIR}/oh_oxidation_fixture_results.json
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(OHOxNegative PROPERTIES DEPENDS OHOxFixture)

    add_executable(test_manifest_reader tests/test_manifest_reader.cpp)
    target_link_libraries(test_manifest_reader PRIVATE MOM)
    target_compile_options(test_manifest_reader PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME ManifestReader COMMAND test_manifest_reader)

    add_executable(test_mechanism_adapter tests/test_mechanism_adapter.cpp)
    target_link_libraries(test_mechanism_adapter PRIVATE MOM)
    target_compile_options(test_mechanism_adapter PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME MechanismAdapter COMMAND test_mechanism_adapter)

    add_executable(test_mom_context tests/test_mom_context.cpp)
    target_link_libraries(test_mom_context PRIVATE MOM)
    target_compile_options(test_mom_context PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
    add_test(NAME MomContext COMMAND test_mom_context)

    if(MOM_USE_DICTIONARY)
        add_executable(test_dictionary_setup tests/test_dictionary_setup.cpp)
        target_link_libraries(test_dictionary_setup PRIVATE MOM)
        target_compile_options(test_dictionary_setup PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
        add_test(NAME DictionarySetup COMMAND test_dictionary_setup)
    endif()

endif()

# Install rules --------------------------------------------------------------
if(MOM_INSTALL)
    include(GNUInstallDirs)
    install(TARGETS MOM MOM_headers
        EXPORT MOMTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING
            PATTERN "*.hpp"
            PATTERN "*.tpp"
            PATTERN "*.h"
    )

    install(EXPORT MOMTargets
        FILE MOMTargets.cmake
        NAMESPACE MOM::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOM
    )

    include(CMakePackageConfigHelpers)
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/MOMConfigVersion.cmake"
        VERSION 0.1.0
        COMPATIBILITY AnyNewerVersion
    )

    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MOMConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/MOMConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOM
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/MOMConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/MOMConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOM
    )
endif()

# Configuration summary ------------------------------------------------------
message(STATUS "")
message(STATUS "MOM library configuration:")
message(STATUS "  Library type     : ${_MOM_LIB_TYPE}")
message(STATUS "  C++ standard     : C++20 (C++23 when MOM_USE_DICTIONARY=ON)")
message(STATUS "  Dictionary       : ${MOM_USE_DICTIONARY}")
if(MOM_USE_DICTIONARY)
    message(STATUS "    Dictionary dir : ${DICTIONARY_DIR}")
endif()
message(STATUS "  Build tests      : ${MOM_BUILD_TESTS}")
message(STATUS "  Install rules    : ${MOM_INSTALL}")
message(STATUS "")
