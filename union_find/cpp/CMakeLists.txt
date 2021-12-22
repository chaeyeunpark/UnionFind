cmake_minimum_required(VERSION 3.14)
project(UnionFindCPP)

set(CMAKE_CXX_STANDARD 20)

option(CLANG_TIDY OFF)

if (CLANG_TIDY)
    set(CMAKE_CXX_CLANG_TIDY clang-tidy --header-filter="${PROJECT_SOURCE_DIR}/include/.*" --extra-arg=-std=c++20)
endif()

# Download nlohmann/json
set(JSON_URL "https://github.com/nlohmann/json/releases/download/v3.10.4/json.hpp")
set(JSON_DOWNLOAD_DIR "${PROJECT_SOURCE_DIR}/externals/nlohmann")
set(JSON_DOWNLOAD_PATH "${JSON_DOWNLOAD_DIR}/json.hpp")

if (NOT EXISTS "${JSON_DOWNLOAD_PATH}")
    message("Downloading nlohmann/json.")
    file(MAKE_DIRECTORY "${JSON_DOWNLOAD_DIR}")
    file(DOWNLOAD "${JSON_URL}" "${JSON_DOWNLOAD_PATH}" STATUS JSON_DOWNLOAD_STATUS)

    # Separate the returned status code, and error message.
    list(GET JSON_DOWNLOAD_STATUS 0 STATUS_CODE)
    list(GET JSON_DOWNLOAD_STATUS 1 ERROR_MESSAGE)

    # Check if download was successful.
    if(${STATUS_CODE} EQUAL 0)
        message(STATUS "Download completed successfully!")
    else()
        message(FATAL_ERROR "Error occurred during download: ${ERROR_MESSAGE}")
    endif()
endif()

# Download Catch
set(CATCH_URL "https://github.com/catchorg/Catch2/releases/download/v2.13.7/catch.hpp")
set(CATCH_DOWNLOAD_DIR "${PROJECT_SOURCE_DIR}/externals")
set(CATCH_DOWNLOAD_PATH "${CATCH_DOWNLOAD_DIR}/catch.hpp")

if (NOT EXISTS "${CATCH_DOWNLOAD_PATH}")
    message("Downloading Catch2.")
    file(DOWNLOAD "${CATCH_URL}" "${CATCH_DOWNLOAD_PATH}" STATUS CATCH_DOWNLOAD_STATUS)

    # Separate the returned status code, and error message.
    list(GET CATCH_DOWNLOAD_STATUS 0 STATUS_CODE)
    list(GET CATCH_DOWNLOAD_STATUS 1 ERROR_MESSAGE)

    # Check if download was successful.
    if(${STATUS_CODE} EQUAL 0)
        message(STATUS "Download completed successfully!")
    else()
        message(FATAL_ERROR "Error occurred during download: ${ERROR_MESSAGE}")
    endif()

endif()

# Add robin-map
add_subdirectory(externals/robin-map)

add_library(union_find_cpp_dependency INTERFACE)
target_include_directories(union_find_cpp_dependency INTERFACE "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/externals")
target_link_libraries(union_find_cpp_dependency INTERFACE tsl::robin_map)
target_sources(union_find_cpp_dependency INTERFACE "src/utility.cpp")

if(BUILD_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif()
