project(UnionFindCPP)

set(CMAKE_CXX_STANDARD 20)

set(JSON_URL "https://github.com/nlohmann/json/releases/download/v3.10.4/json.hpp")
set(JSON_DOWNLOAD_DIR "${PROJECT_SOURCE_DIR}/externals/nlohmann")
set(JSON_DOWNLOAD_PATH "${JSON_DOWNLOAD_DIR}/json.hpp")

set(CATCH_URL "https://github.com/catchorg/Catch2/releases/download/v2.11.1/catch.hpp")
set(CATCH_DOWNLOAD_DIR "${PROJECT_SOURCE_DIR}/externals")
set(CATCH_DOWNLOAD_PATH "${CATCH_DOWNLOAD_DIR}/catch.hpp")

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

add_subdirectory(externals/robin-map)

add_library(union_find_cpp_include_dirs INTERFACE)
target_include_directories(union_find_cpp_include_dirs INTERFACE "${PROJECT_SOURCE_DIR}/include")

if(BUILD_TESTS)
	add_subdirectory(tests)
endif()
