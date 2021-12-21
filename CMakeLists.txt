cmake_minimum_required(VERSION 3.14)

project(UnionFindCPP)

# read and set union_find version
function(set_union_find_version VERSION_FILE_PATH)
    file(STRINGS ${VERSION_FILE_PATH} VERSION_FILE_STR)
    foreach (LINE IN LISTS VERSION_FILE_STR)
		if("${LINE}" MATCHES "__version__.*")
			set(VERSION_LINE_STR "${LINE}")
		endif()
    endforeach()

    string(REGEX REPLACE "__version__ = \"(.*)\"" "\\1" VERSION_STRING ${VERSION_LINE_STR})
    set(VERSION_STRING ${VERSION_STRING} parent_scope)
endfunction()

set_union_find_version(${PROJECT_SOURCE_DIR}/union_find/_version.py)

message(STATUS "union_find version ${VERSION_STRING}")
set(PROJECT_VERSION ${VERSION_STRING})


if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)
set(LINK_WHAT_YOU_USE TRUE)
set(BUILD_RPATH_USE_ORIGIN TRUE)

option(BUILD_EXECUTABLES "Build C++ executables" OFF)
option(BUILD_TESTS "Build C++ tests" OFF)


add_subdirectory(union_find/src "${CMAKE_BINARY_DIR}")

if(BUILD_TESTS)
	enable_testing()
endif()


if (BUILD_EXECUTABLES)
	add_subdirectory(examples)
endif()

