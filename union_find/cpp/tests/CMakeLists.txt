project(UnionFindCPP_tests)

add_executable(test_LatticeFromParity "test_LatticeFromParity.cpp")
target_link_libraries(test_LatticeFromParity union_find_cpp_dependency)
add_test(NAME test_LatticeFromParity
         COMMAND test_LatticeFromParity)
