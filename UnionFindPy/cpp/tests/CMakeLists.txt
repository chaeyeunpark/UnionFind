project(UnionFindCPP_tests)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(test_LatticeFromParity "test_LatticeFromParity.cpp")
target_link_libraries(test_LatticeFromParity union_find_cpp_dependency Eigen3::Eigen)
add_test(NAME test_LatticeFromParity
         COMMAND test_LatticeFromParity)
