file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(partition STATIC
    ${SRC_FILES}
)

target_include_directories(partition
    PRIVATE ${CMAKE_SOURCE_DIR}/include
)

# OpenMP：编译期注入 -fopenmp（pragma 生效），并链接运行时库
if(ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(partition PUBLIC OpenMP::OpenMP_CXX)
endif()

    