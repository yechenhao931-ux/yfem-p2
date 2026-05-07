# tests/CMakeLists.txt
enable_testing()

# 自动查找所有测试文件
file(GLOB TEST_SOURCES "*.cpp")

# 为每个测试文件创建独立的可执行文件
foreach(test_source ${TEST_SOURCES})
    # 获取文件名（不包含扩展名和路径）
    get_filename_component(test_name ${test_source} NAME_WE)
    
    # 创建可执行文件
    add_executable(${test_name} ${test_source})
    
    # 链接必要的库
    target_link_libraries(${test_name} 
        PRIVATE common io graph partition pfem 
        PRIVATE MPI::MPI_CXX
        PRIVATE ${GMSH_LIBRARY}
    )
    

    target_include_directories(${test_name} 
        PRIVATE ${CMAKE_SOURCE_DIR}/include
        PRIVATE ${GMSH_INCLUDE_DIR}
    )


endforeach()