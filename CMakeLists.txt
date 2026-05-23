cmake_minimum_required(VERSION 3.24)
project(fem CXX)

# --- 1. 定义开关选项 ---
# 默认设置为 OFF，使用时通过 cmake -DENABLE_GMSH=ON .. 开启
option(ENABLE_GMSH "Enable Gmsh support" OFF)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 默认构建类型
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "DEBUG")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Enable Gmsh: ${ENABLE_GMSH}")

# 根据构建类型设置不同编译选项和宏定义
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    add_compile_options(-g -O0 -Wall -Wextra)
    add_compile_definitions(YFEM_DEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
    add_compile_options(-O3 -DNDEBUG -Wall)
    add_compile_definitions(YFEM_REL)
endif()


# --- 2. 条件查找 Gmsh ---
if(ENABLE_GMSH)
    find_path(GMSH_INCLUDE_DIR gmsh.h PATHS ${CMAKE_SOURCE_DIR}/thirdparty/gmsh-sdk/include)
    find_library(GMSH_LIBRARY NAMES gmsh PATHS ${CMAKE_SOURCE_DIR}/thirdparty/gmsh-sdk/lib)
    
    if(NOT GMSH_INCLUDE_DIR OR NOT GMSH_LIBRARY)
        message(FATAL_ERROR "Gmsh enabled but not found!")
    endif()
endif()

add_executable(yfem
    src/main.cpp
)

target_include_directories(yfem
    PRIVATE ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(yfem 
    PRIVATE common graph io partition pfem
)

# --- 3. 针对 target 注入宏和库 ---
if(ENABLE_GMSH)
    target_include_directories(yfem PUBLIC ${GMSH_INCLUDE_DIR})
    target_link_libraries(yfem PUBLIC ${GMSH_LIBRARY})
    # 注入 C/C++ 宏定义
    target_compile_definitions(yfem PUBLIC USE_GMSH)
endif()

add_subdirectory(src)
add_subdirectory(test)
