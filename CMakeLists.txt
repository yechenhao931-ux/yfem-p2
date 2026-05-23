cmake_minimum_required(VERSION 3.24)
project(fem CXX)

# --- 1. 定义开关选项 ---
# 默认设置为 OFF，使用时通过 cmake -DENABLE_GMSH=ON .. 开启
option(ENABLE_GMSH "Enable Gmsh support" OFF)
# OpenMP 并行（默认 ON；用 -DENABLE_OPENMP=OFF 可构建串行基线做对比）
option(ENABLE_OPENMP "Enable OpenMP parallelism" ON)

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


# --- 2a. 条件查找 OpenMP ---
# --- 2a. OpenMP support (Linux/macOS) ---
if(ENABLE_OPENMP)

    # macOS special handling
    if(APPLE)

        # Homebrew LLVM
        set(HOMEBREW_LLVM "/usr/local/opt/llvm")

        # Use Homebrew clang instead of AppleClang
        set(CMAKE_C_COMPILER "${HOMEBREW_LLVM}/bin/clang")
        set(CMAKE_CXX_COMPILER "${HOMEBREW_LLVM}/bin/clang++")

        # libomp paths
        set(OpenMP_ROOT "/usr/local/opt/libomp")

        # include/lib
        include_directories("${OpenMP_ROOT}/include")
        link_directories("${OpenMP_ROOT}/lib")

        # OpenMP flags for macOS
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")

        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_LIB_NAMES "omp")

        set(OpenMP_omp_LIBRARY "${OpenMP_ROOT}/lib/libomp.dylib")

    endif()

    find_package(OpenMP REQUIRED)

    if(OpenMP_CXX_FOUND)

        message(STATUS "OpenMP found")
        message(STATUS "OpenMP CXX flags: ${OpenMP_CXX_FLAGS}")

    else()

        message(FATAL_ERROR "OpenMP requested but not found")

    endif()

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

if(ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(yfem PRIVATE OpenMP::OpenMP_CXX)
endif()

# --- 3. 针对 target 注入宏和库 ---
if(ENABLE_GMSH)
    target_include_directories(yfem PUBLIC ${GMSH_INCLUDE_DIR})
    target_link_libraries(yfem PUBLIC ${GMSH_LIBRARY})
    # 注入 C/C++ 宏定义
    target_compile_definitions(yfem PUBLIC USE_GMSH)
endif()

add_subdirectory(src)
add_subdirectory(test)
