file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(pfem STATIC
    ${SRC_FILES}
)

target_include_directories(pfem
    PRIVATE ${CMAKE_SOURCE_DIR}/include
)