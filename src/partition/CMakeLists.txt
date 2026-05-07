file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(partition STATIC
    ${SRC_FILES}
)

target_include_directories(partition
    PRIVATE ${CMAKE_SOURCE_DIR}/include
)

    