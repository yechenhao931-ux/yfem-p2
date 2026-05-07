file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(io STATIC
    ${SRC_FILES}
)

target_include_directories(io
    PRIVATE ${CMAKE_SOURCE_DIR}/include
    PRIVATE ${GMSH_INCLUDE_DIR}
)

target_link_libraries(io
    PRIVATE ${GMSH_LIBRARY}
)
    

