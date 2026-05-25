file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(partition STATIC
    ${SRC_FILES}
)

target_include_directories(partition
    PRIVATE ${CMAKE_SOURCE_DIR}/include
)

# 原生 METIS：仅 metispart.cpp 需要 metis.h 与 USE_METIS 宏（PRIVATE）；
# 库以 PUBLIC 方式链接 libmetis，使依赖 partition 的目标自动获得 -lmetis。
# 消费者无需该宏，运行期可用 MetisAvailable() 判断。
if(ENABLE_METIS AND METIS_LIBRARY)
    target_include_directories(partition PRIVATE ${METIS_INCLUDE_DIR})
    target_link_libraries(partition PUBLIC ${METIS_LIBRARY})
    target_compile_definitions(partition PRIVATE USE_METIS)
endif()

    