file(GLOB children RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

foreach(child ${children})
    # 2. 判断是否为目录，并且该目录下是否存在 CMakeLists.txt 文件
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${child} AND 
       EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${child}/CMakeLists.txt)
        
        add_subdirectory(${child})
        message(STATUS ${child})
    endif()
endforeach()