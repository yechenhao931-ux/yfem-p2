#ifndef MACRO_HPP
#define MACRO_HPP
#include <iostream>
#include <cassert>
#ifdef YFEM_DEBUG
#define YFEM_PRINT(msg) std::cout << msg << std::endl; 
// 2. 带自定义消息的断言
#define YFEM_ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "Assertion failed: " << #condition \
                     << "\nMessage: " << message \
                     << "\nFile: " << __FILE__ \
                     << "\nLine: " << __LINE__ \
                     << "\nFunction: " << __func__ << std::endl; \
            std::abort(); \
        } \
    } while(0)
#else
#define YFEM_PRINT(msg) ((void)0)

// 2. 带自定义消息的断言
#define YFEM_ASSERT(condition, message) ((void)0)
#endif


// 11. 编译时断言
#define STATIC_ASSERT(condition, message) \
    static_assert((condition), message)
#endif