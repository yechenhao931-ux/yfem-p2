#ifndef TYPE_HPP
#define TYPE_HPP
#include <cfloat>
#include <cstdint>



#ifndef USE_32
#define int_t int64_t
#define real_t double
#else
#define int_t int32_t
#define real_t float

#endif

#endif