#pragma once

// this is the header that clcpp files will use
// it has to be included after autogen.h which defines real
// but the #include can't go in here cuz luajit ffi doesn't like them

#if 1
#include "math.luajit.h"
#endif

#if 0
#include "Tensor/Vector.h"
using real3 = Tensor::vec<real, 3>;
using real4 = Tensor::vec<real, 4>;
#endif
