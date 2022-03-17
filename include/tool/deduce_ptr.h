#pragma once
#include "macro.h"
#include <cstddef>
#include <type_traits>

namespace tinker {
/// \ingroup rc
/// A helper class to get the traits of the given 1D or 2D pointer type.
///
/// Example:
/// | PTR        | DeducePtr<PTR>::type | DeducePtr<PTR>::n |
/// |:----------:|:--------------------:|:-----------------:|
/// | float*     | float                | 1                 |
/// | int (*)[3] | int                  | 3                 |
///
/// \tparam PTR Must be a pointer type.
template <class PTR>
struct DeducePtr;

template <class T>
struct DeducePtr<T*>
{
   typedef T type;
   static constexpr size_t n = 1;
};

template <class T, size_t N>
struct DeducePtr<T (*)[N]>
{
   typedef T type;
   static constexpr size_t n = N;
};
}
