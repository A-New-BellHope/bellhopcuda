#pragma once

#ifndef _BHC_INCLUDED_
#error "This file must be included via #include <bhc/bhc.hpp>!"
#endif

#include <cfloat>

////////////////////////////////////////////////////////////////////////////////
//Real types
////////////////////////////////////////////////////////////////////////////////

#ifdef USE_FLOATS
using real = float;
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
#define REAL_MINPOS FLT_MIN
#define REAL_PI ((float)M_PI)
#define REAL_REL_SNAP (1.0e-5f)
// Must be below abs(bit_cast<float>(0xFEFEFEFEu) == -1.69e38f)
#define DEBUG_LARGEVAL (1.0e30f)
// #define DEBUG_LARGEVAL (1.0e15f)
#else
using real = double;
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define REAL_MINPOS DBL_MIN
#define REAL_PI M_PI
#define REAL_REL_SNAP (1.0e-6)
//Must be below abs(bit_cast<double>(0xFEFEFEFEFEFEFEFEull) == -5.31e303)
#define DEBUG_LARGEVAL (1.0e250)
//#define DEBUG_LARGEVAL (1.0e15)
#endif

// BELLHOP uses mostly normal REAL literals, which are float despite most of the
// program using REAL*8 (double). It occasionally uses double-precision literals
// (look like 1.6D-9). For values exactly expressible in both types, e.g. 0.0,
// 1.0, 2.0, 37.0, 0.375, it doesn't matter which type the literal is--except
// that when running in float mode, double-precision literals may cause the
// entire expression to be promoted to double, causing extremely slow double-
// precision instructions to be emitted on CUDA. However, for values not
// expressable as float, e.g. 0.1, the literal type changes the result:
// double d = bar(); float f = foo();
// assert(d * 0.1 == d * 0.1f); //will fail for all d except 0, inf, etc.
// assert((float)(f * 0.1) == f * 0.1f); //will fail for some f
//
// "Real literal"
#ifdef USE_FLOATS
#define RL(a) (a##f)
#else
#define RL(a) a
#endif
// "Float literal"--This macro is not actually needed, values could just be
// always written as e.g. 0.1f, but this way, every floating-point literal
// should have one or the other macro on it, making it easier to spot errors.
#define FL(a) (a##f)

////////////////////////////////////////////////////////////////////////////////
//Real vectors
////////////////////////////////////////////////////////////////////////////////

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/common.hpp>
#include <glm/geometric.hpp>

using vec2 = glm::vec<2, real, glm::defaultp>;
using vec3 = glm::vec<3, real, glm::defaultp>;
