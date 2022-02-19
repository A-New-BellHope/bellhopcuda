#pragma once

#ifndef _BHC_INCLUDED_
#error "This file must be included via #include <bhc/bhc.hpp>!"
#endif

#include <complex>

////////////////////////////////////////////////////////////////////////////////
//Complex types
////////////////////////////////////////////////////////////////////////////////

using cpx = STD::complex<real>;
using cpxf = STD::complex<float>;
#define J cpx(RL(0.0), RL(1.0))
HOST_DEVICE constexpr inline cpxf Cpx2Cpxf(const cpx &c){
	return cpxf((float)c.real(), (float)c.imag());
}
HOST_DEVICE constexpr inline cpx Cpxf2Cpx(const cpxf &c){
	return cpx((real)c.real(), (real)c.imag());
}

//CUDA::std::cpx<double> does not like operators being applied with float
//literals, due to template type deduction issues.
#ifndef USE_FLOATS
HOST_DEVICE constexpr inline cpx operator-(float a, const cpx &b){
    return cpx(a - b.real(), -b.imag());
}
HOST_DEVICE constexpr inline cpx operator*(const cpx &a, float b){
    return cpx(a.real() * b, a.imag() * b);
}
HOST_DEVICE constexpr inline cpx operator*(float a, const cpx &b){
    return cpx(a * b.real(), a * b.imag());
}
HOST_DEVICE constexpr inline cpx operator/(const cpx &a, float b){
    return cpx(a.real() / b, a.imag() / b);
}
HOST_DEVICE inline cpx operator/(float a, const cpx &b){
    return (double)a / b;
}
#endif
