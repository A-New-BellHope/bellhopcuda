/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2022 The Regents of the University of California
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

#ifndef _BHC_INCLUDED_
#error "This file must be included via #include <bhc/bhc.hpp>!"
#endif

#include <cfloat>
#include <complex>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/mat2x2.hpp>

namespace bhc {

#ifdef BHC_USE_FLOATS
using real = float;
#else
using real = double;
#endif

using vec2   = glm::vec<2, real, glm::defaultp>;
using vec3   = glm::vec<3, real, glm::defaultp>;
using int2   = glm::vec<2, int32_t, glm::defaultp>;
using mat2x2 = glm::mat<2, 2, real, glm::defaultp>;

using cpx  = STD::complex<real>;
using cpxf = STD::complex<float>;

} // namespace bhc
