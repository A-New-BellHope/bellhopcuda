/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
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
#include "paramsmodule.hpp"

namespace bhc {

#if BHC_ENABLE_2D
template class ParamsModule<false, false>;
template class SxSy<false, false>;

#endif
#if BHC_ENABLE_NX2D
template class ParamsModule<true, false>;
template class SxSy<true, false>;

#endif
#if BHC_ENABLE_3D
template class ParamsModule<true, true>;
template class SxSy<true, true>;

#endif

} // namespace bhc
