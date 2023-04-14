/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP(3D) underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

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

////////////////////////////////////////////////////////////////////////////////
// Select which standard library
////////////////////////////////////////////////////////////////////////////////

#ifndef STD
#define BHC_UNDEF_STD_AFTER 1
#define STD std
#endif

////////////////////////////////////////////////////////////////////////////////
// Shared library setup
////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
#if defined(BHC_DLL_EXPORT)
#define BHC_API __declspec(dllexport)
#elif defined(BHC_DLL_IMPORT)
#define BHC_API __declspec(dllimport)
#else
#define BHC_API
#endif
#else
#define BHC_API
#endif
