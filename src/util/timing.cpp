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
#include "../common_setup.hpp"

#if defined(_WIN32) || defined(_WIN64)
// #include <windows.h>
#else
// sched_setscheduler():
#include <sched.h>
#endif

namespace bhc {

void SetupThread()
{
#ifdef BHC_USE_HIGH_PRIORITY_THREADS
#if defined(_WIN32) || defined(_WIN64)
    // std::cout << "Warning, not changing thread priority because on Windows\n";
#else
    struct sched_param sp = {.sched_priority = sched_get_priority_max(SCHED_RR)};
    if(sched_setscheduler(0, SCHED_RR, &sp) < 0) {
        std::cout << "Could not change thread priority (probably not root)\n";
    }
#endif
#endif
}

} // namespace bhc
