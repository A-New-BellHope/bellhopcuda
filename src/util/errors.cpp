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
#include "../common_run.hpp"
#include "../common_setup.hpp"

namespace bhc {

#define ERRBUFSIZE 1024

void ExternalCommon(bhcInternal *internal, const char *format, va_list *args)
{
    char *buf = new char[ERRBUFSIZE];
    vsnprintf(buf, ERRBUFSIZE, format, *args);
    if(internal->outputCallback == nullptr) {
        printf("%s\n", buf);
    } else {
        internal->outputCallback(buf);
    }
    delete[] buf;
}

[[noreturn]] void ExternalError(bhcInternal *internal, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    ExternalCommon(internal, format, &args);
    va_end(args);
    throw std::runtime_error("See previous line(s) above for error message");
}

void ExternalWarning(bhcInternal *internal, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    ExternalCommon(internal, format, &args);
    va_end(args);
}

static const char *const errorDescriptions[BHC_ERR_MAX] = {
    "BHC_ERR_TEMPLATE: Internal error with template system; should not be possible "
    "for user to cause this error, even with invalid input",
    "BHC_ERR_JOBNUM: Job number for ray run out of range",
    "BHC_ERR_RAYINIT: Source or angle indices out of range",
    "BHC_ERR_BOUNDARY_SEG_CONTAINS_NAN: Boundary depth or normal is not finite",
    "BHC_ERR_BOUNDARY_CONDITION_TYPE: Boundary condition type is not one of the "
    "valid choices",
    "BHC_ERR_INVALID_SSP_TYPE: SSP type is not one of the valid choices",
    "BHC_ERR_OUTSIDE_ALTIMETRY: Top altimetry undefined above where the ray has gone",
    "BHC_ERR_OUTSIDE_BATHYMETRY: Bottom bathymetry undefined below where the ray has "
    "gone",
    "BHC_ERR_OUTSIDE_SSP: Ray has gone outside the box where the soundspeed is defined",
    "BHC_ERR_QUAD_ISEG: SSP segment index in quad SSP has become invalid",
    "BHC_ERR_INVALID_IMAGE_INDEX: Cerveny beam image index has become invalid, "
    "will happen if Nimage is invalid (must be 1, 2, or 3)",
};

static const char *const warningDescriptions[BHC_WARN_MAX] = {
    "BHC_WARN_RAYS_OUTOFMEMORY: Ran out of memory for rays (in ray copy mode); "
    "subsequent rays will be discarded",
    "BHC_WARN_ONERAY_OUTOFMEMORY: Ran out of memory for individual ray(s), "
    "those rays have been truncated",
    "BHC_WARN_UNBOUNDED_BEAM: Cerveny beam has imaginary component of gamma > 0",
    "BHC_WARN_TOO_FEW_BEAMS: Nalpha is too small; there may be gaps between the beams",
    "BHC_WARN_SOURCE_OUTSIDE_BOUNDARIES: Terminating the ray trace because the "
    "source is on or outside the boundaries",
    "BHC_WARN_OUTSIDE_REFLCOEF: Reflection coefficient is undefined where the ray "
    "has gone",
    "BHC_WARN_OUTSIDE_ALTIMETRY: Top altimetry undefined above where the ray has gone",
    "BHC_WARN_OUTSIDE_BATHYMETRY: Bottom bathymetry undefined below where the ray has "
    "gone",
    "BHC_WARN_STEP_NEGATIVE_H: Ray needs to step backwards to reach next boundary; "
    "boundary edge case handling has gone wrong, this ray's results may be unreliable",
    "BHC_WARN_TRIDIAG_H_NEGATIVE: Ray needs to step backwards to cross diagonal "
    "of top/bottom triangle; tri diagonal edge case handling has gone wrong, "
    "this ray's results may be unreliable",
    "BHC_WARN_TRIDIAG_H_GROWING: Tentative step crossed tri diagonal, but stepping "
    "to the diagonal resulted in a larger step; tri diagonal edge case handling "
    "has gone wrong, this ray's results may be unreliable",
    "BHC_WARN_WKB_UNIMPLEMENTED_3D: WKB beamwidth beams unimplemented in BELLHOP3D "
    "(Nx2D or 3D), PickEpsilon results will be nonsense",
    "BHC_WARN_CERVENY_WIDTH_BUGGY: Cerveny beamwidth Cerveny beams are not propery "
    "implemented in BELLHOP(3D), PickEpsilon results will be nonsense",
    "BHC_WARN_INVALID_WIDTH_BUGGY: BELLHOP(3D) does not properly handle Cerveny "
    "beams with an invalid beam width type; PickEpsilon results will be nonsense",
    "BHC_WARN_BEAMTYPE_CARETSPACE: BELLHOP(3D) does not properly handle hat Cartesian "
    "runs defined as '^' or ' ' in PickEpsilon",
    "BHC_WARN_INVALID_TYPE_BUGGY: BELLHOP(3D) does not properly handle beams with "
    "an invalid beam type in PickEpsilon",
    "BHC_WARN_CPCHIP_INVALIDXT: Difference in depth between ray and segment is "
    "extremely large in PCHIP SSP, likely bug or garbage input",
    "BHC_WARN_CPCHIP_INVALIDCCOEF: cCoef is extremely large in PCHIP SSP, likely "
    "bug or garbage input",
    "BHC_WARN_OCEANTORAYX_GAVEUP: Failed to transform Nx2D ray 3D -> 2D -> 3D in a "
    "consistent way, edge case issues may result",
};

void CheckReportErrors(bhcInternal *internal, const ErrState *errState)
{
    uint32_t error     = errState->error.load(STD::memory_order_acquire);
    uint32_t warning   = errState->warning.load(STD::memory_order_acquire);
    uint32_t errCount  = errState->errCount.load(STD::memory_order_acquire);
    uint32_t warnCount = errState->warnCount.load(STD::memory_order_acquire);
    if((error != 0) != (errCount != 0) || (warning != 0) != (warnCount != 0)) {
        ExternalError(
            internal, "Internal error with error counts in error tracking system");
    }
    if(warning != 0) {
        ExternalWarning(
            internal, "%d warning(s) thrown of the following type(s):", warnCount);
        for(int32_t i = 0; i < BHC_WARN_MAX; ++i) {
            if((warning & (1u << i))) {
                ExternalWarning(internal, "%s", warningDescriptions[i]);
                warning &= ~(1u << i);
            }
        }
        if(warning != 0) {
            ExternalError(
                internal,
                "Internal error in error tracking system: unknown warning thrown");
        }
    }
    if(error != 0) {
        ExternalWarning(
            internal, "%d error(s) thrown of the following type(s):", errCount);
        for(int32_t i = 0; i < BHC_ERR_MAX; ++i) {
            if((error & (1u << i))) {
                ExternalWarning(internal, "%s", errorDescriptions[i]);
                error &= ~(1u << i);
            }
        }
        if(error != 0) {
            ExternalError(
                internal,
                "Internal error in error tracking system: unknown error thrown");
        }
        ExternalError(internal, "Raising error(s) reported above to caller");
    }
}

} // namespace bhc
