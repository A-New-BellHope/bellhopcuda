//Joe Snider
//12/14/2021
//
//Interface for a test shared library (not general)

#pragma once

#include <string>

/// <summary>
/// Run Bellhop on the given file.
/// TODO: intended primarily for testing and only supports rays.
/// </summary>
/// <param name="cFileRoot">The raw file name 
/// (may use ./xxx if it's in the current directory)</param>
/// <param name="result">will be reinterpreted as a ray2dPt</param>
/// <returns>-1 on failure, or number of rays attempted</returns>
extern "C" int BELLHOPCXX_API RunBellhop(const char* cFileRoot, void * result);
