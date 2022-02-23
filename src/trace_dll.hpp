//Joe Snider
//12/14/2021
//
//Interface for a test dll (not general)

#pragma once

#include <string>

#ifdef BELLHOPCXX_EXPORTS
#define BELLHOPCXX_API __declspec(dllexport)
#else
#define BELLHOPCXX_API __declspec(dllimport)
#endif

extern "C" int BELLHOPCXX_API RunBellhop(std::string FileRoot, void * result);
