#pragma once

#ifdef _WIN32
#ifdef STATIC_GARFIELD
#define IMPORTEXPROT 
#else
#if BUILDING_GARFIELD
#define IMPORTEXPROT __declspec(dllexport)
#else
#define IMPORTEXPROT __declspec(dllimport)
#endif
#endif
#else
#define IMPORTEXPROT
#endif
#define GARFIELD_EXTERNAL_SYMBOL IMPORTEXPROT extern

