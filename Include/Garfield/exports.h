#pragma once

#ifdef _WIN32
#ifdef STATIC_GARFIELD
#define GARFIELD_IMPORTEXPROT 
#define MAGBOLTZ_IMPORTEXPROT 
#else
#if BUILDING_GARFIELD
#define GARFIELD_IMPORTEXPROT __declspec(dllexport)
#else
#define GARFIELD_IMPORTEXPROT __declspec(dllimport)
#endif
#if BUILDING_MAGBOLTZ
#define MAGBOLTZ_IMPORTEXPROT __declspec(dllexport)
#else
#define MAGBOLTZ_IMPORTEXPROT __declspec(dllimport)
#endif
#endif
#else
#define GARFIELD_IMPORTEXPROT
#define MAGBOLTZ_IMPORTEXPROT 
#endif
#define GARFIELD_EXTERNAL_SYMBOL GARFIELD_IMPORTEXPROT extern
#define MAGBOLTZ_EXTERNAL_SYMBOL MAGBOLTZ_IMPORTEXPROT extern
