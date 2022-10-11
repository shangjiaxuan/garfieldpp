#pragma once

#ifdef _WIN32
#ifdef STATIC_GARFIELD
#define GARFIELD_IMPORTEXPROT 
#define MEGABOLTZ_IMPORTEXPROT 
#else
#if BUILDING_GARFIELD
#define GARFIELD_IMPORTEXPROT __declspec(dllexport)
#else
#define GARFIELD_IMPORTEXPROT __declspec(dllimport)
#endif
#if BUILDING_MEGABOLTZ
#define MEGABOLTZ_IMPORTEXPROT __declspec(dllexport)
#else
#define MEGABOLTZ_IMPORTEXPROT __declspec(dllimport)
#endif
#endif
#else
#define GARFIELD_IMPORTEXPROT
#define MEGABOLTZ_IMPORTEXPROT 
#endif
#define GARFIELD_EXTERNAL_SYMBOL GARFIELD_IMPORTEXPROT extern
#define MEGABOLTZ_EXTERNAL_SYMBOL MEGABOLTZ_IMPORTEXPROT extern
