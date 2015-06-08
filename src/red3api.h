#ifndef RED3API_H
#define RED3API_H

#if _WIN32 || _MSC_VER

#ifdef red3_core_EXPORTS
#define RED3_API __declspec(dllexport)
#else
#define RED3_API __declspec(dllimport)
#endif
#else
#define RED3_API
#endif

#ifdef _MSC_VER
#pragma warning(disable: 4251)
#endif

#endif
