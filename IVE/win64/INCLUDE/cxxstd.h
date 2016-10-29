#ifndef UCSF_MSG_CXXSTD_H
#define UCSF_MSG_CXXSTD_H

// @(#) $Id: cxxstd.h,v 1.2 2008/04/15 03:09:29 eric Exp $
// $Name:  $
// Defines macros so code can conditionally compile stuff that depends
// on whether or not the standard C++ headers (no .h) and namespace are
// available.  If HAS_STDCXX is set, the compiler, at least, provides the
// standard C++ headers which are not the standard C headers wrapped in the
// std namespace.  If HAS_STDCXX_C is set, the compiler provides the
// standard C headers wrapped in the std namespace.

#ifdef __sgi
#ifdef _STANDARD_C_PLUS_PLUS
#define HAS_STDCXX
#define HAS_STDCXX_C
#endif
#else
#define HAS_STDCXX
// Version 3.2-3 of pgCC does not include the C headers wrapped in the std
// namespace.  Do not know whether other versions of the Portland Group
// compiler do.
#ifndef __PGI
#define HAS_STDCXX_C
#endif
#endif

#ifdef HAS_STDCXX
#define STDNAMESPACE(x) std::x
#else
#define STDNAMESPACE(x) ::x
#endif

#ifdef HAS_STDCXX_C
#define STDNAMESPACE_C(x) std::x
#else
#define STDNAMESPACE_C(x) ::x
#endif

#endif /* include guard */
