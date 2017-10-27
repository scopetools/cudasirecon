#ifndef UCSF_MSG_IVE_STANDARDS_H
#define UCSF_MSG_IVE_STANDARDS_H

/* Generated with genus = linux */
/*
 * Copyright(C) 2002
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Sets feature-set macros to ensure that the system headers define
 * the needed functions and types.  Doing so after including one or more
 * system and then including other system headers may cause inconsistent
 * definitions:  always include this header first!  As a consequence of this,
 * do not use constructs (i.e. typedefs of system types) which depend on the
 * inclusion of this header unless the user is made aware of the side effects.
 */








/*
 * GNU/Linux, any platform.  Set _XOPEN_SOURCE to at least 600 for isfinite
 * (also it needs to be at least 500 for pread and pwrite).
 * If _GNU_SOURCE is set, a superset of the X/Open environment is selected.
 */
#ifndef _GNU_SOURCE
#ifdef _XOPEN_SOURCE
#if _XOPEN_SOURCE+0 < 600
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif /* _XOPEN_SOURCE+0 < 600 */
#else
#define _XOPEN_SOURCE 600
#endif /* _XOPEN_SOURCE */
#endif /* _GNU_SOURCE */










#endif /* include guard */
