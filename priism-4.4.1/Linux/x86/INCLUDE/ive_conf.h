#ifndef UCSF_MSG_IVE_CONF_H
#define UCSF_MSG_IVE_CONF_H

/* Generated with genus = linux , species = x86 , and IVE_SHM64 =    */
/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Set the macros controlling platform-specific behavior and conditional
 * compilation.  Define typedefs for fixed-size types.
 */


/*
 * Clear any outside settings.
 */
#undef IVECONF_XOPEN
#undef IVECONF_BSD
#undef IVECONF_WIN32
#undef IVECONF_PTHREAD
#undef IVECONF_DEFINES_ETIMEDOUT
#undef IVECONF_DEFINES_ELOOP
#undef IVECONF_DEFINES_EOVERFLOW
#undef IVECONF_DEFINES_SEMUN
#undef IVECONF_DEFINES_MAP_FAILED
#undef IVECONF_POSIX_NAMED_SEM_IN_FS
#undef IVECONF_NO_SEM_GETVALUE
#undef IVECONF_UNIMPL_SEM_GETVALUE
#undef IVECONF_SHM64
#undef IVECONF_SHMMUTEX_USE_PTHREAD
#undef IVECONF_SHMMUTEX_USE_SYSVSEM
#undef IVECONF_SHMMUTEX_USE_FCNTL
#undef IVECONF_SHMHOLD_USE_SYSVSEM
#undef IVECONF_SHMHOLD_USE_FCNTL
#undef IVECONF_SHMRWLK_USE_PTHREAD
#undef IVECONF_SHMRWLK_USE_SYSVSEM
#undef IVECONF_SHMRWLK_USE_FCNTL
#undef IVECONF_SHMSEM_USE_POSIXUNNAMED
#undef IVECONF_SHMSEM_USE_POSIXNAMED
#undef IVECONF_SHMSEM_USE_SYSVSEM
#undef IVECONF_SHMSEM_USE_CYGWIN_HACK
#undef IVECONF_HAS_O_SYNC
#undef IVECONF_HAS_SIGINFO
#undef IVECONF_HAS_VSNPRINTF
#undef IVECONF_USE_ISNANF
#undef IVECONF_FINITE_IN_IEEEFP_H
#undef IVECONF_INTTYPES_SYSTYPES_CLASH
#undef IVECONF_HAS_DLOPEN
#undef IVECONF_HAS_RTLD_LOCAL
#undef IVECONF_HAS_NSMODULE
#undef IVECONF_NO_INTTYPES
#undef IVECONF_CLOSE_BEFORE_REMOVE
#undef IVECONF_HAS_IPV6
#undef IVECONF_HAS_UNSIGNED_TIME_T
#undef IVECONF_NPROC_SYSCONF_NPROCESSORS_ONLN
#undef IVECONF_NPROC_SYSCONF_NPROC_ONLN
#undef IVECONF_NPROC_SYSCTL
#undef IVECONF_HAS_STATVFS
#undef IVECONF_HAS_STATVFS64
#undef IVECONF_HAS_STATFS
#undef IVECONF_STATFS_REQUIRES_PARAM_H
#undef IVECONF_STATFS_REQUIRES_MOUNT_H
#undef IVECONF_STATFS_REQUIRES_VFS_H
#undef IVECONF_STATFS_REQUIRES_STATFS_H
#undef IVECONF_STATFS_AVAIL_FIELD
#undef IVECONF_HAS_NSGETEXECUTABLEPATH
#undef IVECONF_HAS_PROCSELF
#undef IVECONF_HAS_PATH_DEFPATH
#undef IVECONF_DEFAULT_PATH
#undef IVECONF_USE_STAT64










/*
 * GNU/Linux, any platform.  Not fully XOpen, but close enough.
 */
#define IVECONF_XOPEN 1
#define IVECONF_DEFINES_ETIMEDOUT 1
#define IVECONF_DEFINES_ELOOP 1
#define IVECONF_DEFINES_EOVERFLOW 1
/*
 * Older versions of Linux (for instance Debian potato, 2.0.36 kernel) define
 * semun in sys/sem.h; newer versions do not but set _SEM_SEMUN_UNDEFINED to
 * indicate that they do not.
 */
#if (defined(__GNU_LIBRARY__) && ! defined(_SEM_SEMUN_UNDEFINED))
#define IVECONF_DEFINES_SEMUN 1
#endif

#define IVECONF_DEFINES_MAP_FAILED 1
#define IVECONF_SHMMUTEX_USE_FCNTL 1
#define IVECONF_SHMHOLD_USE_FCNTL 1
#define IVECONF_SHMRWLK_USE_FCNTL 1
#define IVECONF_SHMSEM_USE_SYSVSEM 1
#define IVECONF_HAS_O_SYNC
#define IVECONF_HAS_SIGINFO
#define IVECONF_HAS_VSNPRINTF
#define IVECONF_HAS_DLOPEN
#define IVECONF_HAS_IPV6
#define IVECONF_NPROC_SYSCONF_NPROCESSORS_ONLN
/*
 * statfs is present in glibc 2.1.3 and may be available in earlier versions.
 * Later versions (glibc 2.3.6 for instance) have statvfs as well and
 * deprecate statfs.  Use statfs for backward compatibility.
 */
#define IVECONF_HAS_STATFS
#define IVECONF_STATFS_REQUIRES_VFS_H
#define IVECONF_STATFS_AVAIL_FIELD(buf) (buf).f_bavail
#define IVECONF_HAS_PROCSELF
/*
 * This is what the manpage execlp() says in the default path for
 * Debian 5.0 and CentOS 5.5.
 */
#define IVECONF_DEFAULT_PATH ":/bin:/usr/bin"
/* This is a hack for large files on 32-bit Linux systems. */
#ifdef _LARGEFILE64_SOURCE
#define IVECONF_USE_STAT64
#endif
typedef float float32_t;
typedef double float64_t;










#if defined(IVECONF_XOPEN) || defined(IVECONF_BSD)
typedef int IVEDesc;
#ifndef IVECONF_NO_INTTYPES
#ifdef IVECONF_INTTYPES_SYSTYPES_CLASH
/*
 * A kludge for IRIX 5.3 where both inttypes.h and sys/types.h define int8_t.
 * int16_t, and int32_t
 */
#include <sys/types.h>
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
#else
#include <inttypes.h>
#endif
#endif
/*
 * If available, use pthreads as the threading package.
 */
#include <unistd.h>
#ifdef _POSIX_THREADS
#define IVECONF_PTHREAD
#endif /* _POSIX_THREADS */
#endif /* IVECONF_XOPEN || IVECONF_BSD */

#endif /* include guard */
