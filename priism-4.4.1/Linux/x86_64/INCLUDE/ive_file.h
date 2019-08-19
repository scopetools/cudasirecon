#ifndef UCSF_MSG_IVE_FILE_H
#define UCSF_MSG_IVE_FILE_H

/*
 * Copyright(C) 2002
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Wrap OS-specific calls for functions that operate on a file name or
 * return a file name.
 */

#include <stddef.h>  /* size_t */

enum {
    IVE_FILE_CHECK_EXIST = 1 << 0,
    IVE_FILE_CHECK_DIR = 1 << 1,
    IVE_FILE_CHECK_ACCESS = 1 << 2,
    IVE_FILE_CHECK_ALL = (1 << 3) - 1
};

enum {
    IVE_FILE_INDETERM = 1 << 0,
    IVE_FILE_EXIST = 1 << 1,
    IVE_FILE_DIR = 1 << 2,
    IVE_FILE_READ = 1 << 3,
    IVE_FILE_WRITE = 1 << 4,
    IVE_FILE_EXEC = 1 << 5
};


#ifdef __cplusplus
extern "C" {
#endif

int IVECheckFileAttributes(const char* name, int mask);
size_t IVEGetExecutableName(const char* argv0, size_t length, char* name);
int IVEGetFreeSpace(const char* path, int units, int* p_status);
void IVEPermitExecution(const char* name, int* p_status);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
