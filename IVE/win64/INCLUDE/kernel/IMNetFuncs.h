#ifndef UCSF_MSG_IMNETFUNCS_H
#define UCSF_MSG_IMNETFUNCS_H

/* @(#) $Id: IMNetFuncs.h,v 1.2 2003/04/03 22:35:09 eric Exp $ */
/* $Name:  $ */

#include "ive_conf.h"
#include <sys/types.h>  /* ssize_t */


/* Defines the elements common to all messages */
typedef struct remote_comm_hdr {
    uint32_t type;
    uint32_t size;
} remote_comm_hdr;

/* Define possible values for the type field in IM messages. */
#define  NET_IMOPEN       101
#define  NET_IMCLOSE      102
#define  NET_IMRDHDR      103
#define  NET_IMPUTHDR     104
#define  NET_IMWRHDR      105
#define  NET_IMALEXHDR    106
#define  NET_IMRDPAS      107    
#define  NET_IMWRPAS      108   
#define  NET_IMALFMT      109

/* Defines the elements common to all IM message structures. */
typedef struct remote_comm_pars {
    uint32_t type;
    uint32_t size;
    uint32_t istream;
} remote_comm_pars;

typedef struct remote_open_pars {
    uint32_t type;
    uint32_t size;
    uint32_t istream;
    uint32_t attr;
    uint32_t name_length;
} remote_open_pars;

typedef struct remote_wrhdr_pars {
    uint32_t type;
    uint32_t size;
    uint32_t istream;
    float32_t mean;
} remote_wrhdr_pars;

typedef struct remote_alexhdr_pars {
    uint32_t type;
    uint32_t size;
    uint32_t istream;
    int32_t section;
    int32_t space_group;
    union {
	struct { uint32_t nbytes; } all;
	struct { uint16_t nint, nfloat; } part;
    } var;
} remote_alexhdr_pars;

typedef struct remote_rdwr_pars {
    uint32_t type;
    uint32_t size;
    uint32_t istream;
    uint32_t ix;
    uint32_t nx;
    uint32_t ny;
    uint32_t cur_y;
    uint32_t cur_section;
    uint32_t cur_resolution;
    uint16_t display;
    uint16_t mode;
} remote_rdwr_pars;

typedef struct remote_alfmt_pars {
    uint32_t type;
    uint32_t size;
    uint32_t istream;
    uint32_t format;
} remote_alfmt_pars;

/* Define possible values for the attr field in the open message structure. */
#define NET_IMOPEN_NEW 0
#define NET_IMOPEN_OLD 1
#define NET_IMOPEN_READONLY 2
#define NET_IMOPEN_SCRATCH 3


#ifdef __cplusplus
extern "C" {
#endif

int IMget_port_file(char* filename, const char* prefix, int max_length);
int IMhandle_imsubs_request(int sock, remote_comm_pars* p_pars, int swab);
ssize_t IMread_socket(int sock, void* buffer, size_t sz);
int IMinitiate_byte_order_determination(int sock, int* p_different_order);
int IMrespond_byte_order_determination(int sock, int* p_different_order);

#ifdef __cplusplus
}
#endif


#endif /* include guard */
