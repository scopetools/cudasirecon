#ifndef UCSF_MSG_IOMENU_H
#define UCSF_MSG_IOMENU_H

#include "IWInclude.h"
#include "iveqd.h"

#define IOMENU_REG_CTRL_ALL 0
#define IOMENU_REG_CTRL_5D_SIZE 1
#define IOMENU_REG_CTRL_3D_SIZE 2
#define IOMENU_REG_CTRL_XY_ONLY 3
#define IOMENU_REG_CTRL_4D_NO_Z 4
#define IOMENU_REG_CTRL_4D_NO_W 5
#define IOMENU_REG_CTRL_4D_NO_T 6


typedef enum IomenuFileType {
    IOMENU_PRIMARY_IMAGE = 1,
    IOMENU_OTHER_INPUT = 2,
    IOMENU_OUTPUT_DEFAULT = 3,
    IOMENU_OTHER_OUTPUT = 4
} IomenuFileType;

typedef struct IomenuFile {
    char filename[256];
    char filter[80];
    char label[80];
    IomenuFileType type;
} IomenuFile;

typedef enum IomenuParameterType {
    IOMENU_INTEGER = 1,
    IOMENU_REAL = 2,
    IOMENU_ACTIVE_WHEN_NONZERO = 3,
    IOMENU_ACTIVE_WHEN_ZERO = -3,
    IOMENU_STRING = 4
} IomenuParameterType;

typedef struct IomenuParameter {
    char label[80];
    union {
	struct { int* values; int count; } i;
	struct { float* values; int count; } r;
	struct { char** values; int count; } c;
	int* pb;
    } v;
    IomenuParameterType type;
} IomenuParameter;

typedef struct Iomenu {
    IomenuFile* files;
    IomenuParameter* parameters;
    void* closure;
    int file_count;
    int parameter_count;
    char xyztw_size[4 + IW_MAX_WAVE][16];
    int iparm, nxyztw[5], ixyztw[4 + IW_MAX_WAVE][3], mode;
    float dmin, dmax, dmean, delxyz[3];
} Iomenu;

struct IomenuCallbackImpl;
typedef struct IomenuCallbackImpl* IomenuCallback;

#ifdef __cplusplus
extern "C" {
#endif

IomenuCallback iomenu_add_callback_int_property(
    const char* name, void (*callback)(void*), void* closure
);

IVEQDMember iomenu_get_current_queue(void);

int iomenu_get_int_property(const char* name, int* p_i);

int iomenu_set_int_property(const char* name, int i);

void iomenu_set_queue_callback(void (*func)(Iomenu*));

void iomenu_set_region_controls(int style);

void iomenu(
    const char* title,
    const char* submenu_title,
    const char* app_name,
    Iomenu* iom,
    void (*menu_proc)(Iomenu* iom),
    void (*special_menu_proc)(Iomenu* iom),
    void (*handle_file_proc)(Iomenu* iom),
    void (*handle_exit)(Iomenu* iom)
);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
