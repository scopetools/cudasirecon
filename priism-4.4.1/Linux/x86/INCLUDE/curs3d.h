#ifndef UCSF_MSG_CURS3D_H
#define UCSF_MSG_CURS3D_H

#include <sys/types.h>  /* pid_t */

typedef struct cursor_info {
    pid_t pid;
    int istr;
    int nxyz[3];
    float del[3];
    float orig[3];
    int depth;
    int z_mode;
    float xyz[3];
    float stretch;
    float tlt[3];
    float rot;
    float img_mat[3][3];
    float xshift[2];
    int ids[12];
    int app_modify;
    int display;
    int num_prog;
    int timing_flag;
} CURSOR_INFO, *CURSOR_INFO_PTR;

#ifdef __cplusplus
extern "C" {
#endif

int IWEnable3DCursor(void);
int IWDisable3DCursor(void);
int IWUpdate3DCursorMatrix(void);
int IWCheck3DCursorDepth(void);
int IWSignal3DCursorWait(int wait);
int IWRt3DCursorPos(float rtn_pos[3]);
int IWSet3DCursorDisplay(int flag);
int IWSet3DCursorDepth(const float position[3]);
int IWProjectedPosition(float proj_pos[2], const float old_point[3]);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
