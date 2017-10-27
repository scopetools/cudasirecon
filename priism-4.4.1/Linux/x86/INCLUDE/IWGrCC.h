#ifndef UCSF_MSG_IWGRCC_H
#define UCSF_MSG_IWGRCC_H

#include "IWApiConstants.h"
#include "IWShmConstants.h"

#define MAX_OBJ_FILE_NAME_SIZE    48
#define MAX_OBJ_NAME_SIZE         8 

#define IWOBJ_COMPOSITE           0x00000001
#define IWOBJ_BOX                 0x00000002
#define IWOBJ_POLY                0x00000004
#define IWOBJ_CUBE                0x00000008
#define IWOBJ_AREA                0x00000010
#define IWOBJ_VECTOR              0x00000020
#define IWOBJ_CIRCLE              0x00000040
#define IWOBJ_SPHERE              0x00000080
#define IWOBJ_PTS                 0x00000100
#define IWOBJ_HARMONIC            0x00000200
#define IWOBJ_SIMPLE              0x000003FE

/* obj's data type */
#define IWOBJ_REAL_SP             0x00000001
#define IWOBJ_DATA_SP             0x00000002

/* obj's display method */
#define IWOBJ_OUTLINE             0x00000001
#define IWOBJ_SOLID               0x00000002

/* load obj */
#define IWOBJ_UNKNOWN             -1
#define IWOBJ_CANCEL              0x00000000
#define IWOBJ_TRUNCATE            0x00000001
#define IWOBJ_APPEND              0x00000002

#define IWOBJ_RANGE_ALL           -1


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef void* IWGrHeaderID;
typedef void* IWGrObjID;

int IWRtNumObjLists(void);
char** IWRtObjListNames(void);
IWGrHeaderID IWInitObjList(const char* name);
int IWInitObjHdr(
    IWGrHeaderID headerID,
    const int nxyzwt[],
    const float delt[],
    const float tilt[],
    const float orig[],
    const int xyzst[]
);
int IWRtObjHdr(
    IWGrHeaderID headerID,
    int nxyzwt[],
    float delt[],
    float tilt[],
    float orig[],
    int xyzst[]
);
int IWAttachObjToWin(IWGrHeaderID headerID, int istr);
int IWDetachObjFromWin(IWGrHeaderID headerID, int istr);
int IWCloseObjList(IWGrHeaderID headerID);
int IWSaveObjList(IWGrHeaderID headerID, const char* filename);
int IWLoadObjList(IWGrHeaderID headerID, const char* filename, int mode);
int IWRtObjVersion(const char* file_name);
IWGrObjID IWGroupObj(
    IWGrHeaderID headerID, const IWGrObjID obj_list[], int nobj
);
IWGrObjID IWBoxObj(
    IWGrHeaderID headerID,
    const IW_POINT box[2],
    int iwave,
    int itime,
    int type,
    int color
);
IWGrObjID IWPolyObj(
    IWGrHeaderID headerID,
    const IW_POINT poly[],
    int npt,
    int iwave,
    int itime,
    int type,
    int color
);
IWGrObjID IWAreaObj(
    IWGrHeaderID headerID,
    const IW_POINT list[],
    int npt,
    int iwave,
    int itime,
    int type,
    int color
);
IWGrObjID IWPtsObj(
    IWGrHeaderID headerID,
    const IW_POINT list[],
    int npt,
    int iwave,
    int itime,
    int type,
    int color
);
IWGrObjID IWCircleObj(
    IWGrHeaderID headerID,
    IW_POINT center,
    int iwave,
    int itime,
    float radius,
    int type,
    int color
);
IWGrObjID IWSphereObj(
    IWGrHeaderID headerID,
    IW_POINT center,
    int iwave,
    int itime,
    float radius,
    int type,
    int color
);
int IWAddMemObj(
    IWGrHeaderID headerID, IWGrObjID com_objID, IWGrObjID mem_objID
);
int IWDelMemObj(
    IWGrHeaderID headerID, IWGrObjID com_objID, IWGrObjID mem_objID
);
IWGrObjID IWCopyObj(IWGrHeaderID headerID, IWGrObjID objID);
int IWDelObj(IWGrHeaderID headerID, IWGrObjID objID, int delete_child);
float IWRtObjVol(IWGrHeaderID headerID, IWGrObjID objID);
float IWRtObjInt(IWGrHeaderID headerID, IWGrObjID objID, int istr);
int IWRtObjCM(
    IWGrHeaderID headerID, IWGrObjID objID, int istr, int d_type, IW_POINT* pt
);
int IWRtObjGC(
    IWGrHeaderID headerID, IWGrObjID objID, int d_type, IW_POINT* pt
);
int IWRtObjs(
    IWGrHeaderID headerID,
    int type_mask,
    int itype,
    const float range[10],
    IWGrObjID** obj_list,
    int* nobj
);
IWGrObjID IWRtObjSelected(
    IWGrHeaderID headerID, int itype, IW_POINT pt, int wave, int time
);
int IWIsObjSelected(
    IWGrHeaderID headerID,
    IWGrObjID objID,
    int itype,
    IW_POINT pt,
    int wave,
    int time
);
IWGrObjID IWRtObjPtSelected(
    IWGrHeaderID headerID,
    int itype,
    IW_POINT pt,
    int wave,
    int time,
    int max_dist,
    int* pt_num
);
int IWRtMemObjs(
    IWGrHeaderID headerID, IWGrObjID objID, IWGrObjID** obj_list, int* nobj
);
int IWRtObjType(IWGrHeaderID headerID, IWGrObjID objID);
int IWRtObjPts(
    IWGrHeaderID headerID,
    IWGrObjID objID,
    int itype,
    IW_POINT** pts, 
    int* pt_num,
    int* iwave,
    int* itime
);
int IWRtObjRange(
    IWGrHeaderID headerID, IWGrObjID objID, int itype, float range[10]
);
int IWRtObjLnCount(IWGrHeaderID headerID, IWGrObjID objID);
int IWRtObjName(IWGrHeaderID headerID, IWGrObjID objID, char* name);
int IWAlObjDis(IWGrHeaderID headerID, IWGrObjID objID, int iflag);
int IWAlObjDisAttr(
    IWGrHeaderID headerID,
    IWGrObjID objID,
    int color,
    int pattern,
    int thickness,
    int method
);
int IWAlObjPt(IWGrHeaderID headerID, IWGrObjID objID, int pt_num, IW_POINT pt);
int IWAlObjName(IWGrHeaderID headerID, IWGrObjID objID, const char* name);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
