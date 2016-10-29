#ifndef UCSF_MSG_PLOT2D_FUNCS_H
#define UCSF_MSG_PLOT2D_FUNCS_H

#ifdef __cplusplus
extern "C" {
#endif

int Plot2D_Attach(void);
int Plot2D_OpenFile(const char* filename);
int Plot2D_SetGraphNum(int graph);
int Plot2D_SetScale(const float scale[4]);
int Plot2D_SetTicks(const int numticks[2], const int ticklabs[2]);
int Plot2D_SetTicksExtended(
    const int numticks[2], const int ticklabs[2], const int nplaces[2]
);
int Plot2D_LegendOn(int legend);
int Plot2D_CurveAttr(
    int curve_num,
    int marker_type,
    char char_marker,
    int line,
    int color,
    int allplots
);
int Plot2D_ResizeWindow(int xst, int yst, int width, int height);
int Plot2D_ChangeLabel(int label_num, const char* label_str);
int Plot2D_Update(void);
int Plot2D_Quit(void);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
