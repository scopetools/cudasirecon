#ifndef _GLwDrawAP_h
#define _GLwDrawAP_h

#ifdef __GLX_MOTIF
#include <X11/GLw/GLwMDrawA.h>
#else
#include <X11/GLw/GLwDrawA.h>
#endif

typedef struct _GLwDrawingAreaClassPart
{
    caddr_t extension;
} GLwDrawingAreaClassPart;

#ifdef __GLX_MOTIF
typedef struct _GLwMDrawingAreaClassRec {
    CoreClassPart		core_class;
    XmPrimitiveClassPart	primitive_class;
    GLwDrawingAreaClassPart	glwDrawingArea_class;
} GLwMDrawingAreaClassRec;

extern GLwMDrawingAreaClassRec glwMDrawingAreaClassRec;

#else /* not __GLX_MOTIF */

typedef struct _GLwDrawingAreaClassRec {
    CoreClassPart		core_class;
    GLwDrawingAreaClassPart	glwDrawingArea_class;
} GLwDrawingAreaClassRec;

extern GLwDrawingAreaClassRec glwDrawingAreaClassRec;
#endif /* __GLX_MOTIF */

typedef struct {
    /* resources */
    int *		attribList;
    XVisualInfo *	visualInfo;
    Boolean		myList;		/* TRUE if we malloced the attribList*/
    Boolean		myVisual;	/* TRUE if we created the visualInfo*/
    Boolean		installColormap;
    Boolean		allocateBackground;
    Boolean		allocateOtherColors;
    Boolean		installBackground;
    XtCallbackList	ginitCallback;
    XtCallbackList	resizeCallback;
    XtCallbackList	exposeCallback;
    XtCallbackList	inputCallback;
    /* specific attributes; add as we get new attributes */
    int			bufferSize;
    int			level;
    Boolean		rgba;
    Boolean		doublebuffer;
    Boolean		stereo;
    int			auxBuffers;
    int			redSize;
    int			greenSize;
    int			blueSize;
    int			alphaSize;
    int			depthSize;
    int			stencilSize;
    int			accumRedSize;
    int			accumGreenSize;
    int			accumBlueSize;
    int			accumAlphaSize;
} GLwDrawingAreaPart;

#ifdef __GLX_MOTIF
typedef struct _GLwMDrawingAreaRec {
    CorePart		core;
    XmPrimitivePart	primitive;
    GLwDrawingAreaPart	glwDrawingArea;
} GLwMDrawingAreaRec;
#else /* not __GLX_MOTIF */
typedef struct _GLwDrawingAreaRec {
    CorePart		core;
    GLwDrawingAreaPart	glwDrawingArea;
} GLwDrawingAreaRec;
#endif /* __GLX_MOTIF */

#endif /* _GLwDrawP_h */
