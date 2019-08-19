/*
 * As the mouse passes through image windows, this application will overlay
 * the image window with a graph displaying the intensity profile along
 * the horizontal line where the mouse is located.  Uses the IWStartOverlay
 * and IWEndOverlay calls if overlay planes are supported to improve
 * reponsiveness; otherwise uses IWGrnFastGr, IWDisFastGr, and IWEndFastGr.
 */

#include "WMInclude.h"
#include "IWInclude.h"
#include <math.h>      /* atan, fabs, floor, sqrt */
#include <stdlib.h>    /* calloc, exit, free */


typedef enum FastGraphics {
    FAST_NONE, FAST_OVERLAY, FAST_MONITOR
} FastGraphics;

typedef struct IntensityBounds {
    float min;
    float max;
} IntensityBounds;

typedef struct DatasetScale {
    IntensityBounds bounds[IW_MAX_WAVE];
    int initialized;
} DatasetScale;


static const int g_in_stream = 1;
static int g_open_window = 0;
static FastGraphics g_fast_mode = FAST_NONE;
static int g_need_clear = 0;
/* Remember the scaling used for each image window. */
static DatasetScale* g_scales = 0;


static void compute_complex_graph(
    int npts,
    const float* raw,
    float rmin,
    float rmax,
    float smin,
    float smax,
    int iz,
    int mapping,
    IW_POINT* graph
);
static void compute_graph(
    int npts,
    const float* raw,
    float rmin,
    float rmax,
    float smin,
    float smax,
    int iz,
    IW_POINT* graph
);
static void find_visible_region(
    int in_stream,
    int data_nx,
    int data_ny,
    int wave_index,
    int lower_left[2],
    int upper_right[2]
);
static int handle_win_change(int unused, XEvent* event);
static int handle_win_event(int type, XEvent* event);
static int quit(int* unused);
static int quit_wrapper(void);
static void swap_xy(int npts, IW_POINT* graph);
static void update_display(
    int in_stream,
    const IntensityBounds* p_bounds,
    const float coord[3],
    int iw,
    int it,
    FastGraphics fast_mode,
    int* p_need_clear
);


int
main(int argc, char* argv[])
{
#define MSG_COUNT 5
#define MSG_LENGTH 30
    char msgs[5][MSG_LENGTH] = {
	"Displays horizontal profile",
	"while mouse is in an image",
	"window.  Special keys:",
	"0 or Page Down compress plot",
	"9 or Page Up expand plot"
    };
    int imsg;

    WMInit("Cross Profile");
    WMSetLoc(500, 500);
    for (imsg = 0; imsg < MSG_COUNT; ++imsg) {
	WMAddText(msgs[imsg], MSG_LENGTH - 1, 0);
	WMAttachRightSide();
	WMNewRow();
    }
    WMAddFuncButton("Quit", quit, NULL, 0, 0);
    WMAttachRightSide();

    IMAlPrt(0);
    WMSetExitFunction(quit_wrapper);
    WMAddEventHandler(
	IW_ALL_WINDOWS, EnterWindowMask, handle_win_event, EnterWindowMask
    );
    WMAddEventHandler(
	IW_ALL_WINDOWS, LeaveWindowMask, handle_win_event, LeaveWindowMask
    );
    WMProcDisplayChange(IW_ALL_WINDOWS, handle_win_change, 0);
    WMEnableIWLEvent();

    g_scales = calloc(IWRtMaxWin(), sizeof(DatasetScale));
    if (g_scales == 0) {
	quit(0);
	return 1;
    }

    WMDisplay();
    WMAppMainLoop();

    return 1;
}


/*
 * Performs as compute_graph below, but raw is complex so a mapping of
 * the complex value to a real value is performed as specified by
 * mapping (allowed options are IW_COMPLEX_AMPLITUDE, IW_COMPLEX_PHASE,
 * IW_COMPLEX_REAL, and IW_COMPLEX_IMAGINARY).
 */
static
void
compute_complex_graph(
    int npts,
    const float* raw,
    float rmin,
    float rmax,
    float smin,
    float smax,
    int iz,
    int mapping,
    IW_POINT* graph
)
{
    float scale;
    int ipt, iraw;

    if (rmax - rmin == 0.0f) {
	rmin -= 1.0f;
    }
    scale = (smax - smin) / (rmax - rmin);

    switch (mapping) {
    case IW_COMPLEX_AMPLITUDE:
	for (ipt = 0, iraw = 0; ipt < npts; ++ipt, iraw += 2) {
	    float amplitude = sqrt(
		raw[iraw] * raw[iraw] + raw[iraw + 1] * raw[iraw + 1]
	    );

	    graph[ipt].x = ipt + 0.5f;
	    graph[ipt].y = (amplitude - rmin) * scale + smin;
	    graph[ipt].z = iz;
	}
	break;

    case IW_COMPLEX_PHASE:
	for (ipt = 0, iraw = 0; ipt < npts; ++ipt, iraw += 2) {
	    float phase = (raw[iraw] == 0.0f && raw[iraw + 1] == 0.0f) ?
		0.0f : atan2(raw[iraw + 1], raw[iraw]);

	    graph[ipt].x = ipt + 0.5f;
	    graph[ipt].y = (phase - rmin) * scale + smin;
	    graph[ipt].z = iz;
	}
	break;

    case IW_COMPLEX_REAL:
    case IW_COMPLEX_IMAGINARY:
	for (ipt = 0, iraw = (mapping == IW_COMPLEX_REAL) ? 0 : 1;
	     ipt < npts;
	     ++ipt, iraw += 2) {
	    graph[ipt].x = ipt + 0.5f;
	    graph[ipt].y = (raw[iraw] - rmin) * scale + smin;
	    graph[ipt].z = iz;
	}
	break;
    }
}


/*
 * Fills graph with points whose x coordinates vary linearly from
 * 0.5 to npts - 0.5 and whose y coordinates are calculated from the values
 * in raw with a linear scaling applied such that a raw value of rmin has a
 * scaled value of smin and a raw value of rmax has a scaled value of smax.
 * The z coordinates for the points are set to iz.
 */
static
void
compute_graph(
    int npts,
    const float* raw,
    float rmin,
    float rmax,
    float smin,
    float smax,
    int iz,
    IW_POINT* graph
)
{
    float scale;
    int ipt;

    /*
     * If drawing a profile through integer-valued data, this is ok, but
     * is questionable for handling rmax == rmin for arbitrary
     * floating-point data.
     */
    if (rmax - rmin == 0.0f) {
	rmin -= 1.0f;
    }
    scale = (smax - smin) / (rmax - rmin);

    for (ipt = 0; ipt < npts; ++ipt) {
	graph[ipt].x = ipt + 0.5f;
	graph[ipt].y = (raw[ipt] - rmin) * scale + smin;
	graph[ipt].z = iz;
    }
}


/*
 * Determines the lower left and upper right corners (in data coordinates)
 * for the visible portion of an data_nx by data_ny section in wavelength
 * wave_index.
 */
static
void
find_visible_region(
    int in_stream,
    int data_nx,
    int data_ny,
    int wave_index,
    int lower_left[2],
    int upper_right[2]
)
{
    int montage_n_horiz, montage_n_vert;
    int win_x, win_y, win_width, win_height;
    int data_coords[4], win_coords[4];

    IWRtMulDis(in_stream, &montage_n_horiz, &montage_n_vert);
    IWRtDisArea(in_stream, &win_x, &win_y, &win_width, &win_height);

    /*
     * Probe to determine where the full extent of the data would fall in
     * window coordinates.
     */
    data_coords[0] = 0;
    data_coords[1] = 1;
    data_coords[2] = data_nx;
    data_coords[3] = data_ny;
    IWDataToWin(in_stream, wave_index, win_coords, data_coords, 2);

    /* Now clip to the bounds of the image displayed */
    if (win_coords[0] < 0) {
	win_coords[0] = 0;
    }
    if (win_coords[1] < 0) {
	win_coords[1] = 0;
    }
    if (win_coords[2] > win_width / montage_n_horiz) {
	win_coords[2] = win_width / montage_n_horiz; 
    }
    if (win_coords[3] > win_height / montage_n_vert) {
	win_coords[3] = win_height / montage_n_vert;
    }

    /* Map back to data coordinates. */
    IWWinToData(in_stream, wave_index, win_coords, data_coords, 2);

    lower_left[0] = data_coords[0];
    lower_left[1] = data_coords[1];
    upper_right[0] = data_coords[2];
    upper_right[1] = data_coords[3];
}


static
int
handle_win_change(int unused, XEvent* event)
{
    int changed_window = IWIsIWLEvt(event);

    if (changed_window == g_open_window) {
	switch (event->xclient.data.s[IW_CM_CHANGE_OR_SYN]) {
	case IW_WIN_KILLED:
	    (void) handle_win_event(LeaveWindowMask, NULL);
	    break;

	default:
	    if (g_fast_mode == FAST_MONITOR) {
		/*
		 * In this mode, a recorded image is redisplayed and the
		 * the profile is drawn on top.  With the window change,
		 * the previously recorded image may no longer be appropriate
		 * so record a new one.
		 */
		IWBgnFastGr(g_in_stream);
	    }
	}
    }

    if (event->xclient.data.s[IW_CM_CHANGE_OR_SYN] == IW_WIN_KILLED &&
	changed_window > 0) {
	g_scales[changed_window - 1].initialized = 0;
    }

    return 1;
}


static
int
handle_win_event(int type, XEvent* event)
{
    int window;

    switch (type) {
    case EnterWindowMask:
	window = IWIsIWLEvt(event);
	(void) handle_win_event(LeaveWindowMask, 0);
	if (window > 0 && IWAttachWin(g_in_stream, window, "ro") != IW_ERROR) {
	    DatasetScale* p_scale = &g_scales[window - 1];
	    int nxyz[3], mxyz[3], mode;
	    float dmin, dmax, dmean;

	    IMRdHdr(g_in_stream, nxyz, mxyz, &mode, &dmin, &dmax, &dmean);
	    g_open_window = window;
	    if (IWRtNumOVColors() > 0) {
		g_fast_mode = FAST_OVERLAY;
	    } else {
		g_fast_mode = FAST_MONITOR;
		IWBgnFastGr(g_in_stream);
	    }

	    if (! p_scale->initialized) {
		int nz, nw, nt, intlv;
		int iw;

		IMRtZWT(g_in_stream, &nz, &nw, &nt, &intlv);
		for (iw = 0; iw < nw; ++iw) {
		    IMRtWavMM(
			g_in_stream,
			iw,
			&(p_scale->bounds[iw].min),
			&(p_scale->bounds[iw].max)
		    );
		    if (p_scale->bounds[iw].min > p_scale->bounds[iw].max) {
			float tmp = p_scale->bounds[iw].min;

			p_scale->bounds[iw].min = p_scale->bounds[iw].max;
			p_scale->bounds[iw].max = tmp;
		    }
		}
	    }
	}
	WMAddEventHandler(
	    g_in_stream, PointerMotionMask, handle_win_event, PointerMotionMask
	);
	WMAddEventHandler(
	    g_in_stream, KeyReleaseMask, handle_win_event, KeyReleaseMask
	);
	break;

    case LeaveWindowMask:
	if (g_open_window != 0) {
	    if (g_fast_mode == FAST_OVERLAY) {
		if (g_need_clear) {
		    IWClearOverlay(g_in_stream);
		    g_need_clear = 0;
		}
	    } else if (g_fast_mode == FAST_MONITOR) {
		IWEndFastGr(g_in_stream);
		IWDisplay(g_in_stream);
		g_need_clear = 0;
	    }
	    g_fast_mode = FAST_NONE;

	    WMCancelEventHandler(
		g_in_stream, PointerMotionMask | KeyReleaseMask
	    );

	    IMClose(g_in_stream);
	    g_open_window = 0;
	}
	break;

    case PointerMotionMask:
	{
	    int iz, iw, it;
	    float coord[3];

	    IWRtDisSec_ZWT(g_in_stream, &iz, &iw, &it);
	    IWXEvtToData(g_in_stream, iw, event, coord, &it);
	    update_display(
		g_in_stream,
		&(g_scales[g_open_window - 1].bounds[iw]),
		coord,
		iw,
		it,
		g_fast_mode,
		&g_need_clear
	    );
	}
	break;

    case KeyReleaseMask:
	{
	    int change_plot = 0;
	    int iz, iw, it;
	    float range, change;
	    float coord[3];
	    IntensityBounds* p_bounds;

	    IWRtDisSec_ZWT(g_in_stream, &iz, &iw, &it);
	    p_bounds = &(g_scales[g_open_window - 1].bounds[iw]);

	    switch (IWXEvtToKeySym(event)) {
	    case XK_Page_Up:  /* Fall through. */
	    case XK_0:
		range = p_bounds->max - p_bounds->min;
		if (range != 0.0f) {
		    change = 0.05f * range;
		    p_bounds->min += change;
		    p_bounds->max -= change;
		    if (p_bounds->max < p_bounds->min) {
			p_bounds->max = p_bounds->min;
		    }
		    change_plot = 1;
		}
		break;

	    case XK_Page_Down:
	    case XK_9:
		range = p_bounds->max - p_bounds->min;
		if (range == 0.0f) {
		    change = (fabs(p_bounds->min) < 1e-6) ?
			.05 : .05 * fabs(p_bounds->min);
		} else {
		    change = 0.05f * range;
		}
		p_bounds->min -= change;
		p_bounds->max += change;
		change_plot = 1;
		break;
	    }

	    if (change_plot) {
		IWXEvtToData(g_in_stream, iw, event, coord, &it);
		update_display(
		    g_in_stream,
		    p_bounds,
		    coord,
		    iw,
		    it,
		    g_fast_mode,
		    &g_need_clear
		);
	    }
	}
	break;
    }

    return 1;
}


static
int
quit(int* unused)
{
    WMCancelEventHandler(IW_ALL_WINDOWS, EnterWindowMask | LeaveWindowMask);
    WMCancelDisplayChange(IW_ALL_WINDOWS);
    (void) handle_win_event(LeaveWindowMask, NULL);
    if (g_scales != 0) {
	free(g_scales);
    }
    exit(0);
    return 1;
}


static
int
quit_wrapper(void)
{
    return quit(NULL);
}


/*
 * It would be more efficient to generate the horizontal and vertical graphs
 * directly in compute_graph and compute_complex_graph, but by
 * handling one case directly and using this routine to swap x and y for the
 * other, less code has to be written.
 */
static
void
swap_xy(int npts, IW_POINT* graph)
{
    int ipt;

    for (ipt = 0; ipt < npts; ++ipt) {
	float tmp = graph[ipt].x;

	graph[ipt].x = graph[ipt].y;
	graph[ipt].y = tmp;
    }
}


static
void
update_display(
    int in_stream,
    const IntensityBounds* p_bounds,
    const float coord[3],
    int iw,
    int it,
    FastGraphics fast_mode,
    int* p_need_clear
)
{
    int ix, iy, iz;
    int nw, nt, intlv;
    int size[3], subsize[3], start[3];

    ix = floor(coord[0]);
    iy = floor(coord[1]);
    iz = floor(coord[2]);
    IMRtSiz(in_stream, size, subsize, start);
    IMRtZWT(in_stream, &size[2], &nw, &nt, &intlv);
    if (ix >= 0 || ix < size[0] ||
	iy >= 0 || iy < size[1] ||
	iz >= 0 || iz < size[2]) {
	int pixel_fmt, is_complex;
	float* raw;
	IW_POINT* graph;

	IMRtMode(in_stream, &pixel_fmt);
	is_complex = pixel_fmt == IW_COMPLEX || pixel_fmt == IW_COMPLEX_SHORT;
	raw = (float*) malloc(
	    ((size[0] > size[1]) ? size[0] : size[1]) * sizeof(float) *
	    ((is_complex) ? 2 : 1)
	);
	graph = (IW_POINT*) malloc((size[0] + size[1]) * sizeof(IW_POINT));
	if (raw != 0 && graph != 0) {
	    int color = (IWRtColorMode(in_stream) == IW_PSEUDO) ?
		IWRtPsdGrColor(in_stream) : IWRtWaveGrColor(in_stream, iw);
	    int lower_left[2], upper_right[2];
	    int ids[4];
	    float smin, smax;

	    find_visible_region(
		in_stream, size[0], size[1], iw, lower_left, upper_right
	    );

	    
	    /* Compute horizontal profile. */
	    IMPosnZWT(in_stream, iz, iw, it);
	    IMRdPas(in_stream, raw, size[0], 1, 0, size[0] - 1, iy, iy);

	    /*
	     * In the graph a point whose intensity is p_bounds->min will
	     * have a y coordinate of smin and a point whose intensity is
	     * p_bounds->max will have a y coordinate of smax.  smin and
	     * smax are chosen to cover half of the visible image
	     * and to be centered on pixels (the extra 0.5f term does that).
	     */
	    smin = lower_left[1] + 0.5f +
		0.25f * (upper_right[1] - lower_left[1]);
	    smax = upper_right[1] + 0.5f -
		0.25f * (upper_right[1] - lower_left[1]);

	    if (is_complex) {
		compute_complex_graph(
		    size[0],
		    raw,
		    p_bounds->min,
		    p_bounds->max,
		    smin,
		    smax,
		    iz,
		    IWRtComplexDis(in_stream),
		    graph
		);
	    } else {
		compute_graph(
		    size[0],
		    raw,
		    p_bounds->min,
		    p_bounds->max,
		    smin,
		    smax,
		    iz,
		    graph
		);
	    }

	    /* Compute vertical profile */
	    IMPosnZWT(in_stream, iz, iw, it);
	    IMRdPas(in_stream, raw, 1, size[1], ix, ix, 0, size[1] - 1);

	    smin = lower_left[0] + 0.5f +
		0.25f * (upper_right[0] - lower_left[0]);
	    smax = upper_right[0] + 0.5f -
		0.25f * (upper_right[0] - lower_left[0]);

	    if (is_complex) {
		compute_complex_graph(
		    size[1],
		    raw,
		    p_bounds->min,
		    p_bounds->max,
		    smin,
		    smax,
		    iz,
		    IWRtComplexDis(in_stream),
		    graph + size[0]
		);
	    } else {
		compute_graph(
		    size[1],
		    raw,
		    p_bounds->min,
		    p_bounds->max,
		    smin,
		    smax,
		    iz,
		    graph + size[0]
		);
	    }
	    swap_xy(size[1], graph + size[0]);

	    if (fast_mode == FAST_OVERLAY) {
		if (color >= IWRtNumOVColors()) {
		    color = 1;
		}
		IWClearOverlay(in_stream);
		IWStartOverlay(in_stream, color);
	    } else {
		IWDisFastGr(in_stream);
	    }

	    ids[0] = IWGrAddLns5D(
		in_stream, graph, size[0], iw, it, 1, 1, color
	    );
	    ids[1] = IWGrAddLns5D(
		in_stream, graph + size[0], size[1], iw, it, 1, 1, color
	    );
	    /*
             * Draw horizontal and vertical line segments through the point.
	     * The segments span the image (the .001 term was chosen
	     * such that for all possible zooms (up to 32); no gap is visible
	     * at the end of the line).
	     */
	    ids[2] = IWGrAddLn5D(
		in_stream,
		0.0f, coord[1], iz,
		size[0] - 0.001f, coord[1], iz,
		iw,
		it,
		1,
		1,
		color
	    );
	    ids[3] = IWGrAddLn5D(
		in_stream,
		coord[0], 0.0f, iz,
		coord[0], size[1] - 0.001f, iz,
		iw,
		it,
		1,
		1,
		color
	    );
	    IWGrDisGrList(in_stream, ids, 4);

	    if (fast_mode == FAST_OVERLAY) {
		IWEndOverlay(in_stream);
	    }

	    /*
	     * Something else could cause a redraw.  In that case if overlay
	     * planes are used, the plot will be redrawn in the normal planes
	     * where this would not clear them (it is only clearing the
	     * overlay).  If the fast graphics calls via the monitor are used,
	     * this is not an issue, but if the redraw added additonal
	     * graphics and this application does not receive a display
	     * change event (that will cause it to rerecord the background)
	     * then those graphics will be hidden when the recorded background
	     * is displayed with IWDisFastGr.
	     *
	     * In any case remove the graphics, it simplifies bookkeeping for
	     * the monitor fast graphics case and for the overlay case
	     * reduces the chance of clearing artifacts.  If something else
	     * causes a clear / redraw, the plot will not be redrawn unless
	     * the mouse moves; this does not seem to cause problems in
	     * practice.
	     */
	    IWGrRmGrList(in_stream, ids, 4);
	    *p_need_clear = 1;
	} else {
	    if (*p_need_clear) {
		if (fast_mode == FAST_OVERLAY) {
		    IWClearOverlay(in_stream);
		} else {
		    IWDisFastGr(in_stream);
		}
		*p_need_clear = 0;
	    }
	}
	if (graph != 0) {
	    free(graph);
	}
	if (raw != 0) {
	    free(raw);
	}
    } else {
	if (*p_need_clear) {
	    if (fast_mode == FAST_OVERLAY) {
		IWClearOverlay(in_stream);
	    } else {
		IWDisFastGr(in_stream);
	    }
	    *p_need_clear = 0;
	}
    }
}
