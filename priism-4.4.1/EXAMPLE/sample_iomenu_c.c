/*
 * Demonstrates the C/C++ interface to the iomenu utility.
 *
 * Assume that you had a command-line application which expected the following
 * command-line:
 *
 *     ffilter in out optional_otf -x=start:end -y=start:end -z=start:end \
 *         -t=start:end:step -w=wave1:wave2:... -bandpass=l1:l2 -enhance=e \
 *         -smooth=s -fraction=f -apodize=n -nzpad=n -2D -norm -orig -mode2
 *
 * (this if the command-line used by the ffilter application in Priism), and
 * want to write a graphical user interface as a wrapper.  The iomenu utility
 * will automatically handle the region selection arguments (-x, -y, -z, -t,
 * -w); the code below will handle the others.  If you want to see what the
 * corresponding graphical user interface wrapper in Priism looks like, go
 * to the Processing->Fourier Tools->FFilter entry in the Priism menu bar.
 */


#include "iomenu.h"     /* To pick up the declarations of the structures and
		         * iomenu()
			 */
#include "WMInclude.h"  /* To pick up the WM function declarations since the
			 * WM library will be the mechanism for creating the
			 * user interface components.
			 */
#include <limits.h>      /* INT_MAX */
#include <string.h>      /* strcmp, strcpy */


/*
 * An instance of this structure will carry around the values of the
 * parameters associated with the command line switches and any state I need
 * to keep track of for the user interface.
 */
typedef struct FFilterInternals {
    /* Stores the two -bandpass= parameter values. */
    float bandpass[2];
    /* Stores the -enhance= parameter value.  Range limited to 0 to 1. */
    float enhancement;
    /* Stores the -smooth= parameter value.  Range limited to 0 to 1. */
    float smoothing;
    /* Stores the -fraction= parameter value.  Range limited to 0 to 1. */
    float fraction;
    /* Stores the -apodize= parameter value.  Must not be negative. */
    int border;
    /*
     * Stores the -nzpad= parameter value.  Must be greater than or equal to
     * the number of z sections to process.
     */
    int padded_size;
    /* Is nonzero if the -2D option is to be printed. */
    int force2d;
    /* Is nonzero if the -norm option is to be printed. */
    int normalize;
    /* Is nonzero if the -orig option is to be printed. */
    int preserve_mean;
    /* Is nonzero if the -mode2 option is to be printed. */
    int force_float;
    /*
     * Store the minimum allowable padded size.
     */
    int min_padded;
    /*
     * Storage for the user interface component IDs so they can be hidden,
     * displayed, or updated depending on actions by the user.
     */
    Widget bandpass_help, bandpass_field;
    Widget enhance_help, enhance_field;
    Widget smooth_help, smooth_field;
    Widget pad_help, pad_field;
} FFilterInternals;


/* Declare prototypes for callback routines. */
static void customize_menu(Iomenu* iom);
static void customize_special_menu(Iomenu* iom);
static int handle_bandpass(int* arg);
static void handle_file(Iomenu* iom);
static int handle_force2d(int* arg);


int
main(int argc, char* argv[])
{
    FFilterInternals internals = {
	{ 0.0f, 0.0f }, /* By default, the bandpass filter is disabled. */
	0.90f,          /* default enhancement value */
	0.15f,          /* default smoothing value */
	0.00f,          /* default fraction value */
	0,              /* By default, do not apodize edges. */
	0,              /* Do not know z size yet, so initialize padded size
			   to zero. */
	0,              /* By default, do not force 2D processing. */
	1,              /* By default, normalize the CTF. */
	1,              /* By default, preserve the mean value. */
	0,              /* By default, write the output with the same pixel
			   type as the input. */
	0,              /* Do not know z size yet, set minimum padded size to
			   zero. */
	0,              /* Will initialize the Widget values when created. */
	0,
	0,
	0,
	0,
	0,
	0,
	0
    };
    /*
     * The first entry on the command line is an input image file.  I'll use
     * that file to set the default region to process and the default name of
     * the output file so I mark the file as an IOMENU_PRIMARY_IMAGE.  I
     * initialize the file name for the first file to be empty (do not not know
     * value yet), the file filter to be "*" (any file will do), and the button
     * label in the menu to be "InFile".  The second entry on the command line
     * is the output file.  I'll generate the default output name based on the
     * name of the first file so I mark the file as IOMENU_OUTPUT_DEFAULT.  I
     * initialize the initial file name to be empty (when the user selects the
     * input file in the interface the default name will be filled in anyways)
     * and set the file filter to be "_flt" so that the default output name
     * is the input name + "_flt".  The label that will appear in the user
     * interface is "OutFile".  The third entry on the command line is an
     * optional input file which will have the special value, "none", if it
     * is not specified.  If I could change the command-line application, this
     * optional file would be better handled via a command-line switch, but
     * since I can not, I'll declare it as an output (that way iomenu does not
     * check that the file exists which it would if I used IOMENU_OTHER_INPUT).
     * I'll initially set the third file to the special value "none".  The
     * file filter is set to "*" (though the value will not be used),
     * and the label in the user interface will be "OTFfile".
     */
    IomenuFile files[3] = {
	{ "", "*", "InFile", IOMENU_PRIMARY_IMAGE },
	{ "", "_flt", "OutFile", IOMENU_OUTPUT_DEFAULT },
	{ "none", "*", "OTFfile", IOMENU_OTHER_OUTPUT }
    };
    /*
     * The remaining command line entries are ten switches.  The parameters[]
     * array will define them, but because there are unions involved, I'll
     * initialize this below rather than directly in the declaration.
     */
    IomenuParameter parameters[10];
    /* Passed to iomenu.  Filled in below. */
    Iomenu menu_description;
    /*
     * Used as the label for the button to open the special parameters dialog.
     * Because it is used as the label for a button added by WMAddFuncButton
     * and the iomenu utility adds Save and Restore buttons with
     * WMAddSaveButton, safest to use a character array with some extra space
     * rather than a literal string.
     */
    char special_label[32] = "Special Parameters";

    (void) strcpy(parameters[0].label, "-bandpass=");
    parameters[0].v.r.values = internals.bandpass;
    parameters[0].v.r.count = 2;
    parameters[0].type = IOMENU_REAL;

    (void) strcpy(parameters[1].label, "-enhance=");
    parameters[1].v.r.values = &internals.enhancement;
    parameters[1].v.r.count = 1;
    parameters[1].type = IOMENU_REAL;

    (void) strcpy(parameters[2].label, "-smooth=");
    parameters[2].v.r.values = &internals.smoothing;
    parameters[2].v.r.count = 1;
    parameters[2].type = IOMENU_REAL;

    (void) strcpy(parameters[3].label, "-fraction=");
    parameters[3].v.r.values = &internals.fraction;
    parameters[3].v.r.count = 1;
    parameters[3].type = IOMENU_REAL;

    (void) strcpy(parameters[4].label, "-apodize=");
    parameters[4].v.i.values = &internals.border;
    parameters[4].v.i.count = 1;
    parameters[4].type = IOMENU_INTEGER;

    (void) strcpy(parameters[5].label, "-nzpad=");
    parameters[5].v.i.values = &internals.padded_size;
    parameters[5].v.i.count = 1;
    parameters[5].type = IOMENU_INTEGER;

    (void) strcpy(parameters[6].label, "-2D");
    parameters[6].v.pb = &internals.force2d;
    parameters[6].type = IOMENU_ACTIVE_WHEN_NONZERO;

    (void) strcpy(parameters[7].label, "-norm");
    parameters[7].v.pb = &internals.normalize;
    parameters[7].type = IOMENU_ACTIVE_WHEN_NONZERO;

    (void) strcpy(parameters[8].label, "-orig");
    parameters[8].v.pb = &internals.preserve_mean;
    parameters[8].type = IOMENU_ACTIVE_WHEN_NONZERO;

    (void) strcpy(parameters[9].label, "-mode2");
    parameters[9].v.pb = &internals.force_float;
    parameters[9].type = IOMENU_ACTIVE_WHEN_NONZERO;

    menu_description.files = files;
    menu_description.parameters = parameters;
    menu_description.closure = &internals;
    menu_description.file_count = sizeof(files) / sizeof(files[0]);
    menu_description.parameter_count =
	sizeof(parameters) / sizeof(parameters[0]);

    iomenu(
	"2D & 3D Fourier Filter Setup Menu",    /* dialog title */
	special_label,                          /* name for special parameters
						   button */
	"ffilter",                              /* name of command-line
						   application */
	&menu_description,                      /* pointer to the Iomenu
						   structure */
	customize_menu,                         /* routine to build main
						   dialog */
	customize_special_menu,                 /* routine to build special
						   parameters dialog */
	handle_file,                            /* routine to handle file
						   changes */
	0                                       /* this application does not
						   require special handling
						   at exit */
    );

    return 1;
}


static
void
customize_menu(Iomenu* iom)
{
    /* Use static storage for these fixed parameter bounds. */
    static float fzero = 0.0f;
    static float fone = 1.0f;
    static float f1k = 1000.0f;
    static int i1 = 1;
    static int i2 = 2;
    FFilterInternals* p_internals = (FFilterInternals*) iom->closure;

    WMNewRow();
    p_internals->bandpass_help =
	WMAddInfoButton("Bandpass (res min,max)", "ffilter bandpass");
    p_internals->bandpass_field = WMAddFloatField(
	p_internals->bandpass,
	2,
	10,
	20,
	&fzero,
	&f1k,
	&i1,
	handle_bandpass,
	p_internals,
	0,
	0
    );

    WMNewRow();
    p_internals->enhance_help =
	WMAddInfoButton("Enhancement (0-1)", "ffilter enhancement");
    p_internals->enhance_field = WMAddFloatField(
	&p_internals->enhancement, 1, 10, 10, &fzero, &fone, &i2, 0, 0, 0, 0
    );

    WMNewRow();
    p_internals->smooth_help =
	WMAddInfoButton("Smoothing   (0-1)", "ffilter smoothing");
    p_internals->smooth_field = WMAddFloatField(
	&p_internals->smoothing, 1, 10, 10, &fzero, &fone, &i2, 0, 0, 0, 0
    );

    WMNewRow();
    WMAddInfoButton("Fraction original image (0-1)", "ffilter fract");
    WMAddFloatField(
	&p_internals->fraction, 1, 8, 8, &fzero, &fone, &i2, 0, 0, 0, 0
    );

    WMNewRow();
    WMAddInfoButton("Force Processing as 2D data", "ffilter twod");
    WMAddToggleButton(
	"", &p_internals->force2d, handle_force2d, p_internals, 0, 0
    );

    /*
     * Force the special parameters button to be offset from the parameter
     * controls and the left side of the dialog.
     */
    WMSetOffset(100, 0, 20, 0);
}


static
void
customize_special_menu(Iomenu* iom)
{
    /* Use static storage for fixed parameter bounds. */
    static int i0 = 0;
    static int i1 = 1;
    static int i2 = 2;
    static int imax = INT_MAX;
    FFilterInternals* p_internals = (FFilterInternals*) iom->closure;

    WMNewRow();
    WMAddInfoButton("Border rolloff, # pixels", "ffilter apodize");
    WMAddIntField(
	&p_internals->border, 1, 12, 12, &i0, &imax, &i1, 0, 0, 0, 0
    );

    WMNewRow();
    p_internals->pad_help =
	WMAddInfoButton("Size for Z transforms", "ffilter nzpad");
    p_internals->pad_field = WMAddIntField(
	&p_internals->padded_size,
	1,
	12,
	12,
	&p_internals->min_padded,
	&imax,
	&i2,
	0,
	0,
	0,
	0
    );

    WMNewRow();
    WMAddInfoButton("Normalize CTF", "ffilter norm");
    WMAddToggleButton("", &p_internals->normalize, 0, 0, 0, 0);
    WMSetOffset(20, 0, 0, 0);
    WMAddInfoButton("Preserve mean", "ffilter orig");
    WMSetOffset(0, 0, 0, 0);
    WMAddToggleButton("", &p_internals->preserve_mean, 0, 0, 0, 0);

    WMNewRow();
    WMAddInfoButton("Force floating-point", "ffilter mode2");
    WMAddToggleButton("", &p_internals->force_float, 0, 0, 0, 0);

    /*
     * Force the OK button to be offset from the parameter controls and the
     * left side of the dialog.
     */
    WMSetOffset(70, 0, 20, 0);
}


/*
 * If the bandpass filter is enabled, hide the enhancement and smoothing
 * controls (since they are not relevant).  Otherwise, display the enhancement
 * and smoothing controls.
 *
 * It is probably better to grey out the unused controls (i.e. WMDisableField
 * and WMEnableField rather than WMUnpasteWidget and WMPasteWidget) because
 * pasting and unpasting can cause geometry problems if the dialog has not
 * been realized.  Hide and display the controls because that is what ffilter_i
 * has traditionally done.
 * done.
 */
static
int
handle_bandpass(int* arg)
{
    FFilterInternals* p_internals = (FFilterInternals*) arg;

    if (p_internals->bandpass[0] + p_internals->bandpass[1] >= 1e-4f ) {
	WMUnpasteWidget(p_internals->enhance_help);
	WMUnpasteWidget(p_internals->enhance_field);
	WMUnpasteWidget(p_internals->smooth_help);
	WMUnpasteWidget(p_internals->smooth_field);
    } else {
	WMPasteWidget(p_internals->enhance_help);
	WMPasteWidget(p_internals->enhance_field);
	WMPasteWidget(p_internals->smooth_help);
	WMPasteWidget(p_internals->smooth_field);
    }

    return 0;
}


static
void
handle_file(Iomenu* iom)
{
    FFilterInternals* p_internals = (FFilterInternals*) iom->closure;

    /*
     * If the first file (main input) changed or the z region changed, update
     * the padded size.
     */
    if (iom->iparm == 1 || iom->iparm == -3) {
	int nz = iom->ixyztw[2][1] - iom->ixyztw[2][0] + 1;

	p_internals->min_padded = nz;
	p_internals->padded_size = 2;
	while (p_internals->padded_size < 1.3 * nz) {
	    p_internals->padded_size *= 2;
	}
	if (p_internals->pad_field != 0) {
	    WMUpdateField(p_internals->pad_field);
	}
    } else if (iom->iparm == 3) {
	/*
	 * If the third file (optional ctf) changed display the bandpass
	 * widgets if a CTF has not been supplied or hide them if a CTF
	 * has been supplied (the bandpass filter is not used in that case).
	 *
	 * It is probably better to grey out the unused controls (i.e.
	 * WMDisableField and WMEnableField rather than WMUnpasteWidget
	 * and WMPasteWidget) because pasting and unpasting can cause geometry
	 * problems if the dialog has not been realized.  Hide and display
	 * the controls because that is what ffilter_i has traditionally done.
	 */
	if (strcmp(iom->files[2].filename, "none") == 0 ||
	    iom->files[2].filename[0] == '\0') {
	    if (p_internals->bandpass_help != 0) {
		WMPasteWidget(p_internals->bandpass_help);
		WMPasteWidget(p_internals->bandpass_field);
	    }
	} else {
	    if (p_internals->bandpass_help != 0) {
		WMUnpasteWidget(p_internals->bandpass_help);
		WMUnpasteWidget(p_internals->bandpass_field);
	    }
	}
    }
}


/*
 * If forcing 2D processing, hide the padding controls (they are not relevant);
 * otherwise, display the padding controls.
 *
 * It is probably better to grey out the unused controls (i.e. WMDisableField
 * and WMEnableField rather than WMUnpasteWidget and WMPasteWidget) because
 * pasting and unpasting can cause geometry problems if the dialog has not
 * been realized.  Hide and display the controls because that is what ffilter_i
 * has traditionally done.
 * done.
 */
static
int
handle_force2d(int* arg)
{
    FFilterInternals* p_internals = (FFilterInternals*) arg;

    if (p_internals->force2d) {
	if (p_internals->pad_help != 0) {
	    WMUnpasteWidget(p_internals->pad_help);
	    WMUnpasteWidget(p_internals->pad_field);
	}
    } else {
	if (p_internals->pad_help != 0) {
	    WMPasteWidget(p_internals->pad_help);
	    WMPasteWidget(p_internals->pad_field);
	}
    }

    return 0;
}
