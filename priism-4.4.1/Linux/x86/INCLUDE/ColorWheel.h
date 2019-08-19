#ifndef UCSF_MSG_COLORWHEEL_H
#define UCSF_MSG_COLORWHEEL_H

#include "WMInclude.h"

typedef struct ColorWheelImpl* ColorWheel;
typedef struct ColorBarsImpl* ColorBars;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Given a widget, w, which was created by WMAddGLwMDrawWidget, sets up for
 * drawing a color wheel with a bar drawn below it in the currently selected
 * color.  The colors in the color wheel are drawn with the given value for
 * the saturation (saturation must be between 0 (everything black) and 1).
 * cur_sel_rgb holds the red, green, and blue components for the initially
 * selected current color.  If bar_as_gradient is nonzero, the bar drawn
 * in the currently selected color is filled with a color gradient from black
 * to the value of the currently selected color; otherwise the bar is a solid
 * block of the currently selected color.  If delayed update is nonzero,
 * function calls which modify the color wheel's attributes do not
 * automatically trigger a redraw; in that case you can force a redraw with
 * the redraw_color_wheel() function.  Otherwise, the function calls which
 * modify the color wheel's attributes do trigger a redraw.  Returns a non-NULL
 * value if the operation succeeded.  Returns NULL if the operation failed.
 */
ColorWheel create_color_wheel(
    Widget w,
    float saturation,
    const float cur_sel_rgb[3],
    int bar_as_gradient,
    int delayed_update
);

void destroy_color_wheel(ColorWheel cw);

void redraw_color_wheel(ColorWheel cw);
/*
 * Adjust the color wheel size for a change in the size of the drawing canvas.
 * Returns a nonzero value if could not adapt to the new size.
 */
int resize_color_wheel(ColorWheel cw);

/*
 * Returns the saturation (in [0, 1]) level currently used to draw the color
 * wheel
 */
float get_color_wheel_saturation(const ColorWheel cw);
void get_color_wheel_selection_rgb(const ColorWheel cw, float rgb[3]);
/* Returns the saturation level of the currently selected color. */ 
float get_color_wheel_selection_saturation(const ColorWheel cw);
int get_color_wheel_bar_as_gradient(const ColorWheel cw);
int get_color_wheel_delayed_update(const ColorWheel cw);

/*
 * Changes the saturation level used to draw the color wheel.  Leaves
 * the currently selected color unmodified.
 */
void set_color_wheel_saturation(ColorWheel cw, float saturation);
/* Sets the selected color from the X event coordinates. */
void set_color_wheel_selection_from_position(ColorWheel cw,  int x, int y);
void set_color_wheel_selection_rgb(ColorWheel cw, const float rgb[3]);
/* Adjusts the saturation level of the currently selected color. */
void set_color_wheel_selection_saturation(ColorWheel cw, float saturation);
void set_color_wheel_bar_as_gradient(ColorWheel cw, int bar_as_gradient);
void set_color_wheel_delayed_update(ColorWheel cw, int delayed_update);


/*
 * Given a widget, w, which was created by WMAddGLwMDrawWidget, sets up for
 * drawing a horizontal bar, broken into ncolor blocks.  Returns a non-NULL
 * value if the operation succeeded.  Returns NULL if the operation failed.
 */
ColorBars create_color_bars(Widget w, int ncolor);

void destroy_color_bars(ColorBars cb);

void redraw_color_bars(ColorBars cb);
/*
 * Adjust the color bars for a change in the size of the drawing canvas.
 * Returns a nonzero value if could not adapt to the new size.
 */
int resize_color_bars(ColorBars cb);

/*
 * Fills rgba with the red, green, blue, and alpha, components, respectively of
 * the indexth color in the color bar.  The components are in the range of zero
 * to one.  Returns zero if index is greater than or equal to zero and less
 * than or equal to the number of colors.  Otherwise, returns a non-zero value.
 */
int get_color_bars_component(const ColorBars cb, int index, float rgba[4]);

/* Returns the number of components in the color bar. */
int get_color_bars_count(const ColorBars cb);

/*
 * Sets the indexth component of the color bar to have the color (given by
 * red, green, blue, and alpha components each in the range of zero to one)
 * in rgba.  Returns zero if index is greater than or equal to zero and less
 * than or equal to the number of colors and the components of rgba are in
 * range.  Otherwise, returns a non-zero value.
 */
int set_color_bars_component(
    const ColorBars cb, int index, const float rgba[4]
);

/*
 * Modifies the color bars to have ncolor components.  Return zero if the
 * operation succeeded.  Returns a nonzero value if ncolor is less than one
 * or if the color bars could not be modified to accommodate ncolor components.
 */
int set_color_bars_count(const ColorBars cb, int ncolor);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
