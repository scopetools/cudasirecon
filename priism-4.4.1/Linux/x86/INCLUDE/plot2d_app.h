#ifndef UCSF_MSG_PLOT2D_APP_H
#define UCSF_MSG_PLOT2D_APP_H

#define PLOT2D_QUIT    1
#define UPDATE_FILE    2
#define UPDATE_SCALE   3
#define SET_GRAPH_NUM  4
#define SHOW_LEGEND    5
#define SET_CURVE_ATTR 6
#define SET_TICKS      7
#define SET_WINDOW_SIZE 8
#define CHANGE_LABEL   9
#define SET_TICKS_EXT  10

typedef struct plot2d_cmd_new_file {
    int length;
} plot2d_cmd_new_file;

typedef struct plot2d_cmd_set_attr {
    int all_graphs;
    int curve_num;
    int marker_type;
    int char_marker;
    int line;
    int color;
} plot2d_cmd_set_attr;

typedef struct plot2d_cmd_set_graph {
    int graph_number;
} plot2d_cmd_set_graph;

typedef struct plot2d_cmd_set_ticks {
    int num_ticks_x;
    int num_ticks_y;
    int tick_labels_x;
    int tick_labels_y;
} plot2d_cmd_set_ticks;

typedef struct plot2d_cmd_set_ticks_ext {
    int num_ticks_x;
    int num_ticks_y;
    int tick_labels_x;
    int tick_labels_y;
    int nplaces_x;
    int nplaces_y;
} plot2d_cmd_set_ticks_ext;

typedef struct plot2d_cmd_set_scale {
    float x_min;
    float x_max;
    float y_min;
    float y_max;
} plot2d_cmd_set_scale;

typedef struct plot2d_cmd_set_legend {
    int on;
} plot2d_cmd_set_legend;

typedef struct plot2d_cmd_set_label {
    int label_number;
    int length;
} plot2d_cmd_set_label;

typedef struct plot2d_cmd_set_window {
    int x_start;
    int y_start;
    int width;
    int height;
} plot2d_cmd_set_window;

typedef struct plot2d_cmd {
    int type;
    union {
	plot2d_cmd_new_file new_file;
	plot2d_cmd_set_attr set_attr;
	plot2d_cmd_set_graph set_graph;
	plot2d_cmd_set_ticks set_ticks;
	plot2d_cmd_set_scale set_scale;
	plot2d_cmd_set_legend set_legend;
	plot2d_cmd_set_label set_label;
	plot2d_cmd_set_window set_window;
	plot2d_cmd_set_ticks_ext set_ticks_ext;
    } var;
} plot2d_cmd;

#endif /* INCLUDE_GUARD */
