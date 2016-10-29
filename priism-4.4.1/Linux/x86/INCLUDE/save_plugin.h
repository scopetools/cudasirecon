#ifndef UCSF_MSG_SAVE_PLUGIN_H
#define UCSF_MSG_SAVE_PLUGIN_H

#include "pixel_format.h"
#include "WMInclude.h"
#include "ive_plugin.h"


typedef IVEPlugin SavePlugin;

typedef struct AvailableSavePlugin {
    /* User interface name for format */
    char* name;
    /* Short (command-line) name for format */
    char* short_name;
    /* Path to plugin */
    char* path;
    /* Plugin interface to load */
    char* interface_name;
} AvailableSavePlugin;

struct SavedSequenceImpl;
typedef struct SavedSequenceImpl* SavedSequence;

typedef struct SavePluginAttrExt {
    struct SavePluginAttrExt* next;
    void* data;
    int type_code;
} SavePluginAttrExt;

typedef struct SavePluginAttr {
    /* Function pointers set by the plugin for use by the client. */

    /*
     * Creates the plugin-specific structure to hold user-specifiable options.
     * Returns non-zero if the operation was successful and zero if the
     * operation failed.
     *
     * The function pointer may be NULL to indicate that the plugin has no
     * user-specifiable options.
     */
    void* (*options_creator)(struct SavePluginAttr* attr);

    /*
     * Frees the resources associated with a value returned from
     * options_creator.
     *
     * The function pointer may be NULL.
     */
    void (*options_destroyer)(struct SavePluginAttr* attr, void* options);

    /*
     * Intended for use with the ProcFunc interface command-line handling,
     * parses a single argument and updates the current options.  The
     * return value must be one of the following:
     *     0:  arg was successfully parsed and options was updated.
     *     1:  arg is relevant to this plugin but was malformed in some way;
     *         the registered error handler, if any, is called to report
     *         the reason.
     * The function pointer may be NULL to inhibit setting the options via
     * the command line.
     */
    int (*single_arg_parser)(
	struct SavePluginAttr* attr, const char* arg, void* options
    );

    /*
     * Prints the command-line arguments corresponding to the passed
     * plugin options.  The third argument is printed as a prefix to the
     * command line options.  The fourth argument is the function that will be
     * called to print each argument.  Its first argument is the fifth
     * argument to *args_printer, and its second argument is a printf-style
     * format string.
     *
     * The args_printer function pointer may be NULL.
     */
    void (*args_printer)(
	struct SavePluginAttr* attr,
	const void* options,
	const char* prefix,
	void (*printer)(void* completion, const char* fmt, ...),
	void* completion
    );

    /*
     * Prints a description of the command line arguments for the plugin.
     * The second argument is the function that will be called to print each
     * line.  Its first argument is the third argument to *usage_printer,
     * and its second argument is a printf-style format string.
     *
     * The usage_printer function pointer may be NULL.
     */
    void (*usage_printer)(
	struct SavePluginAttr* attr,
	void (*printer)(void* completion, const char* fmt, ...),
	void* completion
    );

    /*
     * Creates the WM based dialog for editing the user-specifiable options
     * and returns the id for the menu or 0 if the call failed.  The menu
     * created is bound to the set of options passed and the lifetime of that
     * set of options must exceed the lifetime of the menu.  The parent
     * is used as a hint for initial placement of the dialog; it may be NULL.
     *
     * The function pointer may be NULL to indicate that the user should
     * not be able to adjust the options interactively.
     */
    int (*options_menu_creator)(
	const char* title_prefix,
	Widget parent,
	struct SavePluginAttr* attr,
	void* options
    );

    /*
     * Coerces the current options and performs updates to the WM user
     * interface (if any) to bring the options in agreement with the
     * specified format.
     */
    void (*options_updater)(
	struct SavePluginAttr* attr, PixelFormat fmt, void* options
    );

    /*
     * Returns zero if the pixel format is acceptable and one if the pixel
     * format is not acceptable.
     *
     * This function pointer must not be NULL.
     */
    int (*pixel_format_verifier)(struct SavePluginAttr* attr, PixelFormat fmt);

    /*
     * Returns the allowed pixel format that is the nearest match to the
     * given format, fmt, or null if no match could be found.
     *
     * This function pointer must not be NULL.
     */
    PixelFormat (*pixel_format_matcher)(
	struct SavePluginAttr* attr, PixelFormat fmt
    );

    /*
     * Returns a nonzero value if the plugin will generate one file per image
     * with the given options or zero if the plugin will generate one file for
     * the entire sequence with the given options.
     *
     * This function pointer must not be NULL.
     */
    int (*image_per_file_getter)(
	struct SavePluginAttr* attr, const void* options
    );

    /*
     * Returns the suggested extension for file names.
     *
     * This function pointer may be null.
     */
    const char* (*extension_getter)(
	struct SavePluginAttr* attr, const void* options, PixelFormat fmt
    );

    /*
     * Initiates the process of saving a sequence of images.  If saving
     * the sequence to a single file, name will be the name of that file;
     * otherwise, name is null.  nx is the width, in pixels, of the images.
     * ny is the height, in pixels of the images.  image_count is the number
     * of images or -1 if the total number of images is not known yet.  fmt
     * is the layout of each pixel.  Returns a handle to the sequence or
     * NULL if the operation failed.  In case of failure, the error handler,
     * if any, is called to report the reason for the failure.
     *
     * This function pointer must not be NULL.
     */
    SavedSequence (*sequence_starter)(
	struct SavePluginAttr* attr,
	const void* options,
	const char* name,
	int nx,
	int ny,
	int image_count,
	PixelFormat fmt
    );

    /*
     * Adds an image to the end of the sequence.  If saving the sequence to a
     * single file, name is null; otherwise it is the name of the file to use
     * for this image.  Assumes that the (0,0) element of the image is at
     * the lower left corner when displayed.  Returns zero if successful and
     * one if not.  In case of failure, the error handler, if any, is called
     * to report the reason for the failure.
     *
     * This function pointer must not be NULL.
     */
    int (*image_adder)(
	SavedSequence sequence, const char* name, const void* image
    );

    /*
     * Completes the process of saving a sequence of images.  Returns zero
     * if successful and a non-zero value if not.  In case of failure,
     * the error handler, if any, is called to report the reason for the
     * failure.
     *
     * This function pointer must not be NULL.
     */
    int (*sequence_ender)(SavedSequence sequence, int actual_count);

    /* NULL-terminated list of extension structures set by the plugin. */
    SavePluginAttrExt* extensions;

    /* Function pointers and values set by the client for use by the plugin. */

    /*
     * Called by the plugin from within *single_arg_parser, *sequence_starter,
     * *image_adder, and *sequence_ender to report the reason for a failure.
     * The first argument passed is the value of error_handler_arg; the second
     * is a printf-style format string.
     */
    void (*error_handler)(void* arg, const char* fmt, ...);
    void* error_handler_arg;

    /*
     * Called by the plugin from within *single_arg_parser, *options_updater,
     * or the handlers for the options menu to indicate change that causes
     * *singlefile_getter to report a different value than it did.  The first
     * argument passed is the value of single_change_handler_arg.
     */
    void (*single_change_handler)(void* arg);
    void* single_change_handler_arg;

    /* Reserved for internal use by the plugin. */
    void* p_internal;
} SavePluginAttr;

typedef struct SavePluginIndirectCallback {
    void (*func)(void*);
    void* arg;
} SavePluginIndirectCallback;

/* Recognized extensions */
enum {
    SAVE_PLUGIN_EXT_NONE = -1
};


#ifdef __cplusplus
extern "C" {
#endif

/*
 * Interface for users of plugins.
 */

/*
 * Reads plugins from the file named config_file and returns them in *p_list
 * with the number of plugins returned in *p_count.  If config_file is not a
 * fully specified path (does not start with / or ./) then for each colon
 * separated entry in search_path, prepends the entry (after environment
 * variable expansion) to config_file and appends the plugins from that file
 * to the list of plugins returned to the caller.  If search_path is NULL,
 * uses "$HOME/.iveprefs:$IVEBASE/CONFIG".
 */
void get_save_plugin_list(
    const char* config_file,
    const char* search_path,
    AvailableSavePlugin** p_list,
    int* p_count
);

void destroy_save_plugin_list(AvailableSavePlugin* list);

SavePlugin load_save_plugin(AvailableSavePlugin* p_available);

void unload_save_plugin(SavePlugin plugin);

SavePluginAttr* get_save_plugin_attr(SavePlugin plugin);

/*
 * Returns type_code if the plugin supports the extension with of the type
 * type_code; otherwise returns SAVE_PLUGIN_EXT_NONE.  If p_extension_data
 * is not NULL and the extension is supported, sets *p_extension_data to the
 * address of the extension's data.
 */
int search_save_plugin_extensions(
    SavePlugin plugin, int type_code, void** p_extension_data
);

/*
 * Interface for plugins.
 */

/*
 * Because the plugins are loaded and unloaded, the addresses of callback
 * functions for user-interface controls can become invalid if those functions
 * are in the plugin.  So provide useful callback functions in the
 * libSavePlugin library which is not loaded and unloaded.
 */

int save_plugin_hide_menu(int* p_menu_id);

/* Expects arg to be a SavePluginAttr** cast to an int*. */
int save_plugin_propagate_singlemultiple_change(int* arg);

/*
 * Provides a generic indirection mechanism; expects arg to be a
 * SavePluginIndirectCallback* cast to an int*.
 */
int save_plugin_redirect_cb(int* arg);

/*
 * A convenience function to call after WMInitSubMenu in order to position
 * the dialog.  parent may be NULL.
 */
void save_plugin_position_menu(Widget parent);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
