#ifndef UCSF_MSG_IVE_PLUGIN_H
#define UCSF_MSG_IVE_PLUGIN_H

/*
 * Copyright(C) 2002
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Provides a simple dynamic loading (plugin) interface.  The plugin has to
 * provide a function function of the form
 *     int your_name_here(
 *         const char* interface_name,
 *         void (*error_reporter)(const char* format, ...),
 *         void** p_attr,
 *         IVEPluginCleanupProc* p_cleanup_func
 *     );
 * which should return 0 if successful and nonzero if not.  If not successful
 * it can report the cause of the error with error_reporter.  If successful,
 * it can store the address of the plugin attributes in *p_attr and
 * the address of the routine to call when unloading the plugin in
 * *p_cleanup_func.  The argument passed to the cleanup routine is the address
 * of the plugin attributes.  The name of the function should be the same as
 * the name the plugin loader passes to IVELoadPlugin as the init_symbol_name
 * argument or initialize_plugin if the plugin loader passes NULL for the
 * init_symbol_name argument.
 */

struct IVEPluginImpl;
typedef struct IVEPluginImpl* IVEPlugin;

typedef void (*IVEPluginCleanupProc)(void* attr);
typedef int (*IVEPluginInitProc)(
    const char* interface_name,
    void (*error_reporter)(const char* format, ...),
    void** p_attr,
    IVEPluginCleanupProc* p_cleanup_func
);
    
#ifdef __cplusplus
extern "C" {
#endif

IVEPlugin IVELoadPlugin(
    const char* path, const char* init_symbol_name, const char* interface_name
);
void* IVEGetPluginAttr(IVEPlugin plugin);
void IVEUnloadPlugin(IVEPlugin plugin);
const char* IVEGetPluginErrorMessage(void);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
