#ifndef UCSF_MSG_RENDERER_PLUGIN_H
#define UCSF_MSG_RENDERER_PLUGIN_H

#include "ive_plugin.h"  /* IVEPlugin */


typedef enum CompositingMethod {
    COMPOSITING_UNKNOWN = -1,
    COMPOSITING_ADDITIVE,
    COMPOSITING_MAX,
    COMPOSITING_PROGRESSIVE,
    COMPOSITING_RGBA,
    COMPOSITING_SLICER
} CompositingMethod;

struct VolumeStorageImpl;
typedef struct VolumeStorageImpl* VolumeStorage;

struct VolumeRenderingContextImpl;
typedef struct VolumeRenderingContextImpl* VolumeRenderingContext;

struct VolumeRenderingResultImp;
typedef struct VolumeRenderingResultImpl* VolumeRenderingResult;

typedef enum VolumeFieldFormat {
    FIELD_UNSIGNED_INT, FIELD_SIGNED_FRACTION, FIELD_UNSIGNED_FRACTION
} VolumeFieldFormat;

typedef struct VolumeLight {
  /* Valid values are between 0 and 1. */
  float intensity;
  /* Ranges from -180 to 180; 0 is positive x axis; 90 is positive y axis. */
  float longitude;
  /*
   * Ranges from -90 to 90; 0 lies in the xy plane; 90 is the positive z axis.
   */
  float latitude;
} VolumeLight;

typedef struct RendererPlugin {
    VolumeStorage (*volume_creator)(
	int nx,
	int ny,
	int nz,
	const float spacing[3],
	int field_count,
	int field_size,
	VolumeFieldFormat format
    );

    void (*volume_destroyer)(VolumeStorage volume);

    void (*volume_loader)(
	VolumeStorage volume,
	int range[3][2],
	int nx_src,
	unsigned char* src[]
    );

    void (*volume_loader_rgba)(
	VolumeStorage volume,
	int range[3][2],
	int nx_src,
	unsigned char* red[],
	unsigned char* green[],
	unsigned char* blue[],
	unsigned char* alpha[]
    );

    VolumeRenderingContext (*context_creator)(
	CompositingMethod compositing_method
    );

    void (*context_destroyer)(VolumeRenderingContext context);

    int (*alpha_lookup_changer)(
	VolumeRenderingContext context,
	CompositingMethod compositing_method,
	const float* table
    );

    int (*light_changer)(
	VolumeRenderingContext context,
	int light_count,
	const VolumeLight lights[],
	int use_world_coord
    );

    int (*reflection_property_changer)(
	VolumeRenderingContext context,
	double diffuse_coeff,
	double specular_coeff,
	double emissive_coeff,
	double specular_exp
    );

    int (*gradient_modulation_changer)(
	VolumeRenderingContext context,
	int is_diffuse_modulated,
	int is_specular_modulated,
	int is_emission_modulated,
	int is_opacity_modulated,
	const double* grad_to_mod_tbl
    );

    int (*ray_termination_changer)(
	VolumeRenderingContext context, CompositingMethod compositing_method
    );

    int (*sampling_changer)(VolumeRenderingContext context, double sampling);

    int (*default_gradient_modulation_table_filler)(double* table);

    VolumeRenderingResult (*result_buffer_creator)(
	CompositingMethod method, int nx_out, int ny_out
    );

    void (*result_buffer_destroyer)(VolumeRenderingResult buffer);

    void (*renderer)(
	VolumeStorage volume,
	VolumeRenderingContext context,
	const int xyzpar[6],
	const float angle[3],
	float inv_zoom,
	CompositingMethod method,
	VolumeRenderingResult buffer
    );

    void (*result_extractor)(
	VolumeRenderingResult buffer,
	CompositingMethod compositing_method,
	float* result,
	float mmm[3],
	float* p_nonzero_min
    );

    void (*result_extractor_rgb)(
	VolumeRenderingResult buffer,
	unsigned char* red,
	float red_mmm[3],
	float* p_red_nonzero_min,
	unsigned char* green,
	float green_mmm[3],
	float* p_green_nonzero_min,
	unsigned char* blue,
	float blue_mmm[3],
	float* p_blue_nonzero_min
    );

    const char* name;
    int max_light_count;
    int alpha_lookup_table_size;
    int gradient_modulation_table_size;
} RendererPlugin;

/*
 * Signature for the routine to in the plugin to initialize the plugin
 * structure.  The function should return zero if successful and non-zero
 * if a failure occurred.  The expected symbol name of the routine is
 * intialize_renderer_plugin.
 */
typedef int (*PluginInitializer)(RendererPlugin* plugin);


#ifdef __cplusplus
extern "C" {
#endif

IVEPlugin load_rendering_plugin(const char* name);
RendererPlugin* get_rendering_plugin_attr(IVEPlugin plugin);
void unload_rendering_plugin(IVEPlugin plugin);
const char* return_rendering_plugin_load_error(void);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
