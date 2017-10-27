#ifndef UCSF_MSG_IWMACRO_H
#define UCSF_MSG_IWMACRO_H

/*-----------------------------------------------------------------------------
   Copyright (C) 1993
    Macromolecular Structure Group of Biochemistry Dept. at University of
    California at San Francisco.
    These coded instructions, statements, and computer programs comprise
    unpublished proprietary information of the Macromolecular Structure Group of
    the Biochemistry Department at University of California at San Francisco,
    and are protected by Federal copyright law.  They may not be disclosed
    to third parties, copied or duplicated in any form - in whole or in part -
    without the prior written consent of Macromolecular Structure Group of
    Biochemistry Department at University of California at San Francisco.
-------------------------------------------------------------------------*/

#include "kernel/IKInterleave.h"
#include "ive_shm.h"
#include "IWApiConstants.h"

/* SHM level */
#define IMSLP                           IW_iwl_ptr
#define IMSLcompat_major()              IMSLP->compat_major
#define IMSLcompat_minor()              IMSLP->compat_minor
#define IMSLWL_array                    IMSLP->WL_array
#define IMSLmax_window()                IMSLP->max_window
#define IMSLtop_window()                IMSLP->top_window
#define IMSLbot_window()                IMSLP->bot_window
#define IMSLnum_window()                IMSLP->num_window
#define IMSLfocus_lock()                IMSLP->focus_lock
#define IMSLfocus_window()              IMSLP->focus_window
#define IMSLnum_displayed()             IMSLP->num_displayed
#define IMSLscreen_left()               IMSLP->screen_left
#define IMSLscreen_bottom()             IMSLP->screen_bottom
#define IMSLscreen_width()              IMSLP->screen_width
#define IMSLscreen_height()             IMSLP->screen_height
#define IMSLscale_max()                 IMSLP->scale_max
#define IMSLscale_min()                 IMSLP->scale_min
#define IMSLgraphics_system()           IMSLP->graphics_system
#define IMSLpseudo_colormap()           IMSLP->pseudo_colormap
#define IMSLind_bk_win()                IMSLP->ind_bk_win
#define IMSLtop_shared_block()          IMSLP->top_shared_block
#define IMSLgraphics_id()               IMSLP->graphics_id
#define IMSLrgb_overlay_color(i)        IMSLP->rgb_overlay_color[i]
#define IMSLimg_mem_size()              IMSLP->img_mem_size
#define IMSLtotal_img_size()            IMSLP->total_img_size
#define IMSLmem_queue(i)                IMSLP->mem_queue[i]
#define IMSLque_pri(i)                  IMSLP->que_pri[i]
#define IMSLmem_lock()                  IMSLP->mem_lock
#define IMSLmem_mode()                  IMSLP->mem_mode
#define IMSLpg_thresh()                 IMSLP->pg_thresh
#define IMSLnum_obj_lists()             IMSLP->num_obj_lists
#define IMSLobj_lists()                 IMSLP->obj_lists
#define IMSLblock_lock()                IMSLP->block_lock

/* Window level */
#define IMWLP(w)                        ((IW_WL_PTR) IVERepToPtr(IMSLWL_array, IW_arena))[w]
#define IMWLprev(w)                     IMWLP(w).prev
#define IMWLnext(w)                     IMWLP(w).next
#define IMWLscr_x_pos(w)                IMWLP(w).scr_x_pos
#define IMWLscr_y_pos(w)                IMWLP(w).scr_y_pos
#define IMWLmax_width(w)                IMWLP(w).max_width
#define IMWLmax_height(w)               IMWLP(w).max_height
#define IMWLwin_x_pos(w)                IMWLP(w).win_x_pos
#define IMWLwin_y_pos(w)                IMWLP(w).win_y_pos
#define IMWLwin_width(w)                IMWLP(w).win_width
#define IMWLwin_height(w)               IMWLP(w).win_height
#define IMWLmodified(w)                 IMWLP(w).modified
/*
#define IMWLLUT(w)                      IMWLP(w).LUT
*/
#define IMWLdis_mode(w)                 IMWLP(w).dis_mode
/*
#define IMWLdis_buffer(w)               IMWLP(w).dis_buffer
#define IMWLwrite_buffer(w)             IMWLP(w).write_buffer
*/
#define IMWLexist(w)                    IMWLP(w).exist
#define IMWLwave_off(w,b)               IMWLP(w).bank_off[b]
#define IMWLwave_color(w,b)             IMWLP(w).bank_color[b]
#define IMWLdata_volume_gil_off(w)      IMWLP(w).data_volume_gil_off
#define IMWLil_array(w)                 IMWLP(w).il_array
#define IMWLpseudo_graphics_color(w)    IMWLP(w).pseudo_graphics_color
#define IMWLgraphics_id(w)              IMWLP(w).graphics_id
#define IMWLgraphics_list(w)            IMWLP(w).graphics_list
#define IMWLgraphics_displayed(w)       IMWLP(w).graphics_displayed
#define IMWLgr_z_range(w,i)             IMWLP(w).gr_z_range[i]
#define IMWLgr_t_range(w,i)             IMWLP(w).gr_t_range[i]
/*
#define IMWLgraphics_added_to_front(w)  IMWLP(w).graphics_added_to_front
*/
#define IMWLclear_bkg(w)                IMWLP(w).clear_bkg
#define IMWLborder_style(w)             (3 & IMWLP(w).border_style)
#define IMSet_border_style(w, style)    (IMWLP(w).border_style = \
    (3 & (style)) | (12 & IMWLP(w).border_style))
#define IMWLtool_style(w)               ((12U & IMWLP(w).border_style) >> 2)
#define IMSet_tool_style(w, style)      (IMWLP(w).border_style = \
    (3 & IMWLP(w).border_style) | (12 & ((style) << 2)))
/*
#define IMWLauto_size(w)                IMWLP(w).auto_size
*/
#define IMWLmd_x_off(w)                 IMWLP(w).md_x_off
#define IMWLmd_y_off(w)                 IMWLP(w).md_y_off
#define IMWLmd_dis_off(w)               IMWLP(w).md_dis_off
#define IMWLiconified(w)                IMWLP(w).iconified
#define IMWLdelay(w)                    IMWLP(w).delay
#define IMWLinitial_delay(w)            IMWLP(w).initial_delay
#define IMWLmonitor_pid(w)              IMWLP(w).monitor_pid   
#define IMWLmonitor_wait(w)             IMWLP(w).monitor_wait   
#define IMWLxwindow(w,i)                IMWLP(w).xwindow[i]
#define IMWLmonitor_font(w)             IMWLP(w).monitor_font
#define IMWLsync_enter_lock(w)          IMWLP(w).sync_enter_lock
#define IMWLdata_sync_lock(w)           IMWLP(w).data_sync_lock
#define IMWLscratch_mode(w)             IMWLP(w).scratch_mode
#define IMWLcur_wave_num(w)             IMWLP(w).cur_wave_num
#define IMWLzoom(w)                     IMWLP(w).zoom
/*
#define IMWLselected(w)                 IMWLP(w).selected
*/
#define IMWLscale_gr_off(w)             IMWLP(w).scale_gr_off
#define IMWLscale_bar_vertical(w)       IMWLP(w).scale_bar_vertical
#define IMWLscale_bar_length(w)         IMWLP(w).scale_bar_length
#define IMWLview_port(w)                IMWLP(w).view_port
#define IMWLfirst_reg_event(w)          IMWLP(w).first_reg_event
#define IMWLfirst_reg_dis_chg(w)        IMWLP(w).first_reg_dis_chg
#define IMWLmouse_button(w,i)           IMWLP(w).mouse_button[i]
#define IMWLmouse_button_ext(w,i)       IMWLP(w).mouse_button_ext[i]
/*
#define IMWLrotate_status(w)            IMWLP(w).rotate_status
*/
#define IMWLstep_type(w)                IMWLP(w).step_type
#define IMWLinterpolation(w)            IMWLP(w).interpolation
#define IMWLstep_inc(w)                 IMWLP(w).step_inc
/*
#define IMWLstride(w)                   IMWLP(w).stride
#define IMWLgap(w)                      IMWLP(w).gap
*/
#define IMWLmult_disp_x(w)              IMWLP(w).mult_disp_x
#define IMWLmult_disp_y(w)              IMWLP(w).mult_disp_y
#define IMWLmat_layout(w)               IMWLP(w).mat_layout
#define IMWLapply_to_wave(w)            IMWLP(w).apply_to_wave
#define IMWLproc_id(w)                  IMWLP(w).proc_id
#define IMWLcomplex_disp(w)             IMWLP(w).complex_disp
#define IMWLscaling_algorithm(w)        IMWLP(w).scaling_algorithm
#define IMWLimg_displayed(w)            IMWLP(w).img_displayed
#define IMWLdisp_size(w)                IMWLP(w).disp_size
#define IMWLmem_size(w)                 IMWLP(w).mem_size
#define IMWLview_flag(w)                IMWLP(w).view_flag
#define IMWLfile_name(w)                IMWLP(w).file_name
#define IMWLlocal_disk(w)               IMWLP(w).local_disk
#define IMWLcur_sec_size(w)             IMWLP(w).cur_sec_size
#define IMWLsave_buf_stat(w)            IMWLP(w).save_buf_stat
#define IMWLgr_obj_header_list(w)       IMWLP(w).gr_obj_header_list
#define IMWLstereo_offset(w)            IMWLP(w).stereo_offset
#define IMWLdisplayed_resolution(w)     IMWLP(w).displayed_resolution
#define IMWLsection_number_style(w)     IMWLP(w).section_number_style
#define IMWLsection_number_color(w)     IMWLP(w).section_number_color

/* Window->geometry structure */
#define IMWGP(w)                        ((IW_GIL_PTR) IVERepToPtr(IMWLdata_volume_gil_off(w), IW_arena))
#define IMWLGLx_size(w)                 IMWGP(w)->x_size
#define IMWLGLy_size(w)                 IMWGP(w)->y_size
#define IMWLGLz_size(w)                 IMWGP(w)->z_size
#define IMWLGLmode(w)                   IMWGP(w)->mode
#define IMWLGLld_x_st(w)                IMWGP(w)->ld_x_st
#define IMWLGLld_y_st(w)                IMWGP(w)->ld_y_st
#define IMWLGLld_z_st(w)                IMWGP(w)->ld_z_st
#define IMWLGLld_t_st(w)                IMWGP(w)->ld_t_st
#define IMWLGLfile_mxyz(w,i)            IMWGP(w)->file_mxyz[i]
#define IMWLGLxyz_len(w,i)              IMWGP(w)->xyz_len[i]
/*
#define IMWLGLx_len(w)                  IMWGP(w)->x_len
#define IMWLGLy_len(w)                  IMWGP(w)->y_len
#define IMWLGLz_len(w)                  IMWGP(w)->z_len
*/
#define IMWLGLfile_angle(w,i)           IMWGP(w)->file_angle[i]
#define IMWLGLfile_mapc(w)              IMWGP(w)->file_mapc
#define IMWLGLfile_mapr(w)              IMWGP(w)->file_mapr
#define IMWLGLfile_maps(w)              IMWGP(w)->file_maps
#define IMWLGLfile_space_group(w)       IMWGP(w)->file_space_group
#define IMWLGLfile_ext_header_size(w)   IMWGP(w)->file_ext_header_size
#define IMWLGLfile_num_ints(w)          IMWGP(w)->file_num_ints
#define IMWLGLfile_num_reals(w)         IMWGP(w)->file_num_reals
#define IMWLGLnum_resolutions(w)        IMWGP(w)->num_resolutions
#define IMWLGLfile_nzfact(w)            IMWGP(w)->file_nzfact
#define IMWLGLfile_type(w)              IMWGP(w)->file_type
#define IMWLGLfile_lens(w)              IMWGP(w)->file_lens
#define IMWLGLtilt_set_axis(w)          IMWGP(w)->tilt_set_axis
#define IMWLGLtilt_set_angle(w)         IMWGP(w)->tilt_set_angle
#define IMWLGLsects_in_stereo(w)        IMWGP(w)->sects_in_stereo
#define IMWLGLsects_between_pair(w)     IMWGP(w)->sects_between_pair
#define IMWLGLx_tilt_angle(w)           IMWGP(w)->x_tilt_angle
#define IMWLGLy_tilt_angle(w)           IMWGP(w)->y_tilt_angle
#define IMWLGLz_tilt_angle(w)           IMWGP(w)->z_tilt_angle
#define IMWLGLnum_labels(w)             IMWGP(w)->num_labels
#define IMWLGLlabel(w)                  IMWGP(w)->label
#define IMWLGLnocopy(w)                 IMWGP(w)->nocopy
#define IMWLGLfile_ext_header(w)        IMWGP(w)->file_ext_header
#define IMWLGLfile_dec_format(w)        IMWGP(w)->file_dec_format
#define IMWLGLhdr_format(w)             IMWGP(w)->hdr_format
#define IMWLGLfile_name(w,wave)         IMWGP(w)->file_name[wave]
#define IMWLGLvalid_ld_info(w)          IMWGP(w)->valid_ld_info
#define IMWLGLld_x_end(w)               IMWGP(w)->ld_x_end
#define IMWLGLld_x_bin(w)               IMWGP(w)->ld_x_bin
#define IMWLGLld_y_end(w)               IMWGP(w)->ld_y_end
#define IMWLGLld_y_bin(w)               IMWGP(w)->ld_y_bin
#define IMWLGLld_z_end(w)               IMWGP(w)->ld_z_end
#define IMWLGLld_z_bin(w)               IMWGP(w)->ld_z_bin
#define IMWLGLld_z_inc(w)               IMWGP(w)->ld_z_inc
#define IMWLGLld_time_st(w)             IMWGP(w)->ld_time_st
#define IMWLGLld_time_end(w)            IMWGP(w)->ld_time_end
#define IMWLGLld_time_bin(w)            IMWGP(w)->ld_time_bin
#define IMWLGLld_time_inc(w)            IMWGP(w)->ld_time_inc
#define IMWLGLld_wave(w,i)              IMWGP(w)->ld_wave[i]
#define IMWLGLscaled(w)                 IMWGP(w)->scaled
#define IMWLGLscale_min(w,wave)         IMWGP(w)->scale_min[wave]
#define IMWLGLscale_max(w,wave)         IMWGP(w)->scale_max[wave]
#define IMWLGLscale_coeff1(w,wave)      IMWGP(w)->scale_coeff1[wave]
#define IMWLGLscale_coeff2(w,wave)      IMWGP(w)->scale_coeff2[wave]
#define IMWLGLld_st_end_erase(w)        IMWGP(w)->ld_st_end_erase
#define IMWLGLld_wave_wid(w,wave)       IMWGP(w)->ld_wave_wid[wave]
#define IMWLGLld_wave_bank(w,wave)      IMWGP(w)->ld_wave_bank[wave]
#define IMWLGLld_wave_bank_pos(w,wave)  IMWGP(w)->ld_wave_bank_pos[wave]
#define IMWLGLprev(w)                   IMWGP(w)->prev
#define IMWLGLnext(w)                   IMWGP(w)->next
#define IMWLGLnum_links(w)              IMWGP(w)->num_links
#define IMWLGLlocal_file(w)             IMWGP(w)->local_file
#define IMWLGLset_local_file(w)         IMWGP(w)->set_local_file
#define IMWLGLinten_thresh(w,wave)      IMWGP(w)->inten_thresh[wave]
#define IMWLGLinten_thresh_on(w,wave)   IMWGP(w)->inten_thresh_on[wave]

/* SB */
#define IMWLSBP(w)                      ((IW_GRL_PTR) IVERepToPtr(IMWLscale_gr_off(w), IW_arena))
#define IMWLSBdisplayed(w)              IMWLSBP(w)->displayed
#define IMWLSBmarked_for_delete(w)      IMWLSBP(w)->marked_for_delete
#define IMWLSBprimitive_ptr(w)          ((IW_GR_LINE*) IVERepToPtr(IMWLSBP(w)->primitive_off, IW_arena))
#define IMWLSBid(w)                     IMWLSBprimitive_ptr(w)->id
#define IMWLSBcolor(w)                  IMWLSBprimitive_ptr(w)->color
#define IMWLSBx1(w)                     IMWLSBprimitive_ptr(w)->x1
#define IMWLSBy1(w)                     IMWLSBprimitive_ptr(w)->y1
#define IMWLSBx2(w)                     IMWLSBprimitive_ptr(w)->x2
#define IMWLSBy2(w)                     IMWLSBprimitive_ptr(w)->y2
#define IMWLSBthick(w)                  IMWLSBprimitive_ptr(w)->thickness

/* Window->view port structure */
#define IMWVP(w)                        IMWLview_port(w)
#define IMWLVProtation_matrix(w)        IMWVP(w).rotation_matrix      
#define IMWLVPoffset_vector(w,o)        IMWVP(w).offset_vector[o]
#define IMWLVPz_incr(w)                 IMWVP(w).z_incr
#define IMWLVPdimensions(w,dim)         IMWVP(w).dimensions[dim]
#define IMWLVPdelta(w,del)              IMWVP(w).delta[del]
#define IMWLVPslice_mode(w)             IMWVP(w).slice_mode

/* bank level */
#define IMBLO(w,b)                      IMWLbank_off(w,b)
#define IMBLP(w,b)                      ((IW_BAL_PTR) IVERepToPtr(IMWLbank_off(w,b), IW_arena))
#define IMBLnum_imgs(w,b)               IMBLP(w,b)->num_imgs
#define IMBLmapped(w,b)                 IMBLP(w,b)->mapped
#define IMBLgr_color(w,b)               IMBLP(w,b)->gr_color
#define IMBLcur_z_num(w,b)              IMBLP(w,b)->cur_z_num
#define IMBLcur_time_num(w,b)           IMBLP(w,b)->cur_time_num
#define IMBLx_offset(w,b)               IMBLP(w,b)->x_offset
#define IMBLy_offset(w,b)               IMBLP(w,b)->y_offset
#define IMBLoff_group(w,b)              IMBLP(w,b)->off_group

/* image level */
#define IMILO(w,i)                      (((IVEPtrRep*) IVERepToPtr(IMWLil_array(w), IW_arena))[i])
#define IMILP(w,i)                      ((IW_IAL_PTR) IVERepToPtr(IMILO(w,i), IW_arena))
#define IMILscale_min(w,i)              IMILP(w,i)->scale_min
#define IMILscale_max(w,i)              IMILP(w,i)->scale_max
#define IMILscale_coeff1(w,i)           IMILP(w,i)->scale_coeff1
#define IMILscale_coeff2(w,i)           IMILP(w,i)->scale_coeff2
#define IMILmin(w,i)                    IMILP(w,i)->min
#define IMILmax(w,i)                    IMILP(w,i)->max
#define IMILmean(w,i)                   IMILP(w,i)->mean
#define IMILbuf_stat(w,i)               IMILP(w,i)->buf_stat
#define IMILimg_mem_off(w,i)            IMILP(w,i)->img_mem_off
#define IMILbuf_mem_off(w,i)            IMILP(w,i)->buf_mem_off
#define IMILimg_mem_node(w,i)           IMILP(w,i)->img_mem_node
#define IMILbuf_mem_node(w,i)           IMILP(w,i)->buf_mem_node
#define IMILsec_loc(w,i)                IMILP(w,i)->sec_loc
#define IMILgraphics_list(w,i)          IMILP(w,i)->graphics_list
#define IMILinten_thresh(w,i)           IMILP(w,i)->inten_thresh
#define IMILinten_thresh_on(w,i)        IMILP(w,i)->inten_thresh_on

/* Misc */

#define IMWLbank_off(w,b)               IMWLP(w).bank_off[b]
#define IMWLbank_color(w,b)             IMWLP(w).bank_color[b]
#define IMWLgraphics_list_size(w)       IMWLP(w).graphics_list_size
#define IMWLnum_registered_events(w)    IMWLP(w).num_registered_events
#define IMWLnum_registered_dis_chgs(w)  IMWLP(w).num_registered_dis_chgs
#define IMWLslice_mode(w)               IMWLP(w).slice_mode
#define IMWLimg_width(w)                IMWLP(w).img_width
#define IMWLimg_height(w)               IMWLP(w).img_height
#define IMWLimg_mode(w)                 IMWLP(w).img_mode

#define IMmd_x_offset(w,b)              (IMBLx_offset(w,b)+IMWLmd_x_off(w))
#define IMmd_y_offset(w,b)              (IMBLy_offset(w,b)+IMWLmd_y_off(w))
#define IMmd_dis_sec(w,b)               (IMBLcur_z_num(w,b)+IMWLmd_dis_off(w))

#define IMisWinID(w)                   ((w > 0) && (w < IMSLmax_window()))
#define IMisWinExist(w)                (IMisWinID(w) && IMWLexist(w))
#define IMisWinInDL(w)                 ((IMWLprev(w) != IW_EMPTY) || (IMWLnext(w) != IW_EMPTY) || (IMSLtop_window() == w))

#define IMAddNode(head,node)           node->prev = NULL; \
                                       node->next = head; \
                                       if (head) \
                                           head->prev = node; \
                                       head = node;
#define IKSetWaveModified(w)           IMBLmodified(w, 0) = IW_TRUE; \
                                       IMBLmodified(w, 1) = IW_TRUE; \
                                       IMBLmodified(w, 2) = IW_TRUE; \
                                       IMBLmodified(w, 3) = IW_TRUE; \
                                       IMBLmodified(w, 4) = IW_TRUE;

#define IMGetZoomedRatio(w)            (int)(1.0 / IMWLzoom(w) + 0.499)
#define IMGetComplexNormalValue(ptr)   sqrtf((ptr)->real*(ptr)->real+(ptr)->imaginary*(ptr)->imaginary)

#define IMScaleImg(src,dst,scl_ar,min,max,s_min,s_max,co1,scl)\
            if (src <= min) \
                dst = s_min; \
            else if (src >= max) \
                dst = s_max; \
            else  \
                if (co1 == 0 || co1 == 1) \
                    dst = (unsigned char)((src-min)*scl+s_min); \
                else  \
                    dst = scl_ar[(int)((src-min)*scl)]; 

#define IMScaleImgA(s_ar,d_ar,scl_ar,min,max,s_min,s_max,co1,scl,i,j,ind,hi,wd)\
    for (i = 0; i < hi; i++) { \
        for (j = 0; j < wd; j++) { \
            ind = i * wd + j; \
            IMScaleImg(s_ar[ind],d_ar[ind],scl_ar,min,max,s_min,s_max,co1,scl);\
        } \
    }
#define IMScaleImgB(s_ar,d_ar,scl_ar,min,max,s_min,s_max,co1,scl,i,j,ind,hi,wd)\
    for (i = 0; i < hi; i++) { \
        for (j = 0; j < wd; j++) { \
            ind = i * wd + j; \
            IMScaleImg(sqrtf((&s_ar[ind])->real*(&s_ar[ind])->real+(&s_ar[ind])->imaginary*(&s_ar[ind])->imaginary),d_ar[ind],scl_ar,min,max,s_min,s_max,co1,scl);\
        } \
    }
#define IMScaleImgC(s_ar,d_ar,scl_ar,min,max,s_min,s_max,co1,scl,i,j,ind,hi,wd)\
    for (i = 0; i < hi; i++) { \
        for (j = 0; j < wd; j++) { \
            ind = i * wd + j; \
            IMScaleImg((&s_ar[ind])->real,d_ar[ind],scl_ar,min,max,s_min,s_max,co1,scl);\
        } \
    }
#define IMScaleImgD(s_ar,d_ar,scl_ar,min,max,s_min,s_max,co1,scl,i,j,ind,hi,wd)\
    for (i = 0; i < hi; i++) { \
        for (j = 0; j < wd; j++) { \
            ind = i * wd + j; \
            IMScaleImg((&(s_ar[ind]))->imaginary,d_ar[ind],scl_ar,min,max,s_min,s_max,co1,scl);\
        } \
    }
#define IMScaleImgE(s_ar,d_ar,scl_ar,min,max,s_min,s_max,co1,scl,i,j,ind,hi,wd)\
    for (i = 0; i < hi; i++) { \
        for (j = 0; j < wd; j++) { \
            ind = i * wd + j; \
            IMScaleImg(fatan2((&s_ar[ind])->imaginary,(&s_ar[ind])->real),d_ar[ind],scl_ar,min,max,s_min,s_max,co1,scl);\
        } \
    }


#endif /* include guard */
