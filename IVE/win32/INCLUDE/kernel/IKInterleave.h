#ifndef UCSF_MSG_IKINTERLEAVE_H
#define UCSF_MSG_IKINTERLEAVE_H

/* @(#) $Id: IKInterleave.h,v 1.1 2001/08/09 23:20:50 eric Exp $ */
/* $Name:  $ */
/*
 * These are kernel / monitor helper functions to go back and forth between
 * zwt indices and single section index.  Other applications should use
 * IMRtSecNum and IMRtZWTNum.
 */


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
 * Returns 1 if intlv (the code for how z, waves, and times are interleaved)
 * is invalid.  Otherwise returns 0 and sets *p_isec to the 1D section index.
 * Assumes that nz, nw, and nt are positive integers and that iz, iw, and
 * it are non-negative integers less than nz, nw, and nt respectively.
 */
int IKZWT2Sec(
    int* p_isec, int iz, int iw, int it, int nz, int nw, int nt, int intlv
);


/*
 * Returns 1, if intlv is invalid.  Otherwise returns 0 and sets *p_iz, *p_iw,
 * *p_it to the z, wave, and time indexes corresponding to isec.  Assumes that
 * nz, nw, and nt are positive integers and that isec is a non-negative
 * integer less than nz * nw * nt.
 */
int IKSec2ZWT(
    int* p_iz,
    int* p_iw,
    int* p_it,
    int isec,
    int nz,
    int nw,
    int nt,
    int intlv
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
