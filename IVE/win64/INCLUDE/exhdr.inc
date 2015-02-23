C @(#) $Id: exhdr.inc,v 1.1 2001/08/09 23:20:47 eric Exp $
C $Name:  $
	integer*4 I_PHOTOSENSOR_READING
	parameter (I_PHOTOSENSOR_READING = 0)
	integer*4 I_TIME_STAMP_SECS
	parameter (I_TIME_STAMP_SECS     = 1)
	integer*4 I_STAGE_X_COORD
	parameter (I_STAGE_X_COORD       = 2)
	integer*4 I_STAGE_Y_COORD
	parameter (I_STAGE_Y_COORD       = 3)
	integer*4 I_STAGE_Z_COORD
	parameter (I_STAGE_Z_COORD       = 4)
	integer*4 I_MIN_INTEN
	parameter (I_MIN_INTEN           = 5)
	integer*4 I_MAX_INTEN
	parameter (I_MAX_INTEN           = 6)
	integer*4 I_MEAN_INTEN
	parameter (I_MEAN_INTEN          = 7)
	integer*4 I_EXP_TIME
	parameter (I_EXP_TIME            = 8)
	integer*4 I_ND_FILTER
	parameter (I_ND_FILTER           = 9)
	integer*4 I_EX_WAVELEN
	parameter (I_EX_WAVELEN          = 10)
	integer*4 I_EM_WAVELEN
	parameter (I_EM_WAVELEN          = 11)
	integer*4 I_INTEN_SCALING
	parameter (I_INTEN_SCALING       = 12)
	integer*4 I_ENERGY_CONV_FACTOR
	parameter (I_ENERGY_CONV_FACTOR  = 13)
        integer*4 I_PHOTOSENSOR_NCONV
        parameter (I_PHOTOSENSOR_NCONV   = 14)
