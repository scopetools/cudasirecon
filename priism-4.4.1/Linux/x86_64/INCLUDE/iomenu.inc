	integer IOMENU_REG_CTRL_ALL, IOMENU_REG_CTRL_5D_SIZE
	integer IOMENU_REG_CTRL_3D_SIZE
	integer IOMENU_REG_CTRL_XY_ONLY
	integer IOMENU_REG_CTRL_4D_NO_Z
	integer IOMENU_REG_CTRL_4D_NO_W
	integer IOMENU_REG_CTRL_4D_NO_T
	parameter (IOMENU_REG_CTRL_ALL = 0)
	parameter (IOMENU_REG_CTRL_5D_SIZE = 1)
	parameter (IOMENU_REG_CTRL_3D_SIZE = 2)
	parameter (IOMENU_REG_CTRL_XY_ONLY = 3)
	parameter (IOMENU_REG_CTRL_4D_NO_Z = 4)
	parameter (IOMENU_REG_CTRL_4D_NO_W = 5)
	parameter (IOMENU_REG_CTRL_4D_NO_T = 6)

	structure /iomenu_struct/
	  integer*4	itypes(2,20),ivals(40),ifiletype(11),iuser(20)
	  real*4	rvals(20),ruser(20)
	  character*256	cvals(10)
	  character*80  labels(20),filelabels(10)
	  character*256	filenames(10)
	  character*80  filefilter(10)
	  character*10	xyztw_size(9)
	  integer*4 	iparm,nxyztw(5),ixyztw(3,9),mode
	  real*4	dmin,dmax,dmean,delxyz(3)
        end structure
        record /iomenu_struct/iom
