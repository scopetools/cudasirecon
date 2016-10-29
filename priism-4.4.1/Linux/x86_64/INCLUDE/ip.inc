c     MAX_IN_FILE, MAX_OUT_FILE, and MAX_FILE are obselete and
c     only retained for backward compatibility.  VIEW_FILE and
c     ENV_CMD are not used.
      integer MAX_IN_FILE, MAX_OUT_FILE, MAX_FILE
      integer OPEN_FILE, REG_DIS_CHG, UN_REG_DIS_CHG
      integer CUS_MENU, INIT_PROC, PROC_FUNC, CLEAN_UP_PROC
      integer CHECK_PARAMS, INTERRUPT, NEW_REGION, OPT_MENU
      integer CUST_DOIT, WAVE_FUNC, NEW_XYAREA, AFTER_DO_IT
      integer ENV_CMD, VIEW_FILE, CHANGE_MODE, UNDO_CHECK
      integer CMD_LINE_FUNC, USAGE_FUNC, WRITE_CMD, PROMPT_FUNC
      integer NUM_CUS_FUNC
      parameter (MAX_IN_FILE       =  7)
      parameter (MAX_OUT_FILE      =  3)
      parameter (MAX_FILE          =  MAX_IN_FILE + MAX_OUT_FILE + 1)
      parameter (OPEN_FILE         =  MAX_FILE + 1)
      parameter (REG_DIS_CHG       =  MAX_FILE + 2)
      parameter (UN_REG_DIS_CHG    =  MAX_FILE + 3)
      parameter (CUS_MENU          =  MAX_FILE + 4)
      parameter (INIT_PROC         =  MAX_FILE + 5)
      parameter (PROC_FUNC         =  MAX_FILE + 6)
      parameter (CLEAN_UP_PROC     =  MAX_FILE + 7)
      parameter (CHECK_PARAMS      =  MAX_FILE + 8)
      parameter (INTERRUPT         =  MAX_FILE + 9)
      parameter (NEW_REGION        =  MAX_FILE + 10)
      parameter (OPT_MENU          =  MAX_FILE + 11)
      parameter (CUST_DOIT         =  MAX_FILE + 12)
      parameter (WAVE_FUNC         =  MAX_FILE + 13)
      parameter (NEW_XYAREA        =  MAX_FILE + 14)
      parameter (AFTER_DO_IT       =  MAX_FILE + 15)
      parameter (ENV_CMD           =  MAX_FILE + 16)
      parameter (VIEW_FILE         =  MAX_FILE + 17)
      parameter (CHANGE_MODE       =  MAX_FILE + 18)
      parameter (UNDO_CHECK        =  MAX_FILE + 19)
      parameter (CMD_LINE_FUNC     =  MAX_FILE + 20)
      parameter (USAGE_FUNC        =  MAX_FILE + 21)
      parameter (WRITE_CMD         =  MAX_FILE + 22)
      parameter (PROMPT_FUNC       =  MAX_FILE + 23)
      parameter (NUM_CUS_FUNC      =  MAX_FILE + 24)

      integer SAME_AS_FIRST_INPUT, TWO_D, THREE_D, FOUR_D, FOUR_D_WAVE
      integer THREE_D_XYT
      parameter (SAME_AS_FIRST_INPUT = -1)
      parameter (TWO_D =  0)
      parameter (THREE_D = 1)
      parameter (FOUR_D = 2)
      parameter (FOUR_D_WAVE = 3)
      parameter (THREE_D_XYT = 4)

      integer IP_LOG_ERROR, IP_LOG_WARN, IP_LOG_INFO, IP_LOG_PROGRESS
      integer IP_LOG_DEBUG
      parameter (IP_LOG_ERROR = 10)
      parameter (IP_LOG_WARN = 11)
      parameter (IP_LOG_INFO = 12)
      parameter (IP_LOG_PROGRESS = 13)
      parameter (IP_LOG_DEBUG = 14)

c     Declare IPGetDataPtr to be return a pointer-sized integer.
      integer*8 IPGetDataPtr
      external IPGetDataPtr
