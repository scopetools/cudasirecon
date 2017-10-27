#ifndef UCSF_MSG_IPCONSTANTS_H
#define UCSF_MSG_IPCONSTANTS_H

/* These are obselete and are only retained for backward compatibility. */
#define MAX_IN_FILE        7
#define MAX_OUT_FILE       3
#define MAX_FILE           (MAX_IN_FILE+MAX_OUT_FILE+1)

#define OPEN_FILE          (MAX_FILE+1)
#define REG_DIS_CHG        (MAX_FILE+2)
#define UN_REG_DIS_CHG     (MAX_FILE+3)
#define CUS_MENU           (MAX_FILE+4)
#define INIT_PROC          (MAX_FILE+5)
#define PROC_FUNC          (MAX_FILE+6)
#define CLEAN_UP_PROC      (MAX_FILE+7)
#define CHECK_PARAMS       (MAX_FILE+8)
#define INTERRUPT          (MAX_FILE+9)
#define NEW_REGION         (MAX_FILE+10)
#define OPT_MENU           (MAX_FILE+11)
#define CUST_DOIT          (MAX_FILE+12)
#define WAVE_FUNC          (MAX_FILE+13)
#define NEW_XYAREA         (MAX_FILE+14)
#define AFTER_DO_IT        (MAX_FILE+15)
/*
 * These two are not used.  ENV_CMD is obsolete.  Have no record of what
 * VIEW_FILE was intended to do.
 */
#define ENV_CMD            (MAX_FILE+16)
#define VIEW_FILE          (MAX_FILE+17)
#define CHANGE_MODE        (MAX_FILE+18)
#define UNDO_CHECK         (MAX_FILE+19)
#define CMD_LINE_FUNC      (MAX_FILE+20)
#define USAGE_FUNC         (MAX_FILE+21)
#define WRITE_CMD          (MAX_FILE+22)
#define PROMPT_FUNC        (MAX_FILE+23)
#define NEW_RES            (MAX_FILE+24)
#define NUM_CUS_FUNC       (MAX_FILE+25)

#define SAME_AS_FIRST_INPUT  -1
#define TWO_D                0
#define THREE_D              1
#define FOUR_D               2
#define FOUR_D_WAVE          3
#define THREE_D_XYT          4

#define IP_LOG_ERROR         10
#define IP_LOG_WARN          11
#define IP_LOG_INFO          12
#define IP_LOG_PROGRESS      13
#define IP_LOG_DEBUG         14

#define True                 1
#define False                0

static int          zero = 0, one = 1, two = 2, three = 3, four = 4, five = 5;
static int          six = 6, seven = 7, eight = 8, nine = 9, ten = 10;


#endif /* include guard */
