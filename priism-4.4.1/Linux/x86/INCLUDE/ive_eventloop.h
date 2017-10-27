#ifndef UCSF_MSG_IVE_EVENTLOOP_H
#define UCSF_MSG_IVE_EVENTLOOP_H

/*
 * Copyright(C) 2001, 2007
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Declare types and functions for an Xt-style event loop.
 */


#include "ive_conf.h"
#include "ive_time.h"


enum { IVE_IOEVENT_READ = 1, IVE_IOEVENT_WRITE = 2, IVE_IOEVENT_EXCPT = 4 };

struct IVEEventContextImpl;
typedef struct IVEEventContextImpl* IVEEventContext;

struct IVEInputIdImpl;
typedef struct IVEInputIdImpl* IVEInputId;
typedef void (*IVEInputProc)(
    IVEDesc desc, void* arg, IVEInputId id, IVEEventContext context
);

typedef int IVEMessage;
#define IVE_MESSAGE_MIN (1)
#define IVE_MESSAGE_MAX (255)
typedef void (*IVEMessageProc)(
    void* arg, IVEMessage msg, IVEEventContext context
);

struct IVESignalIdImpl;
typedef struct IVESignalIdImpl* IVESignalId;
typedef void (*IVESignalProc)(
    void* arg, IVESignalId id, IVEEventContext context
);

struct IVETimerIdImpl;
typedef struct IVETimerIdImpl* IVETimerId;
typedef void (*IVETimerProc)(
    void* arg, IVETimerId id, IVEEventContext context
);

struct IVEWorkProcIdImpl;
typedef struct IVEWorkProcIdImpl* IVEWorkProcId;
typedef void (*IVEWorkProc)(
    void* arg, IVEWorkProcId id, IVEEventContext context
);

enum {
    IVE_WATCH_ADD_READ = 1,
    IVE_WATCH_ADD_WRITE = 2,
    IVE_WATCH_ADD_EXCPT = 4,
    IVE_WATCH_DEL_READ = 8,
    IVE_WATCH_DEL_WRITE = 16,
    IVE_WATCH_DEL_EXCPT = 32
};

typedef void (*IVEWatchProc)(IVEDesc desc, int reason, void* arg);


#ifdef __cplusplus
extern "C" {
#endif

IVEEventContext IVECreateEventContext(void);

IVEEventContext IVECreateThreadedEventContext(void);

void IVEDestroyEventContext(IVEEventContext context);

void IVEEventLoop(IVEEventContext context);

int IVEProcessEvent(int block, IVEEventContext context);

void IVEEventReleaseInterest(IVEEventContext context);

IVEInputId IVEAddInput(
    IVEDesc desc,
    int reason,
    IVEInputProc proc,
    void* arg,
    IVEEventContext context
);

void IVERemoveInput(IVEInputId id, IVEEventContext context);

void IVESetMessageHandler(
    IVEMessage msg, IVEMessageProc proc, void* arg, IVEEventContext context
);

void IVECancelMessageHandler(IVEMessage msg, IVEEventContext context);

int IVEPostMessage(IVEMessage msg, IVEEventContext context);

IVESignalId IVEAddSignal(
    IVESignalProc proc, void* arg, IVEEventContext context
);

void IVERemoveSignal(IVESignalId id, IVEEventContext context);

void IVENoticeSignal(IVESignalId id, IVEEventContext context);

IVETimerId IVEAddTimer(
    IVETimeInterval intvl,
    IVETimerProc proc,
    void* arg,
    IVEEventContext context
);

void IVERemoveTimer(IVETimerId id, IVEEventContext context);

IVEWorkProcId IVEAddWorkProc(
    IVEWorkProc proc,
    void* arg,
    IVEEventContext context
);

void IVERemoveWorkProc(IVEWorkProcId id, IVEEventContext context);

void IVERemoveAllEventTriggers(IVEEventContext context);

void IVESetWatchProc(IVEWatchProc proc, void* arg, IVEEventContext context);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
