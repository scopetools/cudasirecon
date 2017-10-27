#ifndef UCSF_MSG_IVEQD_H
#define UCSF_MSG_IVEQD_H

/*
 * Declares a read-only (at this point) interface for getting the
 * IVE/Priism queue definitions.
 */

#include <stdio.h>  /* FILE */


/*
 * An IVEQD serves as a handle to a collection of queues and the global
 * properties of those queues.  For the Fortran interface, the handle is
 * an integer.
 */
struct IVEQDImpl;
typedef struct IVEQDImpl* IVEQD;

/*
 * An IVEQDMember serves as a handle to a single queue.  For the Fortran
 * interface, the handle is an integer.
 */
struct IVEQDMemberImpl;
typedef struct IVEQDMemberImpl* IVEQDMember;

struct IVEQDExpansionContextImpl;
typedef struct IVEQDExpansionContextImpl* IVEQDExpansionContext;

/*
 * Constructs a collection of queues from the file with the given name
 * (if filename is a null pointer or emptry string will attempt to use
 * a user or system-wide default; if the user and system-wide defaults do
 * not exist will fall back to extracting the information from the old-style
 * environment variables).  Returns a handle to the collection for further
 * operations.  When no longer need the handle, you should pass it to
 * IVEQDDeconstruct() to free any associated resources.  The returned
 * handle will be zero if the operation failed.  The equivalent Fortran
 * interface has the prototype
 *
 *     integer function IVEQDConstruct(filename)
 *     character*(*) filename
 *
 * and returns zero when the operation fails.
 */
IVEQD IVEQDConstruct(const char* filename);

/*
 * Frees the resources associated with a handle returned by IVEQDConstruct()
 * and invalidates the handle.  The equivalent Fortran interface has the
 * prototype
 *
 *    subroutine IVEQDDestroy(qd)
 *    integer qd
 */
void IVEQDDestroy(IVEQD qd);

/*
 * For testing or debugging purposes, writes a human readable description
 * of the contents of the handle.  There is no equivalent Fortran interface.
 */
void IVEQDDump(IVEQD qd, FILE* fo);

/*
 * Returns the number of queues in the collection.  The equivalent Fortran
 * interface is
 *
 *     integer function IVEQDGetCount(qd)
 *     integer qd
 */
int IVEQDGetCount(IVEQD qd);

/*
 * These routines query a collection of queues for their global properties.
 * Properties are identified by name and type.  The type is either boolean,
 * integer, string, or a list of strings.  For properties that are list of
 * strings there is a call to return the number of strings in the list and
 * another call to get a member by its zero-based index.  The return value for
 * any of these calls is one of the following:
 *
 *    0    The property is defined for this collection of queues and has
 *         been returned in the value argument.
 *    1    The property is not defined for this collection of queues but
 *         a default is available and has been returned in the value argument.
 *    2    The property is not defined for this collection of queues and
 *         a default is not available.  The value argument has not been
 *         modified.
 *    3    The list element requested does not exist.  The value argument
 *         has not been modified.
 *   -1    An error occurred processing the request.  The value argument
 *         has not been modified.
 *
 * The equivalent Fortran interfaces are
 *
 *     integer function IVEQDQueryBool(qd, name, bvalue)
 *     integer function IVEQDQueryInt(qd, name, ivalue)
 *     integer function IVEQDQueryStringLength(qd, name, ivalue)
 *     integer function IVEQDQueryString(qd, name, cvalue)
 *     integer function IVEQDQuerySListLength(qd, name, ivalue)
 *     integer function IVEQDQuerySListMemberLength(qd, name, ientry, ivalue)
 *     integer function IVEQDQuerySListMember(qd, name, ientry, cvalue)
 *     integer qd, ientry, ivalue
 *     logical bvalue
 *     character*(*) name, cvalue
 *
 * As of this writing, these are the properties that are defined:
 *
 *     Name                  Type    Meaning
 *     DefaultClusterQueue   S       Holds the name of the designated default
 *                                   queue for parallel jobs that may run on
 *                                   multiple hosts.  An empty string indicates
 *                                   that no default has been designated.
 *     DefaultGPUQueue       S       Holds the name of the designated default
 *                                   queue for jobs expecting GPUs to be
 *                                   available as coprocessors.  An empty
 *                                   string indicates that no default has been
 *                                   designated.
 *     DefaultSerialQueue    S       Holds the name of the designated default
 *                                   queue for serial (i.e. single-processor)
 *                                   tasks.  An empty string indicates that
 *                                   no default has been designated.
 *     DefaultThreadedQueue  S       Holds the name of the designated default
 *                                   queue for parallel jobs that can only run
 *                                   on something that is (or acts like) a
 *                                   single host.
 */
int IVEQDQueryBool(IVEQD qd, const char* name, int* p_value);
int IVEQDQueryInt(IVEQD qd, const char* name, int* p_value);
int IVEQDQueryString(IVEQD qd, const char* name, const char** p_value);
int IVEQDQuerySListLength(IVEQD qd, const char* name, int* p_value);
int IVEQDQuerySListMember(
    IVEQD qd, const char* name, int ientry, const char** p_value
);

/*
 * Returns a handle to the "queue" representing the local system.  The
 * equivalent Fortran interface has the prototype
 *
 *     integer function IVEQDLocalMachine(qd)
 *     integer qd
 */
IVEQDMember IVEQDLocalMachine(IVEQD qd);

/*
 * Returns a handle to the nth queue in the collection where n is greater than
 * or equal to zero and less than the return value of IVEQDGetCount().
 * The equivalent Fortran interface has the prototype
 *
 *     integer function IVEQDGetMember(qd, n);
 *     integer qd, n
 */
IVEQDMember IVEQDGetMember(IVEQD qd, int n);

/*
 * Returns one (true in Fortran) if the given IVEQDMember was returned
 * by IVEQDLocalMachine().  Otherwise returns zero (false in Fortran).
 * The equivalent Fortran interface has the prototype
 *
 *     bool function IVEQDMemberIsLocal(qdm)
 *     integer qd
 */
int IVEQDMemberIsLocal(IVEQDMember qdm);

/*
 * These routines query a queue for its properties.  Properties are identified
 * by name and type.  The type is either boolean, integer, string, or a list
 * of strings.  For string properties, there is a call to return the length of
 * the string (not including the terminating null character).  For properties
 * that are list of strings there are calls to return the number of strings in
 * the list, the length of an element (selected by its zero-based index) in
 * the list (not including the terminating null character), or the value of an
 * element (again selected by its zero-based index) in the list.  The return
 * value for any of these calls is one of the following:
 *
 *    0    The property is defined for this queue and has been returned in
 *         the value argument.
 *    1    The property is not defined for this queue but a default is
 *         available and has been returned in the value argument.
 *    2    The property is not defined for this queue and a default is not
 *         available.  The value argument has not been modified.
 *    3    The list element requested does not exist.  The value argument
 *         has not been modified.
 *   -1    An error occurred processing the request.  The value argument
 *         has not been modified.
 *
 * The equivalent Fortran interfaces are
 *
 *     integer function IVEQDQueryMemberBool(qdm, name, bvalue)
 *     integer function IVEQDQueryMemberInt(qdm, name, ivalue)
 *     integer function IVEQDQueryMemberStringLength(qdm, name, ivalue)
 *     integer function IVEQDQueryMemberString(qdm, name, cvalue)
 *     integer function IVEQDQueryMemberSListLength(qdm, name, ivalue)
 *     integer function IVEQDQueryMemberSListMemberLength(qdm, name, ientry, ivalue)
 *     integer function IVEQDQueryMemberSListMember(qdm, name, ientry, cvalue)
 *     integer qdm, ientry, ivalue
 *     logical bvalue
 *     character*(*) name, cvalue
 *
 * As of this writing, these are the properties that are defined:
 *
 *     Name                     Type  Meaning
 *     AcceptsGPU               B     Indicates whether the queue accepts
 *                                    jobs making use of GPUs as coprocessors.
 *     AcceptsSerial            B     Indicates whether the queue accepts
 *                                    serial (single-processor) jobs.
 *     AcceptsThreaded          B     Indicates that the queue accepts parallel
 *                                    jobs that require something that is (or
 *                                    acts like) a single host.
 *     AcceptsCluster           B     Indicates whether the queue accepts
 *                                    parallel jobs that may run on multiple
 *                                    hosts.
 *     DefaultBackgroundLoginOption
 *                              S     Holds the option to pass to the command
 *                                    given in the DefaultLoginCommand when the
 *                                    task is to be run in the background.
 *     DefaultLoginCommand      S     Is the default command to use to log in
 *                                    and run a task on another compute node.
 *                                    Generally assume the command has a
 *                                    similar interface to rsh or ssh.
 *                                    DefaultLoginCommand should be set for
 *                                    systems or queues that set ProcessControl
 *                                    to "mpich-machines" or make use of
 *                                    GrecsrvHosts.
 *     DefaultMachinesFile      S     Holds the name of the MPICH-style
 *                                    machines file containing listing the
 *                                    hosts to use.  DefaultMachinesFile should
 *                                    be set for systems or queue that set
 *                                    ProcessControl to "mpich-machines".
 *     DefaultMPIRunCommand     S     Holds the command necessary for starting
 *                                    an MPI job.  For systems or queues that
 *                                    set ProcessControl to "openmpi", the
 *                                    command must be OpenMPI's (version 1.3 or
 *                                    later) mpirun or mpiexec
 *     DefaultNodeCheckCommand  S     Is the name of the command to run
 *                                    immediately prior to starting a process
 *                                    when ProcessControl has been set to
 *                                    "mpich-machines" or "openmpi".
 *     DefaultProcessorCount    I     Hold the default number of processors
 *                                    to use during parallel processing.
 *     Description              S     Holds a longer description of the queue.
 *     GrecsrvBusyGPUs          S     Holds a comma-separated list of GPU ids
 *                                    to ignore on each host.
 *     GrecsrvExecutable        S     Holds the path to the GPU reconstruction
 *                                    server
 *     GrecsrvHosts             L     Lists the host name or addresses for
 *                                    hosts to use for GPU reconstruction
 *                                    servers.
 *     GrecsrvMemoryPercentage  I     Holds the percentage of each device's
 *                                    onboard memory that the GPU
 *                                    reconstruction server will attempt to
 *                                    use.
 *     GrecsrvPort              I     Is the port number to use when
 *                                    communicating to the GPU reconstruction
 *                                    server when the host address does not
 *                                    specify an address.
 *     GrecsrvRunning           B     If true, the GPU reconstruction server
 *                                    is assumed to be runnning.  Otherwise,
 *                                    it will be started and stopped on demand.
 *     GrecsrvUseProcessControl B     If true, use the setting from
 *                                    ProcessControl for choosing hosts and
 *                                    starting GPU reconstruction servers.  If
 *                                    false, GrecsrvHosts lists the default
 *                                    hosts to use for GPU reconstruction
 *                                    servers.
 *     Name                     S     Holds the name of the queue appropriate
 *                                    for use in a user interface.
 *     ProcessControl           S     Defines what mechanism should be used by
 *                                    parallel jobs that start their own
 *                                    child processes to determine where the
 *                                    processes should run and how the tasks
 *                                    should be started.  Three values are
 *                                    currently understood:  "simple",
 *                                    "mpich-machines", and "openmpi".
 *                                    "simple" means that all the processes
 *                                    should be started on the same host as
 *                                    the parent process.  "mpich-machines"
 *                                    means that the job will be passed a text
 *                                    file, in the same format as used by MPICH
 *                                    for its machines files, listing the hosts
 *                                    to use.  The job will also be passed
 *                                    the command to use to log in to other
 *                                    hosts.  The AutoGeneratedMachinesFile,
 *                                    DefaultBackgroundLoginOption,
 *                                    DefaultLoginCommand, DefaultMachinesFile,
 *                                    and DefaultNodeCheckCommand properties
 *                                    should be set for systems or queues that
 *                                    use "mpich-machines" for ProcessControl.
 *                                    "openmpi" means that the job will be
 *                                    passed how to invoke OpenMPI's mpirun or
 *                                    mpiexec and should use that when starting
 *                                    processes.  The DefaultMPIRunCommand
 *                                    should be set for systems or queues that
 *                                    use "openmpi" for ProcessControl.
 *     SetProcessorCountForApp  B     If true, the requested number of
 *                                    processors is passed to an application.
 *                                    If false, the requested number of
 *                                    processors is only passed to the command
 *                                    given by SubmissionCommand.  When not
 *                                    set, this property defaults to true.
 *     SubmissionCommand        S     Holds the command necessary for
 *                                    submitting a job to the queue.  Clients
 *                                    assume that the value of
 *                                    SubmissionCommand + whitespace +
 *                                    executable will work to start a job.
 *                                    %p in the command is to be expanded
 *                                    to the number of processors to use.
 *                                    %% in the command is to be expanded to %.
 *     StatusCommand            S     Holds the command necessary for checking
 *                                    the status of a queue.  If empty, the
 *                                    status of the queue can not be queried.
 */
int IVEQDQueryMemberBool(IVEQDMember qdm, const char* name, int* p_value);
int IVEQDQueryMemberInt(IVEQDMember qdm, const char* name, int* p_value);
int IVEQDQueryMemberStringLength(
    IVEQDMember qdm, const char* name, int* p_value
);
int IVEQDQueryMemberString(
    IVEQDMember qdm, const char* name, const char** p_value
);
int IVEQDQueryMemberSListLength(
    IVEQDMember qdm, const char* name, int* p_value
);
int IVEQDQueryMemberSListMemberLength(
    IVEQDMember qdm, const char* name, int ientry, int* p_value
);
int IVEQDQueryMemberSListMember(
    IVEQDMember qdm, const char* name, int ientry, const char** p_value
);

/*
 * If dest is non-null and has enough space to hold the host name (as
 * indicated by ndest, the number of characters that dest can hold) fills
 * dest with the host name.  Returns the length of the host name (not
 * including the terminating null) or -1 if unable to determine the length.
 *
 * The corresponding Fortran interface has the prototype (the ndest argument
 * is implict).
 *
 *     integer function IVEQDGetHostname(dest)
 *     character*(*) dest
 */
int IVEQDGetHostname(int ndest, char* dest);

/*
 * IVEQDConstructExpansionContext sets up a context for expanding a command
 * (namely the SubmissionCommand).  Returns zero if not successful.
 * IVEQDDesroyExpansionContext frees any resources associated with a context.
 * IVEQDExpansionContextSetIntProperty and
 * IVEQDExpansionContextSetStringProperty modify the contents of the context.
 * Those functions return zero if successful and a nonzero value if not.
 * Currently the following properties have meanings:
 *
 *     ProcessorCount           I     The number of processors to use.
 *                                    Defaults to one.
 *
 * IVEQDExpandCommand expands the command.  It puts the result in dest
 * unless dest is null or ndest, the maximum number of characters that
 * dest can hold, is smaller than necessary for the expanded command.
 * IVEQDExpand returns the full length (not including the terminating null)
 * of the expanded command or -1 if it was not possible to determine the
 * length.
 *
 * The corresponding Fortran interfaces are (for IVEQDExpandCommand(),
 * ndest is implicit):
 *
 *     integer function IVEQDConstructExpansionContext()
 *     subroutine IVEQDDestroyExpansionContext(icont)
 *     integer function IVEQDExpansionContextSetIntProperty(icont, name, i)
 *     integer function IVEQDExpansionContextSetStringProperty(icont, name, s)
 *     integer function IVEQDExpandCommand(cmd, icont, dest)
 *     integer icont, i
 *     character(*) name, cmd, s, dest
 */
IVEQDExpansionContext IVEQDConstructExpansionContext(void);
void IVEQDDestroyExpansionContext(IVEQDExpansionContext context);
int IVEQDExpansionContextSetIntProperty(
    IVEQDExpansionContext context, const char* name, int value
);
int IVEQDExpansionContextSetStringProperty(
    IVEQDExpansionContext context, const char* name, const char* value
);
int IVEQDExpandCommand(
    const char* cmd, IVEQDExpansionContext context, int ndest, char* dest
);

#endif /* include guard */
