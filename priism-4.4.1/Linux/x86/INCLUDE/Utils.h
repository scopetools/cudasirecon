#ifndef UCSF_MSG_UTILS_H
#define UCSF_MSG_UTILS_H

int ParseComLine(
    int argc, char** argv, int* wid_ptr, int* bank_ptr, int* sec_ptr
);

int SplitName(char* string1, char* string2, char* string3);
int MergeName(char* dest_buffer, char* file_buffer, char* dir_buffer);


/*
 * Attempts to execute a program whose full path name is given by func.
 * res is passed as the first command-line argument and env as the second
 * command-line argument.
*/
void IVEExecFunc(char* func, char* res, char* env);

/*
 * Like IVEExecFunc, but name is interpreted relative to the IVE executable
 * directory.
 */
void IVEExeAppl(char* name, char* res, char* env);


/*
 * Executes the string
 *
 *   <path to IVE executable directory>/<name>  <cmd>
 *
 * as a shell command and returns when the command completes or is stopped.
 */
void IVESysAppl(char* name, char* cmd);

#endif /* include guard */
