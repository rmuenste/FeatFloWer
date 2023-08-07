#include <unistd.h>

int get_C_pid() {
    int pid = getpid();
//     printf("PID from C: %d\n", pid);  // Print the PID value
    return pid;
}