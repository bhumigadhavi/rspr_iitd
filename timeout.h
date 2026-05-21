#ifndef TIMEOUT_H
#define TIMEOUT_H

#include <csignal>

// This tells the compiler "this variable exists somewhere, but not here."
extern volatile sig_atomic_t timeout_signaled;

// Declare the signal handler function prototype
void signal_handler(int signal);

#endif
