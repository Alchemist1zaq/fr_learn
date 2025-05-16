#pragma once

#include <cstdlib>
#include <execinfo.h>
#include <iostream>
#include <stdio.h>
#include <unistd.h>

//! Prints the error message, the source file and line number, and exits
#define fatalError(s)                                                          \
  {                                                                            \
    printf("Fatal error '%s' at %s:%d\n", s, __FILE__, __LINE__);              \
    exit(1);                                                                   \
  }

//! Prints the error message, the source file ane line number, the full stack
//! trace, and exits
#define fatalErrorST(s)                                                        \
  {                                                                            \
    void *array[10];                                                           \
    size_t size;                                                               \
    size = backtrace(array, 10);                                               \
    printf("Fatal error '%s' at %s:%d\n\n", s, __FILE__, __LINE__);            \
    backtrace_symbols_fd(array, size, STDERR_FILENO);                          \
    exit(1);                                                                   \
  }

#define _(x) cout << #x << ": " << x << endl;
#define _print(x, y) cout << #x << ": " << x << ", " << #y << ": " << y << endl;
