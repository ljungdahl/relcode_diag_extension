! // This file contain C preprocessor macros for
! // controlling the write-logging to stdout.
! // It needs to be included as a header for use.
#pragma once

#ifndef USE_LOG_FATAL
#define USE_LOG_FATAL 1
#endif

#ifndef USE_LOG_WRITE
#define USE_LOG_WRITE 1
#endif

#if USE_LOG_WRITE == 1
#define LOG_WRITE write(*,*)
#else
#define LOG_WRITE
#endif

#if USE_LOG_FATAL == 1
#define LOG_FATAL(msg) write(*,'(a,a,a,a,a,i0)')"FATAL: ",msg,NEW_LINE('a'),__FILE__,", line:",__LINE__; \
    call exit(1)
#else
#define LOG_FATAL(msg)
#endif