! // This file contain C preprocessor macros for
! // controlling the write-logging to stdout.
! // It needs to be included as a header for use.
#pragma once

#ifndef USE_LOG_WRITE
#define USE_LOG_WRITE 1
#endif

#if USE_LOG_WRITE == 1
#define LOG_WRITE write(*,*)
#else
#define LOG_WRITE
#endif