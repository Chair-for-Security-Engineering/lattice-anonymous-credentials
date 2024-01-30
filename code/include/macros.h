#ifndef MACROS_H
#define MACROS_H

#include <assert.h>
#include <execinfo.h>

#ifdef REMOVE_ASSERTS
#define ASSERT_DEBUG(test, msg) ;
#else
#define ASSERT_DEBUG(test, msg) do {\
        if (!(test)) {\
            \
            void* callstack[128];\
            int frames = backtrace(callstack, 128);\
            char** strs = backtrace_symbols(callstack, frames);\
            for (int i = 0; i < frames; ++i) {\
                printf("%s\n", strs[i]);\
            }\
            free(strs);\
            printf("%s\n", msg);\
            fflush(stdout);\
            assert(false);\
        }\
    } while (0);
#endif

#define CHK_UI_OVF_ADDITION(sum, summand) do {uint64_t __watchguard = sum; sum += summand; assert(__watchguard <= sum);} while (0);

#endif /* MACROS_H */
