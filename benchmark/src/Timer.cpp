#include "Timer.h"

#ifdef _WIN32

#include <windows.h>


class QPCStopwatch : public Stopwatch
{
    LARGE_INTEGER last_start;
    long long total;

public:
    QPCStopwatch()
    { reset(); }

    void reset()
    { total = 0; }

    void start()
    { QueryPerformanceCounter(&last_start); }

    void stop()
    {
        LARGE_INTEGER now;
        QueryPerformanceCounter(&now);
        total += now.QuadPart - last_start.QuadPart;
    }

    double getElapsed()
    {
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        return total / (double) freq.QuadPart;
    }
};


#else

#include <time.h>

#ifdef CLOCK_MONOTONIC


class CGTStopwatch : public Stopwatch
{
    timespec last_start;
    long long total;
    
public:
    CGTStopwatch()
    { reset(); }
    
    void reset()
    { total = 0; }
    
    void start()
    { clock_gettime(CLOCK_MONOTONIC, &last_start); }
    
    void stop()
    {
        timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        total += (now.tv_sec - last_start.tv_sec) * 1000000000LL
            + now.tv_nsec - last_start.tv_nsec;
    }
    
    double getElapsed()
    {
        return total / 1e9;
    }
};

#else

#include <sys/time.h>


/* This is a backup timer for old Unix systems.
   It uses the now-obsolescent gettimeofday function,
   which should be available everywhere. */
class GTODStopwatch : public Stopwatch
{
    timeval last_start;
    long long total;
    
public:
    GTODStopwatch()
    { reset(); }
    
    void reset()
    { total = 0; }
    
    void start()
    { gettimeofday(&last_start, NULL); }
    
    void stop()
    {
        timeval now;
        gettimeofday(&now, NULL);
        total += (now.tv_sec - last_start.tv_sec) * 1000000LL
            + now.tv_usec - last_start.tv_usec;
    }
    
    double getElapsed()
    {
        return total / 1e6;
    }
};

#endif

#endif


Stopwatch * createStopwatch()
{
#if defined(_WIN32)
    return new QPCStopwatch();
#elif defined(CLOCK_MONOTONIC)
    return new CGTStopwatch();
#else
    return new GTODStopwatch();
#endif
}
