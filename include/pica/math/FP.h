#ifndef PICA_FP_H
#define PICA_FP_H


namespace pica {

// Switch to use single or double precision as type FP, double by default.
#ifndef PICA_USE_SINGLE_PRECISION
    typedef double FP;
#else
    typedef float FP;
#endif

inline FP sqr(FP x) { return x * x; }

} // namespace pica


#endif
