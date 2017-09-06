// This file provides assertion macro commands similar to ones in boost/assert.hpp
// It intentionally does not have include guards

#undef PICA_ASSERT
#undef PICA_ASSERT_MSG

#if defined(PICA_DISABLE_ASSERTS) || ( defined(PICA_ENABLE_ASSERT_DEBUG_HANDLER) && defined(NDEBUG) )

#define PICA_ASSERT(expr) ((void)0)
#define PICA_ASSERT_MSG(expr, msg) ((void)0)

#elif defined(PICA_ENABLE_ASSERT_HANDLER) || ( defined(PICA_ENABLE_ASSERT_DEBUG_HANDLER) && !defined(NDEBUG) )

// In this case a user needs to implement the following functions in namespace pica:
///namespace pica
///{
///void assertionFailed(const char* expr, const char* file, long line);
///void assertionFailedMsg(const char* expr, const char* msg, const char* file, long line);
///}

#define PICA_ASSERT(expr) ((expr)? ((void)0): ::pica::assertionFailed(#expr, __FILE__, __LINE__))
#define PICA_ASSERT_MSG(expr, msg) ((expr)? ((void)0): ::pica::assertionFailedMsg(#expr, msg, __FILE__, __LINE__))

#else

#include <cassert>

#define PICA_ASSERT(expr) assert(expr)
#define PICA_ASSERT_MSG(expr, msg) assert((expr)&&(msg))

#endif
