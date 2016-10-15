#ifndef __SSE_HXX__
#define __SSE_HXX__

#if defined(__GNUC__)
#define MM_ALIGN16             __attribute__ ((aligned (16)))
#define MM_ALIGN32             __attribute__ ((aligned (32)))
#define MM_ALIGN64             __attribute__ ((aligned (64)))
#elif defined(__MSVC__)
#define MM_ALIGN16             __declspec(align(16))
#define MM_ALIGN32             __declspec(align(32))
#define MM_ALIGN64             __declspec(align(64))
#else
#error Unsupported compiler!
#endif

#define STACK_ALIGN16(t)       reinterpret_cast<float *>((reinterpret_cast<size_t>(t)+0x0F) & ~(size_t) 0x0F)
#define STACK_ALIGN32(t)       reinterpret_cast<float *>((reinterpret_cast<size_t>(t)+0x1F) & ~(size_t) 0x1F)
#define STACK_ALIGN64(t)       reinterpret_cast<float *>((reinterpret_cast<size_t>(t)+0x3F) & ~(size_t) 0x3F)

/* ========= SSE intrinsics ========= */

/* Include SSE intrinsics header file */
#include <xmmintrin.h>
#include <emmintrin.h>


//#ifdef __SSE3__
#if 1
#include <pmmintrin.h>
#endif

//#ifdef __SSSE3__
#if 1
#include <tmmintrin.h>
#endif




/* MSVC intrinsics header (for RDTSC) */
#if defined(__MSVC__)
# include <intrin.h>
# pragma intrinsic(__rdtsc)
#endif

#define splat_ps(ps, i)          _mm_shuffle_ps   ((ps),(ps), (i<<6) | (i<<4) | (i<<2) | i)
#define splat_epi32(ps, i)       _mm_shuffle_epi32((ps), (i<<6) | (i<<4) | (i<<2) | i)
#define mux_ps(sel, op1, op2)    _mm_or_ps   (_mm_and_ps   ((sel), (op1)), _mm_andnot_ps   ((sel), (op2)))
#define mux_epi32(sel, op1, op2) _mm_or_si128(_mm_and_si128((sel), (op1)), _mm_andnot_si128((sel), (op2)))
#define enable_fpexcept_sse()	 _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO))
#define query_fpexcept_sse()	 (~_MM_GET_EXCEPTION_MASK() & (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO))
#define disable_fpexcept_sse()	 _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO))
#define load1_epi32(i)           _mm_shuffle_epi32(_mm_cvtsi32_si128(i), 0)
#define negate_ps(val)           _mm_xor_ps((val), SSEConstants::negation_mask.ps)

#define pstoepi32(ps)            _mm_castps_si128(ps)
#define epi32tops(pi)            _mm_castsi128_ps(pi)

  
#endif //__SSE_HXX__
