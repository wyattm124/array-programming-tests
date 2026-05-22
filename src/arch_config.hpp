#pragma once

// Target architecture detection.
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
#define FFT_ARCH_X86 1
#else
#define FFT_ARCH_X86 0
#endif

#if defined(__aarch64__) || defined(__arm__) || defined(_M_ARM) || defined(_M_ARM64)
#define FFT_ARCH_ARM 1
#else
#define FFT_ARCH_ARM 0
#endif

// ISA / SIMD feature detection for the current translation unit.
#if FFT_ARCH_X86 && defined(__AVX2__)
#define FFT_HAS_AVX2 1
#else
#define FFT_HAS_AVX2 0
#endif

#if FFT_ARCH_ARM && (defined(__ARM_NEON) || defined(__ARM_NEON__))
#define FFT_HAS_NEON 1
#else
#define FFT_HAS_NEON 0
#endif
