/* yFoil: the C library headers CORE-MATH includes, for wasm32 (see build.rs).
 *
 * WebAssembly has one rounding mode, round-to-nearest-even, and no floating-point exception
 * flags (WebAssembly core spec, "Floating-Point Operations"). This header says exactly that:
 * fegetround() is FE_TONEAREST, the flags are never raised and always read as clear.
 *
 * CORE-MATH needs it because the sources name every rounding mode and every flag (wasi-libc's
 * fenv.h defines only FE_TONEAREST and FE_ALL_EXCEPT, so they do not compile against it). No
 * returned value depends on it: fegetround() is FE_TONEAREST natively too (nothing in yFoil
 * changes the rounding mode), and the flag functions only set, save and restore the flags
 * (atan2's fetestexcept results decide errno and which flags to clear, never the result).
 * The macro values are arbitrary distinct bits; only FE_TONEAREST is ever returned.
 */
#ifndef YFOIL_WASM_FENV_H
#define YFOIL_WASM_FENV_H

#define FE_INVALID    0x01
#define FE_DIVBYZERO  0x02
#define FE_OVERFLOW   0x04
#define FE_UNDERFLOW  0x08
#define FE_INEXACT    0x10
#define FE_ALL_EXCEPT 0x1f

#define FE_TONEAREST  0
#define FE_DOWNWARD   1
#define FE_UPWARD     2
#define FE_TOWARDZERO 3

typedef unsigned int fexcept_t;
typedef struct { unsigned int flags; } fenv_t;

static inline int fegetround(void) { return FE_TONEAREST; }
static inline int fetestexcept(int excepts) { (void)excepts; return 0; }
static inline int feclearexcept(int excepts) { (void)excepts; return 0; }
static inline int feraiseexcept(int excepts) { (void)excepts; return 0; }
static inline int fegetexceptflag(fexcept_t *flag, int excepts) { (void)excepts; *flag = 0; return 0; }
static inline int fesetexceptflag(const fexcept_t *flag, int excepts) { (void)flag; (void)excepts; return 0; }
static inline int feholdexcept(fenv_t *env) { env->flags = 0; return 0; }
static inline int feupdateenv(const fenv_t *env) { (void)env; return 0; }

#endif
