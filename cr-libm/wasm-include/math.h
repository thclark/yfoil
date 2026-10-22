/* yFoil: the C library headers CORE-MATH includes, for wasm32 (see build.rs and fenv.h).
 * CORE-MATH takes only NAN and INFINITY from <math.h>; every operation it uses is a
 * __builtin_ (fma, floor, round, ...), which clang lowers to a wasm instruction or, for fma
 * and round, a call resolved at link time from Rust's compiler-builtins (wasm32-unknown-unknown)
 * or wasi-libc (wasm32-wasip1), both ports of musl's correctly rounded implementations. */
#ifndef YFOIL_WASM_MATH_H
#define YFOIL_WASM_MATH_H

#define NAN      __builtin_nanf("")
#define INFINITY __builtin_inff()

#endif
