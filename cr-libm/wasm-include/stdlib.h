/* yFoil: the C library headers CORE-MATH includes, for wasm32 (see build.rs and fenv.h).
 * exit is called only after an "unexpected worst case" report (stdio.h); it traps, which
 * aborts the call as the native exit(1) aborts the process. */
#ifndef YFOIL_WASM_STDLIB_H
#define YFOIL_WASM_STDLIB_H

static inline void exit(int status) { (void)status; __builtin_trap(); }

#endif
