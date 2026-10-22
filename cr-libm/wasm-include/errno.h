/* yFoil: the C library headers CORE-MATH includes, for wasm32 (see build.rs and fenv.h).
 * errno is only written under CORE_MATH_SUPPORT_ERRNO, which the build does not define; this
 * header exists so the sources that include it unconditionally compile. Each translation unit
 * gets its own write-only errno, so no symbol is exported or imported. */
#ifndef YFOIL_WASM_ERRNO_H
#define YFOIL_WASM_ERRNO_H

#define EDOM   33
#define ERANGE 34

static int yfoil_wasm_errno __attribute__((unused));
#define errno yfoil_wasm_errno

#endif
