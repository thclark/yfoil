/* yFoil: the C library headers CORE-MATH includes, for wasm32 (see build.rs and fenv.h).
 * Only the printf format macros of its debug printers are used; int64_t is long long on wasm32. */
#ifndef YFOIL_WASM_INTTYPES_H
#define YFOIL_WASM_INTTYPES_H

#include <stdint.h>

#define PRId64 "lld"
#define PRIx64 "llx"

#endif
