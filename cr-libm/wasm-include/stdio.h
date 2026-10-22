/* yFoil: the C library headers CORE-MATH includes, for wasm32 (see build.rs and fenv.h).
 * printf is called by debug printers nothing calls and before the exit(1) of the "unexpected
 * worst case" reports in atan2 and pow (cases CORE-MATH's proofs exclude). The browser has no
 * standard output, so it prints nothing; the exit that follows traps (stdlib.h). A macro, not a
 * function, so that the compiler cannot turn the call into a puts the module would import. */
#ifndef YFOIL_WASM_STDIO_H
#define YFOIL_WASM_STDIO_H

#define printf(...) 0

#endif
