/* Linked into the gcov build of the reference only (scripts/xfoil-build.sh --gcov).
 *
 * A gcov-instrumented program writes its branch counters (*.gcda) at normal exit; a run the
 * fixture watchdog has to end — XFOIL hung in its plot label on a non-finite value, or still
 * marching after the time limit — would otherwise lose every count it made. On SIGTERM this
 * writes the counters and exits, so the watchdog's TERM keeps them. The solver is untouched.
 */
#include <signal.h>
#include <unistd.h>

extern void __gcov_dump(void);

static void flush_and_exit(int sig) {
    (void)sig;
    __gcov_dump();
    _exit(0);
}

__attribute__((constructor)) static void install_flush(void) {
    signal(SIGTERM, flush_and_exit);
}
