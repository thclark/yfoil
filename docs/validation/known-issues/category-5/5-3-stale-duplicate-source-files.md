# 5.3 Stale duplicate source files — Out of scope

`dplot1.f` is an older snapshot of `dplot.f` (missing `FBLGET`); `xqdes1.f`, `modify1/2/3.f`,
`frplot0.f` and `xpanel.new` are the same pattern. None is in the build (`bin/Makefile` links
`dplot.o` only) and none is in `coverage.toml`'s file list.
