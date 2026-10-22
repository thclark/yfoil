# Instrumentation patches and logging

The instrumented build is the pristine reference plus a series of patches that only write: each
dumps some of XFOIL's state, at full precision, at a point in its execution where a test needs it.
The patches are `xfoil/instrumentation/instrument/*.patch`, applied in the order of
`xfoil/instrumentation/series.instrument` to a staged copy of the pristine source
([how it is built](index.md#layout-and-how-it-is-built)); that they change nothing is
[proved on every build](proofs-that-the-reference-is-sound.md#proof-that-instrumentation-changes-nothing).
Two further patches are inert unless a case asks for them: the branch events that
[validity flags](../execution-equivalence/validity-flags-against-xfoils-events.md) are checked
against, and the stops on a non-finite state or after N SETBL calls that the
[fixtures](../fixtures.md#case-definitions-and-groups) of the long studies use.

## Where output goes, and in what format

Every dump is written to the **working directory** — never `/tmp` — so a run is self-contained and
the fixture pipeline can keep a whole case in one directory. Rather than list them here, where the
list would rot: the authoritative set is the `FILE=` targets in
`xfoil/instrumentation/instrument/*.patch`, and any tracked case directory under
`tests/fixtures/xfoil/` shows what a run of that shape actually produces. They fall into five groups:

- **per operating point** — `specal_points.dat`, `speccl_points.dat`, `viscal_points.dat`
- **per VISCAL iteration** — `viscal_iter.dat`, `viscal_iters_all.dat`, `viscal_final.dat`,
  `viscal_inviscid.dat`, `viscal_nonfinite.dat`
- **per SETBL call** (written for the calls listed in `dump_calls.txt`, default 1–3; a case's
  `dump_calls`, `step_calls` and the branch cover fill it) — `mrchdu_input_<k>.dat`,
  `mrchdu_output_<k>.dat`, `setbl_output_<k>.dat`, `update_output_<k>.dat`, `mrchdu_trace_<k>.dat`
  and `trchek2_<k>.dat`; a case tracks only those its tests read
- **one-shot state** — `xfoil_dij.dat`, `xfoil_pointers.dat`, `xfoil_uinv.dat`, `xfoil_panels.dat`,
  `xfoil_inviscid.dat`, `xfoil_pangen.dat`, `xfoil_tgap.dat`, the BLSOLV inputs and traces
- **branch events and call traces** — `events.dat` ([validity flags against XFOIL's events](../execution-equivalence/validity-flags-against-xfoils-events.md)),
  `xfoil_trace.log`, `xfoil_subroutine_log.dat`, `dump_calls.txt`

**All of it at `ES24.16`** — 17 significant figures, the minimum that round-trips an f64. XFOIL's own
debug formats (`E24.16`, `E18.10`, `E12.4`) carry 16, 10 and 4 figures respectively and were the
reason earlier fixtures were not bit-exact; `E24.16` is *not* enough for a double. Prefer the
`xlog.f` JSON harness (`LOGENTER`/`LOGINT`/`LOGREAL`/`LOGBOOL`/`LOGARR1`/`LOGARR2`) for new work;
the per-file `.dat` dumps are legacy and are being migrated stage by stage.

No test compares against XFOIL's *formatted* output files — `PSAVE` is `G15.7`, `.pol` is `F9.4/F10.5`,
`DUMP` is `F9.5/F10.6`, `CPWR` is `F11.5`, none of which can support a comparison below ~1e-5.

## Adding new instrumentation

Edit a staged tree and regenerate the patch; do not edit `target/xfoil-ref/` in place, and do not
touch the pristine source.

```bash
# stage a tree with the full series applied
W=.tmp/instr; rm -rf $W; mkdir -p $W
cp -R xfoil/third-party/xfoil-6.99 $W/tree
for s in series.build series.instrument; do
  while IFS= read -r p; do [ -z "$p" ] && continue
    patch -s -p1 -d $W/tree < xfoil/instrumentation/$p
  done < xfoil/instrumentation/$s
done
cp -R $W/tree $W/base          # snapshot to diff against

# ... edit $W/tree/src/<file>.f ...

diff -u $W/base/src/xfoil.f $W/tree/src/xfoil.f \
  | sed '1s|.*|--- a/src/xfoil.f|;2s|.*|+++ b/src/xfoil.f|' \
  > xfoil/instrumentation/instrument/NN-my-change.patch
echo "instrument/NN-my-change.patch" >> xfoil/instrumentation/series.instrument
scripts/xfoil-build.sh --verify         # must still prove the instrumentation inert
```

Three rules for the edit itself:

- **Write only.** Instrumentation may not change control flow or any stored value. `--verify` is the
  proof ([proofs that the reference is sound](proofs-that-the-reference-is-sound.md)), and it must still pass.
- **Declare your locals.** XFOIL has no `IMPLICIT NONE`, so a name beginning `A-H` or `O-Z` is `REAL`
  and `I-N` is `INTEGER` by default. Declare anything you add explicitly.
- **Branch events go to `events.dat`** via `EVLOG(NAME, IS, IBL, R)`, not to a new file. That is how
  the rarely-taken fallbacks are recorded, and `tests/execution/events.rs` asserts them
  one call at a time.

Cases are then driven by `cargo xtask fixtures`, which generates panels with yFoil, runs the
instrumented build via `LOAD`, asserts the bitwise geometry handoff and keeps the dumps with a
`manifest.json` ([fixtures](../fixtures.md)).
