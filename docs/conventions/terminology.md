# Terminology

Words that carry a precise meaning in yFoil's validation, tests and studies. Use them exactly
as defined here; when a document, a test or a study needs one of these ideas, it uses this word
and no other.

## Deterministic, not nondeterministic

XFOIL is **deterministic**: the same input bits, on the same host, give the same output bits,
every time. So is yFoil. Never describe either as *nondeterministic*, *random* or *noisy*.

Two qualifications, neither of which is nondeterminism:

- **Host dependence.** Transcendentals (`**`, `EXP`, `LOG`, `ATAN2`, `TANH`) come from the host
  libm in both codes. Apple libSystem and glibc differ by one ULP on a fraction of inputs, so
  the same input can give different bits on different hosts. Deterministic per host; an ULP
  budget across hosts (CLAUDE.md Rule 1).
- **Ill-conditioned solutions** and **divergent comparisons**, below.

## Ill-conditioned solution

A **solution is ill-conditioned** where a tiny change of input moves the output a lot: near
stall, in deep separation, near the Kármán–Tsien pole, at a Mach number pushed towards 0.99. The
numerical method can no longer follow the physics there. It is a property of the problem at that
point, met by **each code on its own**: XFOIL disagrees with itself under a one-ULP perturbation
of its input, and so does yFoil.

Ill-conditioned solutions are not to be avoided in testing. Most of XFOIL's rarely-taken branches
— the garbage extrapolation of a failed station, the relaxation limits, the negative-Ue islands,
the stagnation-point fallbacks — exist precisely to cope with them, so the tests must reach them
to cover those branches.

A result can be ill-conditioned and still valid: conditioning is Class C in
[solution validity](../guide/validity.md), which a single run cannot detect.

## Divergent comparison

A **comparison is divergent** when XFOIL and yFoil, run on the same input, take different routes
through the code: at some real comparison `IF (A .GT. B)` the operands are equal to within
rounding, and the last-bit difference between the two codes (different association, a host libm
ULP) falls on opposite sides. XFOIL has about 140 real inequality branches; the ones met in
practice include `AMPL > ACRIT` (the transition station), `RMSBL < EPS1` (convergence) and the
station at which UPDATE's relaxation factor is limited.

A divergent comparison is a *consequence* of an ill-conditioned solution — the control flow of a
well-conditioned solution does not sit on a knife edge — but it is not itself ill-conditioning.
It says nothing about whether the translation is right, so **no test is ever a divergent
comparison** (`testing.md`). Whether a comparison diverges is decided by observing both codes'
routes through it, not by the size of any difference in the results.

In validation, as opposed to testing, divergent comparisons are expected in a polar that contains
ill-conditioned solutions: transition may move one panel at α = 10° in XFOIL and at α = 11° in
yFoil, giving much the same polar by a different path. The twins and sensitivity studies explain
such differences; they are not translation errors.

This replaces the earlier term *threshold-straddling*.

## Twins

A **twin** is a rerun of the reference with every panel coordinate jogged by −1, 0 or +1 ULP,
x and y independently, drawn reproducibly from a seed. Five seeded twins per case are the
instrument that measures how ill-conditioned a solution is: where a twin's branch trace differs
from the reference's, the reference's own route sits on a knife edge there. They are made by
`cargo xtask fixtures` for every case with `twins = true` — every validation run — or by
`cargo xtask twins` (`target/fixtures/<case>/ulp<seed>/`). A twin describes the run it is the
twin of: its conditioning says nothing about the same section at another panelling or reached
along another path.

Moving every coordinate the same way is not a twin: it translates the aerofoil and measures
almost nothing.

## Noise floor

The **noise floor** of a recorded value is the largest difference between the reference and any
of its twins. It is the reference's own spread under a one-ULP perturbation. It is a study
quantity, used to derive the named tolerances in `tests/common/utilities/tolerances.rs` and in
the studies, which read it from the `noise_floor.json` written with the twins beside each
case's work directory. Tests do not read it (see
[testing.md](testing.md), including its Status section).

## Case cover and step cover

A **cover** is a set cover: a collection of runs that together take every branch of XFOIL's
translated subroutines that some candidate takes. The project has two, and they are not
interchangeable.

The **case cover** is the fewest whole reference cases that together take every branch of the
analysis path a finite run can reach. It is chosen by the branch-coverage study
(`scripts/branch-coverage`) from each candidate case's gcov counts, and its members are the cases
with `group = "branch-coverage"` in `xtask/fixtures-config/cases.toml`. A **case-cover member** is
one of those cases (for example `naca0012_n60_a8_re1e4_damp`). It is a study's claim that this set
of runs exercises the analysis path; the tests do not gate branches through it.

The **step cover** is the fewest single steps — one Newton iteration of one VISCAL call, or one
inviscid call — that together take every branch any candidate step takes. It is chosen by
`cargo xtask steps` from every tracked case's steps, only among steps through which
`cargo xtask route` observed yFoil take XFOIL's route, and is recorded in
`xtask/fixtures-config/step-cover.toml`. A **step-cover member** is one such step, replayed by one
test in `tests/execution/branches.rs`. This is the cover the tests gate branches with
([testing.md](testing.md), "How branch coverage is established").

Say which cover is meant; "the cover" alone is ambiguous.
