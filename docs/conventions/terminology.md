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
- **Ill-conditioning**, below.

## Ill-conditioned

A **branch is ill-conditioned** when the two operands of the real comparison that decides it are
equal to within rounding, so that the decision depends on the last bit of the input. A one-ULP
change to the input panels flips the decision, and the run that follows differs by O(1).
XFOIL has about 140 real inequality branches; the ones met in practice include

- `AMPL > ACRIT`, which picks the transition station (`ITRAN`);
- `RMSBL < EPS1`, which decides whether a VISCAL call has converged;
- the station at which UPDATE's relaxation factor `RLX` is limited.

An **ill-conditioned iteration** or **ill-conditioned result** is one in which yFoil follows a
different code path from the reference because a rounding-level difference fell on such a
branch — not because of a translation error. A step or a point that is not ill-conditioned is
**well-conditioned**.

Ill-conditioning is a property of XFOIL's problem at that point, not of either code: XFOIL
disagrees with *itself* there under a one-ULP perturbation. It is invisible in a single run and
is found by perturbing the input (the twins, below).

Consequences:

- An XFOIL–yFoil difference after an ill-conditioned branch is not, by itself, evidence of a bug.
  Larger differences are expected in a polar that contains ill-conditioned solutions: transition
  may move one panel at α = 10° in XFOIL and at α = 11° in yFoil, giving much the same polar by a
  different path.
- Tests never gate an ill-conditioned step (see [testing.md](testing.md)); the studies report
  where they occur.
- A result can be ill-conditioned and still valid: conditioning is Class C in
  [solution validity](../guide/validity.md), which a single run cannot detect.

This replaces the earlier term *threshold-straddling*.

## Twins

A **twin** is a rerun of the reference with every panel coordinate jogged by −1, 0 or +1 ULP,
x and y independently, drawn reproducibly from a seed. Five seeded twins per case are the
instrument that finds ill-conditioned branches: where a twin's branch trace differs from the
reference's, the reference is ill-conditioned there.

Moving every coordinate the same way is not a twin: it translates the aerofoil and measures
almost nothing.

## Noise floor

The **noise floor** of a recorded value is the largest difference between the reference and any
of its twins. It is the reference's own spread under a one-ULP perturbation. It is a study
quantity, used to choose which steps are well-conditioned enough to become test cases and to
derive the named tolerances in `tests/utilities/tolerances.rs`. Tests do not read it (see
[testing.md](testing.md), including its Status section).
