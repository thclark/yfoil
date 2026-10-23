---
icon: lucide/circle-dot
---

# GitHub issues

How issues on [thclark/yfoil](https://github.com/thclark/yfoil/issues) are labelled, and why. This document
holds the label set and the rules for applying it. Branches and pull requests are covered in
[git-branching.md](git-branching.md) and [git-pull-requests.md](git-pull-requests.md).

## How the label scheme is shaped

Every label is written `class: value`, so an issue's labels read as a short record of what it is, where it stands
and what it touches, and the label list sorts by class.

Each class takes one [Tailwind](https://tailwindcss.com/docs/colors) colour family, and no two classes use
families that sit next to each other in the palette, so a class is recognisable by hue at a glance. Within a
class the shades are ordered to mean something: priority darkens as it rises, flow darkens as an issue moves
towards being ready to work on, and crates darken along the pipeline from the user-facing crate to the outputs.
The tables below list each class from its lightest shade to its darkest.

A new class takes red, amber or rose, the families still free and not adjacent to one in use. Blue is reserved
for a `group:` class, to be added only if a family of related issues ever needs one.

A new crate in the Cargo workspace gets a `crate:` label named after the crate, with any repository prefix
dropped, described "Mostly related to the &lt;crate&gt; &lt;purpose&gt; crate" in at most 100 characters
(GitHub's limit on a label description), and a fuchsia shade placed by where the crate sits in the pipeline.

## Labels

### Kind

What sort of work the issue is. Tailwind emerald.

| Label | Colour | Description |
|---|---|---|
| `kind: devops` | `#d1fae5` | Build, CI, release or repository infrastructure |
| `kind: enhancement` | `#6ee7b7` | An improvement to an existing capability |
| `kind: feature` | `#10b981` | A new capability |
| `kind: epic` | `#047857` | A group of related issues delivered together |
| `kind: bug` | `#064e3b` | Something isn't working |

### Priority

How soon it should be done. Tailwind lime, darkest for the highest priority.

| Label | Colour | Description |
|---|---|---|
| `priority: low` | `#bef264` | Do eventually or wait for user need/urgency |
| `priority: medium` | `#84cc16` | Do in the normal course of work |
| `priority: high` | `#4d7c0f` | Do next |

### Flow

Where the issue stands on its way to being worked on. Tailwind cyan.

| Label | Colour | Description |
|---|---|---|
| `flow: wont-fix` | `#a5f3fc` | Decided not to do this |
| `flow: information-required` | `#67e8f9` | Cannot progress until the reporter or a maintainer supplies more information |
| `flow: decision-required` | `#22d3ee` | Not yet decided whether to do this - awaiting more information or experience |
| `flow: plan-required` | `#0891b2` | A plan is required before this can progress to ai-ready |
| `flow: ai-ready` | `#155e75` | Sufficiently described that an agent can take this and write a PR |

### Area

What the work is for. Tailwind violet.

| Label | Colour | Description |
|---|---|---|
| `area: documentation` | `#ede9fe` | Improvements or additions to documentation |
| `area: validation` | `#ddd6fe` | Engineering assessment of results (series cases, aerofoil selection) |
| `area: testing` | `#c4b5fd` | Tests, fixtures and test infrastructure |
| `area: library` | `#a78bfa` | The analysis library: panel method, boundary layer and solver |
| `area: geometry` | `#8b5cf6` | Aerofoil geometry: series generators, panelling and conversion |
| `area: cli` | `#7c3aed` | The yfoil command-line interface and its input/output |
| `area: plot` | `#6d28d9` | Plotting and figures |
| `area: automation` | `#5b21b6` | Scripts, workflows and agent automation |
| `area: marketing` | `#4c1d95` | Website, examples, showcases and outreach |

`area: validation` and `area: testing` follow the repository's split between the engineering assessment and the
consistency tracking ([testing.md](testing.md), "Validation and consistency"): `validation` is the assessment of
yFoil's results for aerodynamicists (`docs/validation/`), `testing` is the tests, the fixtures and consistency with
XFOIL (`docs/development/testing/`).

### Crate

Where the work lands: one label per Cargo workspace member, in pipeline order from the user-facing crate to the
outputs. Tailwind fuchsia.

| Label | Colour | Description |
|---|---|---|
| `crate: yfoil` | `#fae8ff` | Mostly related to the yfoil aerofoil analysis library and CLI crate |
| `crate: cr-libm` | `#f5d0fe` | Mostly related to the cr-libm correctly rounded transcendentals (CORE-MATH) crate |
| `crate: xtask` | `#f0abfc` | Mostly related to the xtask reference build, fixture and report automation crate |
| `crate: xfoil-instrumentation-check` | `#e879f9` | Mostly related to the xfoil-instrumentation-check instrumented-vs-pristine reference crate |
| `crate: input-sensitivity` | `#d946ef` | Mostly related to the input-sensitivity XFOIL input-perturbation study crate |
| `crate: branch-coverage` | `#c026d3` | Mostly related to the branch-coverage minimal branch-covering case study crate |
| `crate: series-cases` | `#a21caf` | Mostly related to the series-cases ±30° polar study crate |
| `crate: known-issues` | `#86198f` | Mostly related to the known-issues XFOIL known-issues study crate |
| `crate: aerofoil-series` | `#701a75` | Mostly related to the aerofoil-series generator-family validation study crate |

## Applying labels

### Open issues

An open issue carries one `kind`, at least one `area`, one `priority` and one `flow` label, plus a `crate` label for
each crate it is expected to change.

`area` says what the work is for and `crate` says where it lands, so the two combine freely: work on a study crate
done for the showcase is `crate: series-cases` and `area: marketing` as well as whatever else it is. An issue that
changes no crate, such as CI work or an ideas inbox, carries no `crate` label.

`kind: epic` marks an issue that groups related issues.

### Flow

An issue enters the flow at one of two stages, then advances:

1. `flow: decision-required` when it is not yet decided whether to do it at all, or `flow: information-required` when
   it cannot progress without more detail from the reporter or a maintainer;
2. `flow: plan-required` once it is to be done and needs a plan;
3. `flow: ai-ready` once it is described well enough that an agent can take it and write a pull request.

`flow: wont-fix` ends an issue. Anything that adds a new data structure or public entity starts at
`flow: plan-required`, however clear it looks.

A feature that would deliberately diverge from XFOIL starts at `flow: decision-required`: yFoil translates XFOIL and
never modifies its algorithm, and a deliberate divergence is allowed only with a row in the divergence table of
CLAUDE.md Rule 2, which needs the same evidence as a bug report against XFOIL.

### Closed issues

A closed issue keeps its `kind`, `area` and `crate` labels and loses its `priority` and `flow` labels: priority and
flow stage mean nothing once an issue is closed. The labels describe what the issue is about, not what was done
about it, so an issue about the CLI that is decided against still carries `area: cli` and `crate: yfoil`.

An issue closed as out of scope also carries `flow: wont-fix`, the one `flow` label that outlives the issue, and is
closed with GitHub's "not planned" reason.
