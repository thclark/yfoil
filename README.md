# yFoil - an xFoil-based aerofoil solver

An aerofoil solver based on (and **exactly matching** the output of) XFOIL, with some useful extensions and up-to-date cross-platform builds.

**About this project.** XFOIL is an obvious choice as the most widely validated and universally adopted tool out there. But, it has limitations (see below). This project delivers:
- A core library with solver whose results match XFOIL 6.99 outputs identically (*~1e-09 maximum variation across extensive fuzz testing that covers all logical code branches*).
- Intuitive CLI and plotting tools built around the core.
- Well-defined data schema for i/o, to support automation workflows, FAIR data practices and customised plotting.
- Polar sweep shortcuts in the CLI supporting best convergence.
- An extensive test framework.
- Adoption of continuous integration/deployment tools for consistent production of build artifacts.
- Guaranteed memory-safe operation over millions of runs (built in rust).
- WASM target for the library, making it trivial to call from any other language or in-browser.

**Who are you?** I (the lead maintainer of this repo) am Thomas Clark, part of a small group of researchers led by Sarah Barber at Ostschweizer Fachochschule (OST, the Eastern Switzerland University of Applied Sciences) Institut für Energietechnik (IET). We are experimental and computational fluid dynamicists who work in Wind Energy.

**Who funded this?** We're working with a group at Technische Universität München (TUM) on a joint project "Novel experimental methods for investigating wind turbine rotor blade tip vortex formation and breakdown in atmospheric turbulence​". Our part of the work is funded by the [Swiss National Science Foundation](https://www.snf.ch/en).

**Limitations of XFOIL.** Every aerodynamicist in the world has run XFOIL, or something based on it. However, in 2026, the program as originally built is extremely difficult to use in a professional environment:
 - File formats and automation processes are not compliant with modern FAIR data practices.
 - No version control system exists.
 - No exhaustive consistency testing or validation between builds and platforms, with minimal maintenance of the build process:
   - Ongoing upgrades to libraries like libgfortran require ever-more-difficult rebuilds of the xfoil executables.
   - Different processor architectures like RISC (Apple M-series and Snapdragon) require customisations to build. The community supports such rebuilds but availability is disparate.
 - Memory management. Process escape on non-convergence can be extremely unreliable on some architectures we've used, leading to memory leaks. Run once? No problem. Run a hundred thousand times with solutions not converging? Crashed server.

## Future Improvements

We aim to always maintain consistency with XFOIL whilst introducing:
- A warning system which includes explicit raises of several silently-dropped concerns and edge cases handled in the XFOIL logic.
- Improved robustness of boundary layer calculations, rotational corrections and better estimation of CD developed by [the `RFOIL` program](https://repository.tno.nl/SingleDoc?docId=68350) also based on `XFOIL`.
- Full 360 degree polar correction based on adjusted flat plate technique and convergence metrics ("always get an answer regardless of alpha")
- An unsteady solver based on vortex shedding, giving a much better post-stall behaviour.

## Citing yFoil

If you use yFoil in any research or work leading to a publication or technical report, please cite both yFoil *and* XFOIL.

Years of work went into building and validating XFOIL. Authors, maintainers and users of yFoil owe a greate debt of thanks to Prof. Drela, Harold Youngren and any other contributors to the original XFOIL.

See [citations](docs/citations.md) and the [license](LICENSE) for how to properly credit our months of work (yFoil) and their years (XFOIL).

## License

This repo `yFoil` is licenced under `GPL v3 or later` (SPDX `GPL-3.0-or-later`), compatible with the `GPL v2 or later` license used by XFOIL. This (or another compatible license) is required because:
- translations are considered derivative works under the GPL,
- we modify the original source for the purpose of instrumentation (machine-precision level cross validation fo all subroutines), and
- we include the exact package of source code used for validation within the version-controlled repo.

## Development

### Contributions and the use of AI

All contributions must be made and owned by a human (*automated PRs from bots will be rejected instantly*).

All contributions are subject to "the silicon rule": if you contribute code, you own it and you are responsible for its correctness.

No contributions may contain any 'Co-Authored-By' body which refers to any AI or other tool (see [CLAUDE.md](CLAUDE.md) for the reasons why).

Use of AI is not only acceptable but has been essential to the project. Whether it's appropriate is context-sensitive:
- ✅ **accelerate the normal workflow of a human scientist**, by doing one well-bounded task at a time (*example: instrumenting xfoil and extracting its program flow would have taken weeks for an individual but was straightforward for an LLM*)
- ✅ **share *how* you used AI**, so others can learn and improve
- ✅ **use it in a highly-contrained way**, to improve accuracy and relevance (*example: when adding instrumentation, the workflow orchestrator was instructed to also compile a pristine copy, and cross validate outputs from instrumented xfoil, rejecting any instrumentation that alters outputs*)
- ⚠️ **take great care when using AI to write docs**. Any documentation for normal human consumption, like this README and upper level files in `docs/`, shall be written directly by a human because LLMs are terrible at making material relevant and concise. By contrast, using AIs to maintain records that must be routinely updated, like Architecture Decision Records and the known issues log,  is the difference between useless out-of-date docs and a record that's incredibly valuable for capturing knowledge.
- ❌ **do not vibe-code new features and submit them for review without scientifically rigorous validation** first taking place (example: we built  that's the bit that takes the time.
- ❌ **do not try to hide that you used AI** there's no shame in using it - most of the project has been touched by AI - but let's be transparent.

Claude Fable 5.1 was used heavily for the instrumentation of xfoil, translation of subroutines, and construction of the test harness. But not all at once and therein lies the tale.


### Developer tools

**Rust stack.** You'll need to make sure you have the latest stable version of rust and cargo tools installed.

**Precommit.** Make sure the precommit tools are installed. They're a useful QA gate prior to committing and pushing code. If using an AI agent, ensure it's instructed to run these checks prior to making any commit (this helps reduce token consumption through churn, and helps constrain the agent to improve performance, so double win).

```
pre-commit install && pre-commit install -t commit-msg
pre-commit install-hooks
```
