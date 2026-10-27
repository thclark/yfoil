//! yFoil's own twins: the conditioning check every analysed point carries unless `--no-twins`.
//!
//! A **twin** of a point is the same solve — the same OPER command from the same initial state —
//! on panels jittered by −1, 0 or +1 ULP per coordinate (`crate::geometry::jitter`). Five twins
//! per point, seeds 1–5, the instrument the validation studies use on the reference
//! (`docs/development/conventions/terminology.md`). Where a twin lands far from the point, the
//! point is balanced on a knife edge: it is *ill-conditioned*, and the output layer withholds it
//! (`crate::output::conditioning`).
//!
//! **Per point, given its initial state.** A twin starts from a copy of the session as the point
//! finds it — the predecessor's converged boundary layer, pointers and flags — with only the
//! airfoil panels swapped for the jittered ones and the quantities built from them on entry
//! (GGCALC's inviscid system, the DIJ) rebuilt. So a point's verdict depends on exactly what its
//! content id hashes (inputs and predecessor), never on twin histories that parted at an earlier
//! point of a sweep: it is the same whether the point is reached in a polar or restarted from
//! its predecessor. The base solve is never steered by its twins.

use std::sync::Arc;

use crate::geometry::jitter::{jitter_geometry, TWIN_SEEDS};
use crate::geometry::{panel_foil, Geometry, PanelledFoil};

/// The jittered panels of a geometry's five twins, built once per analysis.
#[derive(Debug, Clone)]
pub struct Twins {
    /// One (seed, panels) pair per twin
    pub(crate) panels: Arc<Vec<(u64, PanelledFoil)>>,
}

impl Twins {
    /// The five twins of `geometry` (the coordinates the analysis panels with `panel_foil`).
    pub fn new(geometry: &Geometry) -> Self {
        Self::with_seeds(geometry, &TWIN_SEEDS)
    }

    /// Twins under the given seeds. Seed 0 is not special: an *identity* twin, one whose panels
    /// are those of the base, is made with [`Twins::identity`].
    pub fn with_seeds(geometry: &Geometry, seeds: &[u64]) -> Self {
        let panels = seeds
            .iter()
            .map(|&seed| (seed, panel_foil(&jitter_geometry(geometry, seed))))
            .collect();
        Self {
            panels: Arc::new(panels),
        }
    }

    /// One twin with the base's own panels: it must reproduce the base bit for bit, which is the
    /// proof that a twin session carries the whole of the state it is copied from.
    pub fn identity(geometry: &Geometry) -> Self {
        Self {
            panels: Arc::new(vec![(0, panel_foil(geometry))]),
        }
    }

    /// The seeds, in order
    pub fn seeds(&self) -> Vec<u64> {
        self.panels.iter().map(|(s, _)| *s).collect()
    }
}
