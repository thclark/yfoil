//! The minimum set cover of the finite-reachable analysis branches by the candidate cases.
//!
//! Exact minimum set cover is NP-hard in general; here the structure is small (tens of cases,
//! ~1000 branches) and dominated by *essential* cases — a branch taken by exactly one candidate
//! forces that candidate in. After the essentials, the remainder is covered greedily by marginal
//! gain (ties broken towards the case with fewer VISCAL iterations, i.e. the cheaper gate), then
//! pruned: any case whose branches are all covered by the others is dropped, largest first. The
//! result is minimal (no case can be removed); the report says so rather than "minimum".

use crate::Candidate;
use std::collections::{BTreeMap, BTreeSet};

pub fn minimum_cover(candidates: &[&Candidate], target: &BTreeSet<String>) -> Vec<String> {
    // who takes each branch
    let mut takers: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (i, c) in candidates.iter().enumerate() {
        for id in &c.taken {
            if target.contains(id) {
                takers.entry(id.as_str()).or_default().push(i);
            }
        }
    }
    let mut chosen: BTreeSet<usize> = takers.values().filter(|v| v.len() == 1).map(|v| v[0]).collect();
    let covered = |chosen: &BTreeSet<usize>| -> BTreeSet<&str> {
        chosen
            .iter()
            .flat_map(|&i| candidates[i].taken.iter().map(String::as_str))
            .collect()
    };
    // greedy on the remainder
    loop {
        let cov = covered(&chosen);
        let remaining: Vec<&str> = target
            .iter()
            .map(String::as_str)
            .filter(|id| !cov.contains(id))
            .collect();
        if remaining.is_empty() {
            break;
        }
        let best = (0..candidates.len())
            .filter(|i| !chosen.contains(i))
            .map(|i| {
                let gain = remaining.iter().filter(|id| candidates[i].taken.contains(**id)).count();
                (gain, std::cmp::Reverse(cost(candidates[i])), i)
            })
            .max()
            .expect("a candidate covers every target branch by construction");
        assert!(best.0 > 0, "target branch not taken by any candidate");
        chosen.insert(best.2);
    }
    // prune: drop any case whose branches are all covered by the rest, most expensive first
    let mut order: Vec<usize> = chosen.iter().copied().collect();
    order.sort_by_key(|&i| std::cmp::Reverse(cost(candidates[i])));
    for i in order {
        let mut without = chosen.clone();
        without.remove(&i);
        let cov = covered(&without);
        if target.iter().all(|id| cov.contains(id.as_str())) {
            chosen.remove(&i);
        }
    }
    let mut names: Vec<String> = chosen.iter().map(|&i| candidates[i].case.name.clone()).collect();
    names.sort();
    names
}

/// Gate cost: the reference's VISCAL calls × iteration budget (an inviscid case is cheapest).
fn cost(c: &Candidate) -> usize {
    let points = c.case.alphas.len() + c.case.alphas_after_reinit.len() + c.case.cls.len();
    if c.case.inviscid {
        points
    } else {
        points * c.case.max_iterations
    }
}
