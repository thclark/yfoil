//! The study definition: one entry per generator family, its default section (the validation
//! set's member, drawn black on every panel), and one panel per parameter with the members
//! that vary it. Every member is a real, addressable section (`yfoil geometry naca <d>` or
//! `yfoil geometry karman-trefftz ...`), so the figure is reproducible from the CLI.

use yfoil::geometry::{Geometry, KarmanTrefftz, Section};

pub const N_NODES: usize = 160;

/// One curve of a panel
pub struct Member {
    /// The value of the varied parameter (orders the colour scale)
    pub value: f64,
    /// Legend entry, e.g. "m = 0.04"
    pub label: String,
    /// The section's designation (or Kármán–Trefftz parameters)
    pub designation: String,
    /// Filesystem-safe name of the geometry file
    pub slug: String,
    pub geometry: Geometry,
}

/// One panel: one parameter varied, the other parameters fixed
pub struct Panel {
    pub slug: &'static str,
    /// Title, e.g. "maximum camber m (p = 0.4, t = 0.12)"
    pub title: String,
    pub members: Vec<Member>,
}

/// One generator family
pub struct SeriesStudy {
    pub slug: &'static str,
    pub title: &'static str,
    /// The series default: the validation set's member, drawn black on top of every panel
    pub default: Member,
    pub panels: Vec<Panel>,
    /// Parameter names of the family, in designation order
    pub parameters: &'static str,
    /// What the family was designed for, one paragraph
    pub description: &'static str,
    /// Why the default is in the validation set, one paragraph
    pub validation_use: &'static str,
    /// `docs/references.bib` keys
    pub references: &'static [&'static str],
}

fn slugify(designation: &str) -> String {
    designation
        .to_ascii_lowercase()
        .replace("naca ", "naca")
        .replace(", a = ", "-a")
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() || c == '-' { c } else { '_' })
        .collect()
}

fn naca(designation: &str, value: f64, label: String) -> Member {
    let s = Section::from_designation(designation).unwrap_or_else(|e| panic!("{designation}: {e}"));
    Member {
        value,
        label,
        slug: slugify(&s.designation),
        designation: s.designation.clone(),
        geometry: s.geometry(N_NODES),
    }
}

fn naca_a(designation: &str, a: f64, value: f64, label: String) -> Member {
    let s = Section::from_designation(designation)
        .unwrap()
        .with_a(a)
        .unwrap_or_else(|e| panic!("{designation} a = {a}: {e}"));
    Member {
        value,
        label,
        slug: slugify(&s.designation),
        designation: s.designation.clone(),
        geometry: s.geometry(N_NODES),
    }
}

fn kt(x_centre: f64, y_centre: f64, te_angle: f64, value: f64, label: String) -> Member {
    let k = KarmanTrefftz::new(x_centre, y_centre, te_angle).unwrap();
    Member {
        value,
        label,
        slug: format!("kt_xc{x_centre}_yc{y_centre}_tau{te_angle}").replace('-', "m"),
        designation: k.designation(),
        geometry: k.geometry(N_NODES),
    }
}

fn panel(slug: &'static str, title: impl Into<String>, mut members: Vec<Member>) -> Panel {
    members.sort_by(|a, b| a.value.partial_cmp(&b.value).unwrap());
    Panel {
        slug,
        title: title.into(),
        members,
    }
}

pub fn series() -> Vec<SeriesStudy> {
    vec![
        SeriesStudy {
            slug: "naca-4-digit",
            title: "NACA 4-digit",
            default: naca("0012", 0.0, "0012".into()),
            parameters: "m (maximum camber), p (its position), t (maximum thickness)",
            description: "The 1933 family of Report 460: a four-term polynomial thickness form on a two-parabola \
                mean line, chosen to fit the best sections of the day (Clark Y, Göttingen 398). It remains the \
                most used section family in general aviation and in codes' self-tests, and XFOIL's `NACA` command \
                generates it.",
            validation_use: "NACA 0012 is the symmetric baseline (Rule 6's XFOIL-independent invariants: \
                C_L = C_M = 0 at α = 0, mirror symmetry at −α) and the CI reference case; NACA 4412 is the \
                cambered 4-digit member and a real wing section (Lednicer lists the Champion 7EC Traveler).",
            references: &["jacobs1933", "abbott1959"],
            panels: vec![
                panel(
                    "camber",
                    "maximum camber m (p = 0.4, t = 0.12)",
                    (0..=9)
                        .map(|m| naca(&format!("{m}412"), m as f64 / 100.0, format!("m = 0.0{m}")))
                        .collect(),
                ),
                panel(
                    "camber-position",
                    "position of maximum camber p (m = 0.04, t = 0.12)",
                    (1..=9)
                        .map(|p| naca(&format!("4{p}12"), p as f64 / 10.0, format!("p = 0.{p}")))
                        .collect(),
                ),
                panel(
                    "thickness",
                    "maximum thickness t (m = 0)",
                    [6, 9, 12, 15, 18, 21, 24, 30]
                        .iter()
                        .map(|t| naca(&format!("00{t:02}"), *t as f64 / 100.0, format!("t = 0.{t:02}")))
                        .collect(),
                ),
            ],
        },
        SeriesStudy {
            slug: "naca-4-digit-modified",
            title: "NACA 4-digit modified",
            default: naca("0012-34", 0.0, "0012-34".into()),
            parameters: "m, p (the 2-digit mean line), t, I (leading-edge radius index), M (position of maximum thickness)",
            description: "Report 492 (1934) varied the leading-edge radius and the position of maximum thickness of \
                the 4-digit form independently, at speeds up to and beyond the compressibility burble: index I = 6 \
                is the 4-digit radius, 0 a sharp nose; M is in tenths of chord. The aft part is a cubic in (1 − x) \
                with a fixed trailing-edge thickness.",
            validation_use: "NACA 0012-34 has a sharper nose (I = 3) and its maximum thickness further aft than \
                the 0012, with the same t: the pair separates the effect of leading-edge radius on transition and \
                stagnation-point behaviour from that of thickness.",
            references: &["stack1934", "abbott1959", "ladson1975"],
            panels: vec![
                panel(
                    "le-radius-index",
                    "leading-edge radius index I (M = 0.4, t = 0.12)",
                    (0..=9)
                        .map(|i| naca(&format!("0012-{i}4"), i as f64, format!("I = {i}")))
                        .collect(),
                ),
                panel(
                    "max-thickness-position",
                    "position of maximum thickness M (I = 3, t = 0.12)",
                    (2..=6)
                        .map(|m| naca(&format!("0012-3{m}"), m as f64 / 10.0, format!("M = 0.{m}")))
                        .collect(),
                ),
                panel(
                    "thickness",
                    "maximum thickness t (I = 3, M = 0.4)",
                    [6, 9, 12, 15, 18, 21]
                        .iter()
                        .map(|t| naca(&format!("00{t:02}-34"), *t as f64 / 100.0, format!("t = 0.{t:02}")))
                        .collect(),
                ),
            ],
        },
        SeriesStudy {
            slug: "naca-5-digit",
            title: "NACA 5-digit",
            default: naca("23018", 0.0, "23018".into()),
            parameters: "cl (design lift coefficient), p (position of maximum camber), reflex, t",
            description: "Report 537 (1935) put the maximum camber unusually far forward — a cubic mean line \
                joined to a straight (or, in the reflexed lines, a second cubic) aft part — to raise the maximum \
                lift with little pitching moment. The 230 line (cl = 0.3, p = 0.15) is by far the most used.",
            validation_use: "NACA 23018 is a real aircraft section (Lednicer lists the Beech 90 King Air root; \
                the 23012 is its tip) and, with the camber concentrated at the nose, loads the leading edge \
                heavily: transition and separation are exercised at the nose rather than mid-chord.",
            references: &["jacobs1935", "abbott1959"],
            panels: vec![
                panel(
                    "design-cl",
                    "design lift coefficient cl (p = 0.15, t = 0.18)",
                    (1..=6)
                        .map(|c| naca(&format!("{c}3018"), c as f64 * 0.15, format!("cl = {:.2}", c as f64 * 0.15)))
                        .collect(),
                ),
                panel(
                    "camber-position",
                    "position of maximum camber p (cl = 0.3, t = 0.18)",
                    (1..=5)
                        .map(|p| naca(&format!("2{p}018"), p as f64 / 20.0, format!("p = {:.2}", p as f64 / 20.0)))
                        .collect(),
                ),
                panel(
                    "reflex",
                    "reflex mean line, position p (cl = 0.3, t = 0.18)",
                    (2..=5)
                        .map(|p| naca(&format!("2{p}118"), p as f64 / 20.0, format!("p = {:.2}", p as f64 / 20.0)))
                        .collect(),
                ),
                panel(
                    "thickness",
                    "maximum thickness t (230 mean line)",
                    [6, 9, 12, 15, 18, 21, 24]
                        .iter()
                        .map(|t| naca(&format!("230{t:02}"), *t as f64 / 100.0, format!("t = 0.{t:02}")))
                        .collect(),
                ),
            ],
        },
        SeriesStudy {
            slug: "naca-16",
            title: "NACA 16-series",
            default: naca("16-212", 0.0, "16-212".into()),
            parameters: "cl (design lift coefficient of the a = 1 mean line), t",
            description: "Stack's 1943 sections (Report 763) were designed to delay the compressibility burble: \
                the 4-digit modified form with I = 4 and maximum thickness at mid-chord (a 16-0tt is a 00tt-45), \
                on the uniform-load a = 1 mean line. They were developed for propellers, and 24 of them were \
                tested at Mach 0.3–0.8 (TN 1546).",
            validation_use: "NACA 16-212 is the compressible-flow member of the set: a section designed and \
                measured at subcritical Mach numbers, where XFOIL's Kármán–Tsien correction applies.",
            references: &["stack1943", "lindsey1948", "ladson1975"],
            panels: vec![
                panel(
                    "design-cl",
                    "design lift coefficient cl (t = 0.12)",
                    (0..=9)
                        .map(|c| naca(&format!("16-{c}12"), c as f64 / 10.0, format!("cl = 0.{c}")))
                        .collect(),
                ),
                panel(
                    "thickness",
                    "maximum thickness t (cl = 0.2)",
                    [6, 9, 12, 15, 18, 21, 24, 30]
                        .iter()
                        .map(|t| naca(&format!("16-2{t:02}"), *t as f64 / 100.0, format!("t = 0.{t:02}")))
                        .collect(),
                ),
            ],
        },
        SeriesStudy {
            slug: "naca-6",
            title: "NACA 6-series",
            default: naca("63-415", 0.0, "63-415".into()),
            parameters: "family (63…67: position of minimum pressure), cl (design lift coefficient), a (extent of uniform loading), t",
            description: "The laminar-flow sections of Report 824 (1945): thickness forms derived by conformal \
                mapping for a prescribed pressure distribution — the family digit is the position of minimum \
                pressure in tenths of chord — on the a mean line, which carries uniform load to x = a and falls \
                linearly to zero at the trailing edge. There is no closed form: the sections are defined by the \
                tabulated ε and ψ mapping functions of TM-4741.",
            validation_use: "NACA 63-415 was the standard wind-turbine section before dedicated families \
                appeared; the Risø catalogue (Bertagnolio et al. 2001) compares measurements with both \
                Navier–Stokes and XFOIL results for it. Its trailing edge closes, so it runs XFOIL's sharp-TE \
                branch, and the 6-series is the family whose generator is *not* a formula.",
            references: &["abbott1945", "abbott1959", "ladson1974", "ladson1996", "carmichael2001", "bertagnolio2001"],
            panels: vec![
                panel(
                    "family",
                    "family: position of minimum pressure (cl = 0.4, t = 0.15)",
                    (3..=7)
                        .map(|f| naca(&format!("6{f}-415"), f as f64 / 10.0, format!("6{f}")))
                        .collect(),
                ),
                panel(
                    "design-cl",
                    "design lift coefficient cl (63, t = 0.15)",
                    [0, 2, 4, 6, 8]
                        .iter()
                        .map(|c| naca(&format!("63-{c}15"), *c as f64 / 10.0, format!("cl = 0.{c}")))
                        .collect(),
                ),
                panel(
                    "loading-extent",
                    "extent of uniform loading a (63-415)",
                    [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
                        .iter()
                        .map(|a| naca_a("63-415", *a, *a, format!("a = {a}")))
                        .collect(),
                ),
                panel(
                    "thickness",
                    "maximum thickness t (63, cl = 0.4)",
                    [6, 9, 12, 15, 18, 21]
                        .iter()
                        .map(|t| naca(&format!("63-4{t:02}"), *t as f64 / 100.0, format!("t = 0.{t:02}")))
                        .collect(),
                ),
            ],
        },
        SeriesStudy {
            slug: "naca-6a",
            title: "NACA 6A-series",
            default: naca("64A010", 0.0, "64A010".into()),
            parameters: "family (63A, 64A, 65A), cl (design lift coefficient of the 6A mean line), t",
            description: "Loftin's 1948 revision (Report 903) of the 6-series removed the cusped trailing edge: \
                the 6A thickness forms are straight-sided from 0.8 chord, and the modified a = 0.8 mean line is \
                straight aft of 0.86 chord, so the cambered sections are too. The trailing-edge angle is finite \
                and the sections are easier to build.",
            validation_use: "NACA 64A010 is the classic transonic and unsteady-aerodynamics test section (the \
                AGARD standard cases): symmetric, sharp trailing edge with a finite angle, straight-sided aft.",
            references: &["loftin1948", "ladson1974", "ladson1996"],
            panels: vec![
                panel(
                    "family",
                    "family (cl = 0, t = 0.10)",
                    (3..=5)
                        .map(|f| naca(&format!("6{f}A010"), f as f64 / 10.0, format!("6{f}A")))
                        .collect(),
                ),
                panel(
                    "design-cl",
                    "design lift coefficient cl (64A, t = 0.10)",
                    [0, 2, 4, 6, 8]
                        .iter()
                        .map(|c| naca(&format!("64A{c}10"), *c as f64 / 10.0, format!("cl = 0.{c}")))
                        .collect(),
                ),
                panel(
                    "thickness",
                    "maximum thickness t (64A, cl = 0)",
                    [6, 8, 10, 12, 15, 18, 21]
                        .iter()
                        .map(|t| naca(&format!("64A0{t:02}"), *t as f64 / 100.0, format!("t = 0.{t:02}")))
                        .collect(),
                ),
            ],
        },
        SeriesStudy {
            slug: "karman-trefftz",
            title: "Kármán–Trefftz",
            default: kt(-0.1, 0.05, 10.0, 0.0, "xc = −0.1, yc = 0.05, τ = 10°".into()),
            parameters: "x_centre (thickness), y_centre (camber), τ (trailing-edge angle; 0 is the Joukowski cusp)",
            description: "The conformal map of a circle (von Kármán and Trefftz 1918), generalising Joukowski's \
                1910 transformation to a finite trailing-edge angle. The potential-flow solution is known in \
                closed form, so the section checks a panel method against an exact answer without a boundary \
                layer in the way. It is not a wing section anyone builds.",
            validation_use: "The analytic member: its exact inviscid pressure distribution is an \
                XFOIL-independent reference for the panel solver, and its sharp trailing edge with a finite \
                angle exercises the SHARP branch on a smooth, non-NACA shape.",
            references: &["karman1918", "joukowski1910"],
            panels: vec![
                panel(
                    "thickness",
                    "circle centre x (yc = 0.05, τ = 10°)",
                    [-0.02, -0.05, -0.1, -0.15, -0.2, -0.25]
                        .iter()
                        .map(|x| kt(*x, 0.05, 10.0, -x, format!("xc = {x}")))
                        .collect(),
                ),
                panel(
                    "camber",
                    "circle centre y (xc = −0.1, τ = 10°)",
                    [0.0, 0.025, 0.05, 0.1, 0.15, 0.2]
                        .iter()
                        .map(|y| kt(-0.1, *y, 10.0, *y, format!("yc = {y}")))
                        .collect(),
                ),
                panel(
                    "te-angle",
                    "trailing-edge angle τ (xc = −0.1, yc = 0.05)",
                    [0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0]
                        .iter()
                        .map(|t| kt(-0.1, 0.05, *t, *t, format!("τ = {t}°")))
                        .collect(),
                ),
            ],
        },
    ]
}

/// The validation set: the sections the XFOIL-equivalence cases are run on, one per family
/// (plus the cambered 4-digit), with why each is there
pub struct ValidationPick {
    pub designation: &'static str,
    pub series: &'static str,
    pub cli: &'static str,
    pub why: &'static str,
    pub references: &'static [&'static str],
}

pub fn validation_set() -> Vec<ValidationPick> {
    vec![
        ValidationPick {
            designation: "NACA 0012",
            series: "4-digit",
            cli: "yfoil geometry naca 0012",
            why: "Symmetric baseline and CI reference case: the XFOIL-independent invariants (C_L = C_M = 0 at α = 0, \
                mirror symmetry at −α) and every solver stage's first fixture.",
            references: &["jacobs1933"],
        },
        ValidationPick {
            designation: "NACA 4412",
            series: "4-digit",
            cli: "yfoil geometry naca 4412",
            why: "The cambered 4-digit member, a real wing section, and the ±15° polar case that runs past the \
                limits of convergence.",
            references: &["jacobs1933", "lednicer"],
        },
        ValidationPick {
            designation: "NACA 0012-34",
            series: "4-digit modified",
            cli: "yfoil geometry naca 0012-34",
            why: "The 0012 thickness with a sharper nose (I = 3) and its maximum thickness aft (M = 0.4): \
                leading-edge radius separated from thickness.",
            references: &["stack1934"],
        },
        ValidationPick {
            designation: "NACA 23018",
            series: "5-digit",
            cli: "yfoil geometry naca 23018",
            why: "A real aircraft root section (Beech 90 King Air) with the camber concentrated at the nose.",
            references: &["jacobs1935", "lednicer"],
        },
        ValidationPick {
            designation: "NACA 16-212",
            series: "16-series",
            cli: "yfoil geometry naca 16-212",
            why: "Designed for high subsonic speed and tested at Mach 0.3–0.8: the compressible (Kármán–Tsien) \
                case.",
            references: &["stack1943", "lindsey1948"],
        },
        ValidationPick {
            designation: "NACA 63-415",
            series: "6-series",
            cli: "yfoil geometry naca 63-415",
            why: "The classic wind-turbine section, with published XFOIL and Navier–Stokes comparisons; a closed \
                trailing edge; the family without a closed-form generator.",
            references: &["abbott1945", "bertagnolio2001"],
        },
        ValidationPick {
            designation: "NACA 64A010",
            series: "6A-series",
            cli: "yfoil geometry naca 64A010",
            why: "The standard transonic test section: symmetric, straight-sided aft, sharp trailing edge with a \
                finite angle.",
            references: &["loftin1948"],
        },
        ValidationPick {
            designation: "Kármán–Trefftz (xc = −0.1, yc = 0.05, τ = 10°)",
            series: "analytic",
            cli: "yfoil geometry karman-trefftz --x-centre -0.1 --y-centre 0.05 --te-angle 10",
            why: "Exact potential-flow solution: the panel method checked against a closed-form answer.",
            references: &["karman1918", "joukowski1910"],
        },
    ]
}
