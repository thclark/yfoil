//! Publication figure conventions shared by the validation scripts: the figure is drawn at its final physical size so that it
//! goes into a LaTeX document at natural size (`\includegraphics{...}` with no `width=` and no
//! `\resizebox`), and the text in the SVG is the text that ends up on the page.
//!
//! - Width is `TEXT_WIDTH_PT`, the `\textwidth` of the target template, or less.
//! - Every string is `FONT` at `TICK_PT` (8 pt, the journal minimum) or larger; the SVG carries `font-family` so the SVG
//!   viewer and the PDF converter both use it.
//! - Symbols are typeset the way LaTeX math would set them — italic letters, subscripts,
//!   Greek — by [`Label`], a run of styled segments, so "$C_L$" reads as C with a subscript L,
//!   not "CL". Units are signed ("°"), not spelt out.
//! - [`svg_to_pdf`] converts the SVG with `rsvg-convert`, which embeds the font, so the PDF is
//!   what `\includegraphics` takes without Inkscape in the loop. With the `svg` package use
//!   `\includesvg[inkscapelatex=false]` for the same result.
//!
//! plotters works in CSS px (96 per inch, 1 pt = 4/3 px) and its font "size" is 1.24 × the em
//! it renders (both the SVG `font-size` and its layout use `size / 1.24`), so the two helpers
//! below convert once and nothing else in the crate mentions px.

use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use std::path::Path;

/// `\textwidth` of a single-column AIAA journal page in pt: A4 (210 mm) with 1 in margins is
/// 451.3 pt, letter with the same margins is 468 pt; the smaller fits both and is centred on the
/// larger. Check with `\the\textwidth` in the template and change this one number.
pub const TEXT_WIDTH_PT: f64 = 451.0;
pub const FONT: &str = "Times New Roman";
/// Text sizes, pt: tick values and legends at the journal's minimum for figures (8 pt), axis
/// labels one point larger
pub const TICK_PT: f64 = 8.0;
pub const LEGEND_PT: f64 = 8.0;
pub const AXIS_LABEL_PT: f64 = 9.0;
/// Subscripts and superscripts relative to the base size
const SCRIPT_SCALE: f64 = 0.7;

/// pt → plotters px
pub fn px(pt: f64) -> f64 {
    pt * 4.0 / 3.0
}
pub fn pxi(pt: f64) -> i32 {
    px(pt).round() as i32
}
pub fn pxu(pt: f64) -> u32 {
    px(pt).round() as u32
}
/// pt → the plotters font size that renders an em of that many pt
pub fn font_size(pt: f64) -> f64 {
    px(pt) * 1.24
}
/// Regular text at `pt`
pub fn text(pt: f64) -> TextStyle<'static> {
    (FONT, font_size(pt)).into_font().color(&BLACK)
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Script {
    Normal,
    Sub,
    Sup,
}

/// One run of a label: text in one style
#[derive(Clone)]
pub struct Seg {
    pub text: String,
    pub italic: bool,
    pub script: Script,
}

/// A label made of styled runs, e.g. `C` italic + `L` italic subscript
#[derive(Clone, Default)]
pub struct Label(pub Vec<Seg>);

impl Label {
    pub fn new() -> Self {
        Self(Vec::new())
    }
    /// Upright text
    pub fn t(mut self, s: &str) -> Self {
        self.0.push(Seg {
            text: s.into(),
            italic: false,
            script: Script::Normal,
        });
        self
    }
    /// Italic (a math symbol)
    pub fn i(mut self, s: &str) -> Self {
        self.0.push(Seg {
            text: s.into(),
            italic: true,
            script: Script::Normal,
        });
        self
    }
    /// Italic subscript (a symbol)
    pub fn subi(mut self, s: &str) -> Self {
        self.0.push(Seg {
            text: s.into(),
            italic: true,
            script: Script::Sub,
        });
        self
    }
    /// Upright subscript (a word, not a symbol)
    pub fn subt(mut self, s: &str) -> Self {
        self.0.push(Seg {
            text: s.into(),
            italic: false,
            script: Script::Sub,
        });
        self
    }
    /// Upright superscript
    pub fn sup(mut self, s: &str) -> Self {
        self.0.push(Seg {
            text: s.into(),
            italic: false,
            script: Script::Sup,
        });
        self
    }

    fn style(seg: &Seg, pt: f64) -> FontDesc<'static> {
        let size = match seg.script {
            Script::Normal => pt,
            _ => pt * SCRIPT_SCALE,
        };
        let f = (FONT, font_size(size)).into_font();
        if seg.italic {
            f.style(FontStyle::Italic)
        } else {
            f
        }
    }

    /// A segment's advance in px: its ink width plus a space's worth (0.25 em) for each
    /// leading or trailing space, which SVG would otherwise collapse
    fn advance(seg: &Seg, pt: f64) -> (i32, i32, i32) {
        let font = Self::style(seg, pt);
        let trimmed = seg.text.trim();
        let ink = font.box_size(trimmed).map(|(w, _)| w as i32).unwrap_or(0) + 1;
        let space = (px(pt) * 0.25).round() as i32;
        let lead = (seg.text.len() - seg.text.trim_start().len()) as i32 * space;
        let trail = (seg.text.len() - seg.text.trim_end().len()) as i32 * space;
        (lead, ink, trail)
    }

    /// Width in px at base size `pt`
    pub fn width(&self, pt: f64) -> i32 {
        self.0
            .iter()
            .map(|s| {
                let (a, b, c) = Self::advance(s, pt);
                a + b + c
            })
            .sum()
    }

    /// Height of the base text in px (the ascent above the baseline, near enough one em)
    pub fn height(pt: f64) -> i32 {
        px(pt).round() as i32
    }

    /// Draw with the text's left end (reading direction) at `origin`, baseline-anchored: the
    /// glyphs rise above `origin.1`. `rotated` draws it bottom-to-top for a y axis: the reading
    /// direction is −y, the baseline is at `origin.0` and the glyphs extend to its left.
    pub fn draw<DB: DrawingBackend>(
        &self,
        area: &DrawingArea<DB, plotters::coord::Shift>,
        origin: (i32, i32),
        pt: f64,
        rotated: bool,
    ) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
        let (mut x, mut y) = origin;
        for seg in &self.0 {
            let font = Self::style(seg, pt);
            let (lead, ink, trail) = Self::advance(seg, pt);
            let text = seg.text.trim().to_string();
            if rotated {
                y -= lead;
            } else {
                x += lead;
            }
            let w = ink + trail;
            // baseline shift of scripts, in px, positive downwards
            let shift = match seg.script {
                Script::Normal => 0.0,
                Script::Sub => px(pt) * 0.22,
                Script::Sup => -px(pt) * 0.38,
            }
            .round() as i32;
            let pos = Pos::new(HPos::Left, VPos::Bottom);
            if rotated {
                let style = TextStyle::from(font.transform(FontTransform::Rotate270))
                    .pos(pos)
                    .color(&BLACK);
                area.draw(&Text::new(text, (x + shift, y), style))?;
                y -= w;
            } else {
                let style = TextStyle::from(font).pos(pos).color(&BLACK);
                area.draw(&Text::new(text, (x, y + shift), style))?;
                x += w;
            }
        }
        Ok(())
    }

    /// Draw centred on `centre` (a horizontal label centred at x, or a rotated one at y)
    pub fn draw_centred<DB: DrawingBackend>(
        &self,
        area: &DrawingArea<DB, plotters::coord::Shift>,
        centre: (i32, i32),
        pt: f64,
        rotated: bool,
    ) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
        let w = self.width(pt);
        let origin = if rotated {
            (centre.0, centre.1 + w / 2)
        } else {
            (centre.0 - w / 2, centre.1)
        };
        self.draw(area, origin, pt, rotated)
    }
}

/// `rsvg-convert -f pdf`: the PDF next to the SVG, fonts embedded, page = the SVG's size
/// (1 px = 0.75 pt). Returns the PDF path, or the reason when the converter is missing or failed.
pub fn svg_to_pdf(svg: &Path) -> Result<std::path::PathBuf, String> {
    let pdf = svg.with_extension("pdf");
    let out = std::process::Command::new("rsvg-convert")
        .args(["-f", "pdf", "-o"])
        .arg(&pdf)
        .arg(svg)
        .output()
        .map_err(|e| format!("rsvg-convert not run ({e}); install librsvg for the PDF"))?;
    if out.status.success() {
        Ok(pdf)
    } else {
        Err(String::from_utf8_lossy(&out.stderr).trim().to_string())
    }
}
