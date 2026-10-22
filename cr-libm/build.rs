//! Compiles the vendored CORE-MATH subset (`third-party/core-math`, upstream commit in
//! `third-party/core-math.ORIGIN`). The XFOIL reference build (`scripts/xfoil-build.sh`) compiles
//! the same files with the same flags, so yFoil and the reference call the same functions.

/// Must equal `CORE_MATH_CFLAGS` in `scripts/xfoil-build.sh`. `-ffp-contract=off` keeps the
/// compiler from fusing multiply-adds the sources do not write as `fma`.
const CORE_MATH_OPT_LEVEL: u32 = 2;
const CORE_MATH_CFLAGS: &[&str] = &["-ffp-contract=off"];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // the `metallic` trial takes the functions from Rust: no C to compile
    if std::env::var_os("CARGO_FEATURE_METALLIC").is_some() {
        return Ok(());
    }
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../third-party/core-math/src/binary64");
    let mut files = Vec::new();
    for entry in std::fs::read_dir(&root)? {
        let dir = entry?.path();
        let name = dir
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or("bad directory name")?
            .to_owned();
        let file = dir.join(format!("{name}.c"));
        if file.is_file() {
            println!("cargo:rerun-if-changed={}", dir.display());
            files.push(file);
        }
    }
    files.sort();
    let mut build = cc::Build::new();
    build.files(files).opt_level(CORE_MATH_OPT_LEVEL).cargo_warnings(false);
    for flag in CORE_MATH_CFLAGS {
        build.flag(flag);
    }
    build.try_compile("core-math")?;
    Ok(())
}
