//! The native side of the bit-identity check (`scripts/wasm.sh check`).
//!
//! - `native evaluate <function> <in.bin> <out.bin>`: `yfoil_wasm::evaluate_native` on the
//!   inputs in `in.bin` (little-endian u32 count, then the x, y and z arrays of f64), writing the
//!   results as little-endian f64.
//! - `native bench <geometry.json> <options.json> <repeats>`: times `analyse_json` and prints
//!   the mean milliseconds per call.

use std::time::Instant;

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.iter().map(String::as_str).collect::<Vec<_>>().as_slice() {
        ["evaluate", function, input, output] => {
            let bytes = std::fs::read(input).map_err(|e| format!("{input}: {e}"))?;
            let n = u32::from_le_bytes(bytes[0..4].try_into().map_err(|_| "short input")?) as usize;
            let values: Vec<f64> = bytes[4..]
                .as_chunks::<8>()
                .0
                .iter()
                .map(|c| f64::from_le_bytes(*c))
                .collect();
            if values.len() != 3 * n {
                return Err(format!("{input}: expected {} values, found {}", 3 * n, values.len()));
            }
            let result = yfoil_wasm::evaluate_native(function, &values[..n], &values[n..2 * n], &values[2 * n..])?;
            let out: Vec<u8> = result.iter().flat_map(|v| v.to_le_bytes()).collect();
            std::fs::write(output, out).map_err(|e| format!("{output}: {e}"))
        }
        ["bench", geometry, options, repeats] => {
            let geometry = std::fs::read_to_string(geometry).map_err(|e| e.to_string())?;
            let options = std::fs::read_to_string(options).map_err(|e| e.to_string())?;
            let repeats: usize = repeats.parse().map_err(|e| format!("repeats: {e}"))?;
            yfoil_wasm::analyse_json(&geometry, &options)?;
            let start = Instant::now();
            for _ in 0..repeats {
                yfoil_wasm::analyse_json(&geometry, &options)?;
            }
            println!("{:.3}", start.elapsed().as_secs_f64() * 1e3 / repeats as f64);
            Ok(())
        }
        _ => Err(
            "usage: native evaluate <function> <in.bin> <out.bin> | native bench <geometry> <options> <repeats>".into(),
        ),
    }
}
