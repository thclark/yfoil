//! LUDCMP / BAKSUB (xsolve.f): dense LU factorisation with implicit scaling and partial
//! pivoting, ported verbatim so that AIJ's factorisation and every back-substitution through
//! it (GGCALC, QDCALC) follow XFOIL's operation order exactly. 1-based storage.

/// LU factors of an n×n matrix stored 1-based as `a[i][j]`, with the pivot index `indx`.
#[derive(Debug, Clone)]
pub struct LuFactors {
    pub n: usize,
    pub lu: Vec<Vec<f64>>,
    pub pivots: Vec<usize>,
}

/// LUDCMP: factor `a` (1-based, n×n, modified in place) into LU form.
#[doc(alias = "LUDCMP")]
pub fn lu_decompose(n: usize, mut a: Vec<Vec<f64>>) -> LuFactors {
    let mut vv = vec![0.0; n + 1];
    let mut indx = vec![0usize; n + 1];

    for i in 1..=n {
        let mut aamax = 0.0_f64;
        for j in 1..=n {
            aamax = a[i][j].abs().max(aamax);
        }
        vv[i] = 1.0 / aamax;
    }

    let mut imax = 1;
    for j in 1..=n {
        for i in 1..j {
            let mut sum = a[i][j];
            for k in 1..i {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
        }

        let mut aamax = 0.0_f64;
        for i in j..=n {
            let mut sum = a[i][j];
            for k in 1..j {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
            let dum = vv[i] * sum.abs();
            if dum >= aamax {
                imax = i;
                aamax = dum;
            }
        }

        if j != imax {
            for k in 1..=n {
                let dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if j != n {
            let dum = 1.0 / a[j][j];
            for i in (j + 1)..=n {
                a[i][j] *= dum;
            }
        }
    }
    LuFactors { n, lu: a, pivots: indx }
}

/// BAKSUB: solve L·U·x = b in place (`b` is 1-based, length n+1).
#[doc(alias = "BAKSUB")]
pub fn lu_back_substitute(lu: &LuFactors, b: &mut [f64]) {
    let n = lu.n;
    let a = &lu.lu;
    let mut ii = 0usize;
    for i in 1..=n {
        let ll = lu.pivots[i];
        let mut sum = b[ll];
        b[ll] = b[i];
        if ii != 0 {
            for j in ii..i {
                sum -= a[i][j] * b[j];
            }
        } else if sum != 0.0 {
            ii = i;
        }
        b[i] = sum;
    }
    for i in (1..=n).rev() {
        let mut sum = b[i];
        if i < n {
            for j in (i + 1)..=n {
                sum -= a[i][j] * b[j];
            }
        }
        b[i] = sum / a[i][i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_a_small_system() {
        // 3x3 with pivoting needed
        let a = vec![
            vec![0.0; 4],
            vec![0.0, 0.0, 2.0, 1.0],
            vec![0.0, 1.0, 1.0, 1.0],
            vec![0.0, 4.0, 3.0, 2.0],
        ];
        let lu = lu_decompose(3, a);
        let mut b = vec![0.0, 5.0, 6.0, 20.0]; // x = (3, 2, 1)
        lu_back_substitute(&lu, &mut b);
        assert!(
            (b[1] - 3.0).abs() < 1e-12 && (b[2] - 2.0).abs() < 1e-12 && (b[3] - 1.0).abs() < 1e-12,
            "{b:?}"
        );
    }
}
