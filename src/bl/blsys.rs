//! BLSYS and TESYS (xblsys.f): assemble the local 4×5 Newton system for the current interval
//! from the "1" and "2" station states, exactly as XFOIL sequences BLVAR/BLMID/TRDIF/BLDIF,
//! the similarity-station folding, and the conversion of the Ue columns to incompressible Uei.

use crate::bl::system::{BLLocalSystem, BLStationState, FlowParameters, FlowRegime, MidpointCf, TransitionLocation};

/// XFOIL's interval flags (XBL.INC): SIMI, TRAN, TURB, WAKE.
#[derive(Debug, Clone, Copy, Default)]
pub struct IntervalFlags {
    pub simi: bool,
    pub tran: bool,
    pub turb: bool,
    pub wake: bool,
}

/// BLSYS. `s1` is mutable because at the similarity station XFOIL copies COM1 = COM2.
/// `trans` must be `Some` when `flags.tran`.
pub fn blsys(
    sys: &mut BLLocalSystem,
    s1: &mut BLStationState,
    s2: &mut BLStationState,
    flags: IntervalFlags,
    trans: Option<&TransitionLocation>,
    acrit: f64,
    params: &FlowParameters,
) {
    // calculate secondary BL variables and their sensitivities
    let ityp = if flags.wake {
        FlowRegime::Wake
    } else if flags.turb || flags.tran {
        FlowRegime::Turbulent
    } else {
        FlowRegime::Laminar
    };
    s2.blvar(ityp, params);

    // for the similarity station, "1" and "2" variables are the same
    if flags.simi {
        *s1 = s2.clone();
    }
    // BLMID (midpoint Cf) — reads the "1" state, which at SIMI is now the "2" state
    let cfm = MidpointCf::compute(s1, s2, ityp, flags.simi);

    // set up appropriate finite difference system for current interval
    if flags.tran {
        sys.trdif(
            s1,
            s2,
            trans.expect("TRAN requires the transition location"),
            acrit,
            params,
        );
    } else if flags.simi {
        sys.bldif(s1, s2, &cfm, FlowRegime::Laminar, true, acrit, params.idampv);
    // BLDIF(0)
    } else if !flags.turb {
        sys.bldif(s1, s2, &cfm, FlowRegime::Laminar, false, acrit, params.idampv);
    // BLDIF(1)
    } else if flags.wake {
        sys.bldif(s1, s2, &cfm, FlowRegime::Wake, false, acrit, params.idampv); // BLDIF(3)
    } else {
        sys.bldif(s1, s2, &cfm, FlowRegime::Turbulent, false, acrit, params.idampv);
        // BLDIF(2)
    }

    if flags.simi {
        // at similarity station, "1" variables are really "2" variables
        for k in 0..4 {
            for l in 0..5 {
                sys.vs2[k][l] += sys.vs1[k][l];
                sys.vs1[k][l] = 0.0;
            }
        }
    }

    // change system over into incompressible Uei and Mach
    for k in 0..4 {
        let res_u1 = sys.vs1[k][3];
        let res_u2 = sys.vs2[k][3];
        let res_ms = sys.vsm[k];
        sys.vs1[k][3] = res_u1 * s1.u_uei;
        sys.vs2[k][3] = res_u2 * s2.u_uei;
        sys.vsm[k] = res_u1 * s1.u_ms + res_u2 * s2.u_ms + res_ms;
    }
}

/// TESYS(CTE, TTE, DTE): the "dummy" system between the airfoil TE point and the first wake
/// point. Calls BLVAR(3) first, as XFOIL does; no Uei conversion is applied.
pub fn tesys(sys: &mut BLLocalSystem, s2: &mut BLStationState, cte: f64, tte: f64, dte: f64, params: &FlowParameters) {
    for k in 0..4 {
        sys.vsrez[k] = 0.0;
        sys.vsm[k] = 0.0;
        sys.vsr[k] = 0.0;
        sys.vsx[k] = 0.0;
        for l in 0..5 {
            sys.vs1[k][l] = 0.0;
            sys.vs2[k][l] = 0.0;
        }
    }
    s2.blvar(FlowRegime::Wake, params);

    sys.vs1[0][0] = -1.0;
    sys.vs2[0][0] = 1.0;
    sys.vsrez[0] = cte - s2.ctau;

    sys.vs1[1][1] = -1.0;
    sys.vs2[1][1] = 1.0;
    sys.vsrez[1] = tte - s2.theta;

    sys.vs1[2][2] = -1.0;
    sys.vs2[2][2] = 1.0;
    sys.vsrez[2] = dte - s2.dstar - s2.dw;
}
