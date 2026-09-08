//! BLSYS and TESYS (xblsys.f): assemble the local 4×5 Newton system for the current interval
//! from the "1" and "2" station states, exactly as XFOIL sequences BLVAR/BLMID/TRDIF/BLDIF,
//! the similarity-station folding, and the conversion of the Ue columns to incompressible Uei.

use crate::bl::system::{FlowParameters, FlowRegime, IntervalSystem, MidpointCf, StationState, Transition};

/// Mirrors XFOIL's `V_INT`.
///
/// XFOIL's interval flags (XBL.INC): SIMI, TRAN, TURB, WAKE.
#[derive(Debug, Clone, Copy, Default)]
#[doc(alias = "V_INT")]
pub struct IntervalFlags {
    pub similarity: bool,
    pub transition: bool,
    pub turbulent: bool,
    pub wake: bool,
}

/// BLSYS. `s1` is mutable because at the similarity station XFOIL copies COM1 = COM2.
/// `trans` must be `Some` when `flags.tran`.
#[doc(alias = "BLSYS")]
pub fn assemble_interval_system(
    sys: &mut IntervalSystem,
    s1: &mut StationState,
    s2: &mut StationState,
    flags: IntervalFlags,
    trans: Option<&Transition>,
    acrit: f64,
    params: &FlowParameters,
) {
    // calculate secondary BL variables and their sensitivities
    let ityp = if flags.wake {
        FlowRegime::Wake
    } else if flags.turbulent || flags.transition {
        FlowRegime::Turbulent
    } else {
        FlowRegime::Laminar
    };
    s2.set_closure_variables(ityp, params);

    // for the similarity station, "1" and "2" variables are the same
    if flags.similarity {
        *s1 = s2.clone();
    }
    // BLMID (midpoint Cf) — reads the "1" state, which at SIMI is now the "2" state
    let cfm = MidpointCf::compute(s1, s2, ityp, flags.similarity);

    // set up appropriate finite difference system for current interval
    if flags.transition {
        sys.assemble_transition_equations(
            s1,
            s2,
            trans.expect("TRAN requires the transition location"),
            acrit,
            params,
        );
    } else if flags.similarity {
        sys.assemble_interval_equations(
            s1,
            s2,
            &cfm,
            FlowRegime::Laminar,
            true,
            acrit,
            params.amplification_model,
        );
    // BLDIF(0)
    } else if !flags.turbulent {
        sys.assemble_interval_equations(
            s1,
            s2,
            &cfm,
            FlowRegime::Laminar,
            false,
            acrit,
            params.amplification_model,
        );
    // BLDIF(1)
    } else if flags.wake {
        sys.assemble_interval_equations(s1, s2, &cfm, FlowRegime::Wake, false, acrit, params.amplification_model);
    // BLDIF(3)
    } else {
        sys.assemble_interval_equations(
            s1,
            s2,
            &cfm,
            FlowRegime::Turbulent,
            false,
            acrit,
            params.amplification_model,
        );
        // BLDIF(2)
    }

    if flags.similarity {
        // at similarity station, "1" variables are really "2" variables
        for k in 0..4 {
            for l in 0..5 {
                sys.jacobian_station2[k][l] += sys.jacobian_station1[k][l];
                sys.jacobian_station1[k][l] = 0.0;
            }
        }
    }

    // change system over into incompressible Uei and Mach
    for k in 0..4 {
        let res_u1 = sys.jacobian_station1[k][3];
        let res_u2 = sys.jacobian_station2[k][3];
        let res_ms = sys.residual_d_machsqd[k];
        sys.jacobian_station1[k][3] = res_u1 * s1.ue_d_uei;
        sys.jacobian_station2[k][3] = res_u2 * s2.ue_d_uei;
        sys.residual_d_machsqd[k] = res_u1 * s1.ue_d_machsqd + res_u2 * s2.ue_d_machsqd + res_ms;
    }
}

/// TESYS(CTE, TTE, DTE): the "dummy" system between the airfoil TE point and the first wake
/// point. Calls BLVAR(3) first, as XFOIL does; no Uei conversion is applied.
#[doc(alias = "TESYS")]
pub fn assemble_te_system(
    sys: &mut IntervalSystem,
    s2: &mut StationState,
    cte: f64,
    tte: f64,
    dte: f64,
    params: &FlowParameters,
) {
    for k in 0..4 {
        sys.residual[k] = 0.0;
        sys.residual_d_machsqd[k] = 0.0;
        sys.residual_d_re[k] = 0.0;
        sys.residual_d_xi[k] = 0.0;
        for l in 0..5 {
            sys.jacobian_station1[k][l] = 0.0;
            sys.jacobian_station2[k][l] = 0.0;
        }
    }
    s2.set_closure_variables(FlowRegime::Wake, params);

    sys.jacobian_station1[0][0] = -1.0;
    sys.jacobian_station2[0][0] = 1.0;
    sys.residual[0] = cte - s2.sqrtctau;

    sys.jacobian_station1[1][1] = -1.0;
    sys.jacobian_station2[1][1] = 1.0;
    sys.residual[1] = tte - s2.theta;

    sys.jacobian_station1[2][2] = -1.0;
    sys.jacobian_station2[2][2] = 1.0;
    sys.residual[2] = dte - s2.dstar - s2.wake_gap;
}
