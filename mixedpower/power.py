import numpy as np
from typing import Union
from scipy.stats import t as students_t # type: ignore
from scipy.stats import nct # type: ignore
from scipy.optimize import minimize # type: ignore
import math

from mixedpower.constants import DESIGNS, VARIABLES

def power(
        design:str='CCC', 
        cohens_d:Union[int, np.float32]=0.5,
        resid:Union[int, np.float32]=1.0,
        target_intercept:Union[int, np.float32]=0.05,
        participant_intercept:Union[int, np.float32]=0.05,
        participant_x_target:Union[int, np.float32]=0.05,
        target_slope:Union[int, np.float32]=0.05,
        participant_slope:Union[int, np.float32]=0.05,
        n_participants:Union[int, np.int32]=100,
        n_targets:Union[int, np.int32]=100,
        code:Union[int, np.int32]=1.00,
        alpha:Union[int, np.float32]=0.05
        ):
    """
    Calculate the power of a CCC design.
    
    Parameters
    ----------
    design : str
        The design to calculate power for. Must be one of ['CCC'].
    cohens_d : float
    resid : float
    target_intercept : float
    participant_intercept : float
    participant_x_target : float
    target_slope : float
    participant_slope : float
    n_participants : int
    n_targets : int
    code : float 
    alpha : float
    """
    
    assert design in DESIGNS, f"Design must be one of {DESIGNS}"

    if design in ['CCC', 'CCC'.lower()]:
        return power_ccc(
            cohens_d=cohens_d,
            resid=resid,
            target_intercept=target_intercept,
            participant_intercept=participant_intercept,
            participant_x_target=participant_x_target,
            target_slope=target_slope,
            participant_slope=participant_slope,
            n_participants=n_participants,
            n_targets=n_targets,
            code=code,
            alpha=alpha
        )

def solve(variable:str='n_participants', **kwargs):
    """
    Solve for the number of participants or targets required to achieve a desired power.

    Parameters
    ----------
    variable : str
        The variable to solve for. Must be one of ['n_participants', 'n_targets'].
    """

    assert variable in VARIABLES, f"Variable must be one of {VARIABLES}"
    
    if variable in ['n_participants', 'n_participants'.lower()]:
        return solve_n_participants(**kwargs)
    elif variable in ['n_targets', 'n_targets'.lower()]:
        return solve_n_targets(**kwargs)

def power_ccc(
        cohens_d:Union[int, np.float32]=0.5,
        resid:Union[int, np.float32]=1.0,
        target_intercept:Union[int, np.float32]=0.05,
        participant_intercept:Union[int, np.float32]=0.05,
        participant_x_target:Union[int, np.float32]=0.05,
        target_slope:Union[int, np.float32]=0.05,
        participant_slope:Union[int, np.float32]=0.05,
        n_participants:Union[int, np.int32]=100,
        n_targets:Union[int, np.int32]=100,
        code:Union[int, np.int32]=1.00,
        alpha:Union[int, np.float32]=0.05
        ):
    
    # non-centrality parameter
    ncp_denom = np.sqrt(2.0 * (resid/(n_participants*n_targets) + 
                   2*code**2*target_slope/n_targets + 
                   2*code**2*participant_x_target/n_participants))
    ncp = cohens_d / ncp_denom

    # DoF
    MS_e  = resid
    MS_sc = resid + n_participants * target_slope * code**2
    MS_pc = resid + n_targets * participant_slope * code**2
    dof_numer = (MS_pc + MS_sc - MS_e)**2
    dof_denom = (MS_e**2 / ((n_participants-1)*(n_targets-1))) + (
        MS_sc**2 / (n_targets-1)) + (MS_pc**2 / (n_participants-1))
    df = dof_numer / dof_denom if dof_denom != 0 else 1.0  # avoid zero denom

    alpha_half = alpha / 2.0

    # cir values from the *central* t-distribution
    t_crit_upper = students_t.ppf(1.0 - alpha_half, df) 
    t_crit_lower = students_t.ppf(alpha_half,       df)

    cdf_upper = nct.cdf(t_crit_upper, df, ncp)
    cdf_lower = nct.cdf(t_crit_lower, df, ncp)

    p = (1.0 - cdf_upper) + cdf_lower
    v = np.sum([resid, participant_x_target,  
                target_intercept, participant_intercept, 
                ((code**2)*target_slope), 
                ((code**2)*participant_slope)])
    d_stdz = cohens_d / np.sqrt(v)

    results = {
        'power': p,
        'ncp': ncp,
        'dof': df,
        't_crit_upper': t_crit_upper,
        't_crit_lower': t_crit_lower,
        'cdf_upper': cdf_upper,
        'cdf_lower': cdf_lower,
        'v': v,
        'd_stdz': d_stdz,
    }

    return p, results

def solve_n_participants(
        p:Union[int, np.float32]=0.8,
        cohens_d:Union[int, np.float32]=0.5,
        resid:Union[int, np.float32]=1.0,
        target_intercept:Union[int, np.float32]=0.05,
        participant_intercept:Union[int, np.float32]=0.05,
        participant_x_target:Union[int, np.float32]=0.05,
        target_slope:Union[int, np.float32]=0.05,
        participant_slope:Union[int, np.float32]=0.05,
        n_targets:Union[int, np.int32]=100,
        code:Union[int, np.int32]=1.00,
        alpha:Union[int, np.float32]=0.05
):
    
    def cost(p_candidate: float) -> float:
            current_power, _ = power_ccc(
                cohens_d=cohens_d,
                resid=resid,
                target_intercept=target_intercept,
                participant_intercept=participant_intercept,
                participant_x_target=participant_x_target,
                target_slope=target_slope,
                participant_slope=participant_slope,
                n_participants=p_candidate,
                n_targets=n_targets,
                code=code,
                alpha=alpha
            )
            return (current_power - p)**2
    bounds = [(2.0, None)]
    result = minimize(
        cost,
        x0=[1],
        method='L-BFGS-B',
        bounds=bounds
    )

    if not dict(result)['success']:
        return None, result

    return math.ceil(float(dict(result)['x'][0])), result

def solve_n_targets(
        p:Union[int, np.float32]=0.8,
        cohens_d:Union[int, np.float32]=0.5,
        resid:Union[int, np.float32]=1.0,
        target_intercept:Union[int, np.float32]=0.05,
        participant_intercept:Union[int, np.float32]=0.05,
        participant_x_target:Union[int, np.float32]=0.05,
        target_slope:Union[int, np.float32]=0.05,
        participant_slope:Union[int, np.float32]=0.05,
        n_participants:Union[int, np.int32]=100,
        code:Union[int, np.int32]=1.00,
        alpha:Union[int, np.float32]=0.05
):
    
    def cost(t_candidate: float) -> float:
            current_power, _ = power_ccc(
                cohens_d=cohens_d,
                resid=resid,
                target_intercept=target_intercept,
                participant_intercept=participant_intercept,
                participant_x_target=participant_x_target,
                target_slope=target_slope,
                participant_slope=participant_slope,
                n_participants=n_participants,
                n_targets=t_candidate,
                code=code,
                alpha=alpha
            )
            return (current_power - p)**2
    bounds = [(2.0, None)]
    result = minimize(
        cost,
        x0=[1],
        method='L-BFGS-B',
        bounds=bounds
    )

    if not dict(result)['success']:
        return None, result

    return math.ceil(float(dict(result)['x'][0])), result