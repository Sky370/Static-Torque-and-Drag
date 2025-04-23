from lib.constants import Calculations
from lib.circulation import circl
from scipy.integrate import solve_ivp
import numpy as np

def trip_operation(F_b, T_b, GPM, length, operation='POOH'):
    clc = Calculations(depth=length)
    fluid_force = circl(length, flowrate=GPM)

    def dstate_ds(s, F):
        # Interpolate the necessary values for the current s
        DLS     = np.interp(s, clc.MD[:-1], clc.DLS)
        n_z     = np.interp(s, clc.MD[:-1], clc.n_z)
        b_z     = np.interp(s, clc.MD[:-1], clc.b_z)
        t_z     = np.interp(s, clc.MD[:-1], clc.t_z)
        bw_pipe = np.interp(s, clc.MD[:-1], clc.bw_pipe)
        fluid  = np.interp(s, clc.MD[:-1], fluid_force)
        OD = np.interp(s, clc.MD[:-1], clc.global_od_array)



        # Calculate the contact force
        w_c = np.sqrt((F[0] * DLS + bw_pipe * n_z)**2 + (bw_pipe * b_z)**2)

        if operation == 'POOH':
            dFds = -(bw_pipe*t_z + clc.mu*w_c + fluid)
            dTds = 0.0
        elif operation == 'RIH':
            dFds = -(bw_pipe*t_z - clc.mu*w_c + fluid)
            dTds = 0.0
        elif operation == 'ROB':
            dFds = -(bw_pipe*t_z + fluid)
            dTds = 0.0
        elif operation == 'TQ':
            # dFds = -(bw_pipe*t_z + fluid)
            dFds = 0.0
            dTds = -(clc.mu * w_c * OD / 2)
        else:
            raise ValueError(f"Unknown operation {operation!r}")
        return [dFds, dTds]

    s_span = [clc.MD[-1], clc.MD[0]]             # Integration limits (run from bottom MD to top MD)
    y0 = [F_b, T_b]                              # Initial state: [Weight on bit, torque at bit]

    sol = solve_ivp(
                    dstate_ds,
                    s_span,
                    y0,
                    method='RK45',
                    t_eval=np.flip(clc.MD),
                    max_step=5,
                    rtol=1e-3,
                    atol=1e-3
    )

    # sol.y has shape (2, N): row 0 = F(s), row 1 = T(s)
    F_profile = sol.y[0]
    T_profile = sol.y[1]
    return F_profile, T_profile
