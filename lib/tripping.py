from lib.constants import Calculations
from lib.circulation import circl
from scipy.integrate import solve_ivp
import numpy as np


def trip_operation(F_b, flowrate, length, operation='POOH'):
    clc = Calculations(lengths=length)
    fluid_force = circl(flowrate, length) if flowrate != 0 else np.zeros_like(clc.MD[:-1])
    def dFds(s, F):
        # Interpolate the necessary values for the current s
        DLS_interp = np.interp(s, clc.MD[:-1], clc.DLS)
        n_z_interp = np.interp(s, clc.MD[:-1], clc.n_z)
        b_z_interp = np.interp(s, clc.MD[:-1], clc.b_z)
        t_z_interp = np.interp(s, clc.MD[:-1], clc.t_z)
        bw_pipe_interp = np.interp(s, clc.MD[:-1], clc.bw_pipe)
        fluid_interp = np.interp(s, clc.MD[:-1], fluid_force)

        # Calculate the contact force
        w_c = np.sqrt((F * DLS_interp + bw_pipe_interp * n_z_interp)**2 + (bw_pipe_interp * b_z_interp)**2)

        if operation == 'POOH':
            return -(bw_pipe_interp*t_z_interp + clc.ff*w_c + fluid_interp)
        elif operation == 'RIH':
            return -(bw_pipe_interp*t_z_interp - clc.ff*w_c + fluid_interp)
        elif operation == 'ROB':
            return -(bw_pipe_interp*t_z_interp + fluid_interp)
        else:
            raise ValueError("Unknown operation type specified.")

    s_span = [clc.MD[-1], clc.MD[0]]  # Span from the end to the start
    solution = solve_ivp(dFds,
                         s_span,
                         [F_b],
                         method='RK45',
                         t_eval=np.flip(clc.MD),
                         max_step=5,
                         rtol=0.001,
                         atol=0.001,
                         dense_output=True
    ) 
    
    return solution.y.flatten()