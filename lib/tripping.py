from lib.constants import Calculations
from lib.init_xls import get_area
from lib.topdrive import vel_acc
from lib.circulation import circl
from scipy.integrate import solve_ivp
import numpy as np


def trip_operation(F_b, flowrate, length, operation='POOH'):
    clc = Calculations(depth=length)
    fluid_force = circl(flowrate, length) if flowrate != 0 else np.zeros_like(clc.MD[:-1])

    def dFds(s, F, velocity, acceleration):
        cross_area = np.pi*(clc.global_od_array**2 - clc.global_id_array**2)/4
        
        # Interpolate the necessary values for the current s
        DLS_interp = np.interp(s, clc.MD[:-1], clc.DLS)
        n_z_interp = np.interp(s, clc.MD[:-1], clc.n_z)
        b_z_interp = np.interp(s, clc.MD[:-1], clc.b_z)
        t_z_interp = np.interp(s, clc.MD[:-1], clc.t_z)
        bw_pipe_interp = np.interp(s, clc.MD[:-1], clc.bw_pipe)
        fluid_interp = np.interp(s, clc.MD[:-1], fluid_force)

        # Calculate the contact force with new dynamic terms
        vel_term = 12/231*clc.pipe_density * get_area(s, clc.MD, cross_area) * velocity**2                       # Velocity term
        acc_term = 12/231*clc.pipe_density * get_area(s, clc.MD, cross_area) * acceleration                      # Acceleration term
        w_c = np.sqrt((F * DLS_interp + bw_pipe_interp * n_z_interp - vel_term)**2 + (bw_pipe_interp * b_z_interp)**2)

        if operation == 'POOH':
            return -(bw_pipe_interp*t_z_interp + clc.ff*w_c + fluid_interp + acc_term)
        elif operation == 'RIH':
            return -(bw_pipe_interp*t_z_interp - clc.ff*w_c + fluid_interp + acc_term)
        elif operation == 'ROB':
            return -(bw_pipe_interp*t_z_interp + fluid_interp)
        else:
            raise ValueError("Unknown operation type specified.")

    s_span = [clc.MD[-1], clc.MD[0]]  # Span from the end to the start
        
    if operation == 'POOH':
        times = np.arange(0, clc.a3, 0.5)           # Time domain
        vel, acc = vel_acc(times, clc)                   # Computing velocity and acceleration
        
    elif operation == 'RIH':
        times = np.arange(clc.a3, clc.a6, 0.5)      # Time domain
        vel, acc = vel_acc(times, clc)                   # Computing velocity and acceleration

    elif operation == 'ROB':
        times = np.arange(0, 30, 0.5)               # Time domain
        vel, acc = (
            np.zeros(len(times)), 
            np.zeros(len(times))
            )                                       # Computing velocity and acceleration
    else:
        raise ValueError("Unknown operation type specified.")
    
    sol = [solve_ivp(dFds, 
                    s_span, 
                    [F_b], 
                    args=(vel[i], acc[i]), 
                    method='RK45', 
                    t_eval=np.flip(clc.MD),
                    max_step=10,
                    rtol=0.001,
                    atol=0.001,
                    dense_output=True
                    ).y.flatten()[-1] for i,_ in enumerate(times)]
    
    return sol, times