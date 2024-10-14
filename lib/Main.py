from lib.constants import Calculations as clc
from lib.solver import Main_Func
from lib.visualize import visualize_vel
from lib.init_xls import outputFolderPath
import numpy as np
from scipy.integrate import solve_ivp

def simulation(depth, fric_mod):
    # Normal Force Calculation
    P = {}
    clc(depth)
    
    # Forces = np.zeros(clc.noe)
    # for i in range(clc.noe-1, 1, -1)
    #     F0 = Forces[i+1]
    #     function = lambda F: (-F + F0 + (clc.global_length_array[i]))*(clc.bw_pipe[i]*clc.t_z[i]
    #     + clc.mu_s*np.sqrt((F*clc.DLS[i] + clc.bw_pipe[i]*clc.n_z[i])**2 + (clc.bw_pipe[i]*clc.b_z[i])**2))
    #     x = newton(function, F0)
    #     Forces[i] = x if x > 0 else 0.00001

    P['STATIC_CHECK_PREV'] = np.zeros(clc.noe)
    P['DOWNHOLE_WEIGHT'] = []
    P['DOWNHOLE_TORQUE'] = []
    P['DOC'] = []
    P['SOLUTION_TIME'] = []
    P['THETA_PREV'] = [0.0]
    P['HOLE_DEPTH_PREV'] = [clc.HOLE_DEPTH]
    P['FRICTION_FORCE_STORE'] = []
    P['FRICTION_TORQUE_STORE'] = []
    P['STATIC_CHECK_PREV_STORE'] = []
    P['new_force'] = []

    sol = solve_ivp(
        lambda t, x: Main_Func(t, x, P, fric_mod),
        [0, clc.time_val],
        np.zeros(4 * clc.noe),
        max_step=0.001,
        rtol=1,
        atol=1,
        method="RK45",
        dense_output=False,
    )

    time_arr = sol.t
    sol_arr = sol["y"]
    P['TIME_ARRAY'] = time_arr
    P['SOLUTION_MAIN_ARRAY'] = sol_arr

    visualize_vel(time_arr, sol_arr, clc.noe, P, outputFolderPath, fric_mod)