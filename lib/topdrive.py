from .init_xls import *

def vel_acc(t, constants):
    ''' 
    This function takes a time array as an input and returns velocity and acceleration at any time.
    
    Input:
    -----------------
        t : np.array,
            time in seconds
    
    Returns:
    -----------------
        v : np.array,
            velocity in ft/s
        a : np.array,
            acceleration in ft/s^2
    '''

    # Time intervals definition
    a1 = constants.a1
    a2 = constants.a2
    a3 = constants.a3
    a4 = constants.a4
    a5 = constants.a5
    a6 = constants.a6
    
    # Velocity definition
    v1 = constants.v1
    v2 = constants.v2

    velocity = np.zeros_like(t)
    acceleration = np.zeros_like(t)

    # Velocity calculation
    velocity[(0 <= t) & (t < a1)] = v1 * (t[(0 <= t) & (t < a1)] / a1)
    velocity[(a1 <= t) & (t < a2)] = v1
    velocity[(a2 <= t) & (t < a3)] = v1 * (a3 - t[(a2 <= t) & (t < a3)]) / (a3 - a2)
    velocity[(a3 <= t) & (t < a4)] = -v2 * (t[(a3 <= t) & (t < a4)] - a3) / (a4 - a3)
    velocity[(a4 <= t) & (t < a5)] = -v2
    velocity[(a5 <= t) & (t <= a6)] = v2 * (t[(a5 <= t) & (t <= a6)] - a6) / (a6 - a5)

    # Acceleration calculation
    acceleration[(0 <= t) & (t < a1)] = v1 / a1
    acceleration[(a1 <= t) & (t < a2)] = 0
    acceleration[(a2 <= t) & (t < a3)] = -v1 / (a3 - a2)
    acceleration[(a3 <= t) & (t < a4)] = -v2 / (a4 - a3)
    acceleration[(a4 <= t) & (t < a5)] = 0
    acceleration[(a5 <= t) & (t <= a6)] = v2 / (a6 - a5)

    return velocity, acceleration