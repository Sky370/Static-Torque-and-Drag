import numpy as np


def Stribeck_friction(z, v, theta, omega, Forcing_F, Forcing_T, static_check_prev, constants):
    """
    Friction model of the system
    =============================
    This function enforces stribeck friction model on the system with following features :
    a. Friction model is coupled i.e, resultant reaction force from axial and tangential direction is taken into account
    b. No movement untill the resultant reaction force reaches static limit
    c. coloumb friction model is implemented when the drill string is in motion
    d. STATIC_CHECK_PREV is used to track the friction state for each element

    Parameters
    ----------
    z : `float`,
        displacement in meters
    v : `float`,
        velocity in m/sec
    theta : `float`,
        rotation angle in radians
    omega : `float`,
        rotational speed in rads/sec
    Forcing_F : `np.array`,
        Forcing_F vector
    Forcing_T : `np.array`,
        Forcing_T vector
    static_check_prev :`np.array`,
        A vector containing state of friction ( 0 being static friction,
        1 representing dynamic friction ) against each element
    p : `dict`,
        A dictionary containing the constant values

    Returns
    -------
    out : tuple of `np.array`
        Friction_force
        Friction_torque

    """
    clc = constants
    Normal_force = clc.Normal_force
    ca_array = clc.global_ca_array
    KA_top_s = clc.global_ka_matrix
    ct_array = clc.global_ct_array
    KT_top_s = clc.global_kt_matrix
    dia_pipe_equiv = clc.DIA_EQ
    mu_static = clc.mu_s
    mu_dynamic = clc.mu_d
    # motor_elem = p[MOTOR_INDEX]
    
    Friction_limit = mu_static * Normal_force
    v_tan = 0.5 * dia_pipe_equiv * omega
    resultant_vel = np.abs(np.sqrt(v**2 + v_tan**2))

    mu_effective = mu_dynamic + (mu_static - mu_dynamic) * np.exp(-resultant_vel / clc.v_cs)
    coloumb_r = np.where(mu_effective * Normal_force < 1e-10, 0, mu_effective * Normal_force)

    # Calculating forces
    Fd_a = -KA_top_s.dot(z) - ca_array * v + Forcing_F
    Fd_t = (-KT_top_s.dot(theta)/(0.5 * dia_pipe_equiv)        # Converting from Torque to Force_tangential
            - ct_array * omega * (0.5 * dia_pipe_equiv)        # Converting to Force_tangential
            + Forcing_T/(0.5 * dia_pipe_equiv)                 # Converting from Torque to Force_tangential
    )
    Fd_resultant = np.sqrt(Fd_a**2 + Fd_t**2)
    comp_a = np.divide(v, resultant_vel, out=np.zeros_like(v), where=(resultant_vel!=0))
    comp_t = np.divide(v_tan, resultant_vel, out=np.zeros_like(v_tan), where=(resultant_vel!=0))
    

    # Update static checks
    static_check_prev[resultant_vel < 1e-06] = 0
    static_check_prev[Fd_resultant > Friction_limit] = 1    

    # Force and torque calculations
    Friction_force = np.where(static_check_prev == 0, Fd_a, coloumb_r * comp_a)
    Friction_torque = np.where(static_check_prev == 0, 
                               Fd_t * 0.5 * dia_pipe_equiv, 
                               coloumb_r * comp_t * 0.5 * dia_pipe_equiv)

    # # Handle motor element
    # if motor_elem != "N":
    #     Friction_force[motor_elem + 1] = 0
    #     Friction_torque[motor_elem + 1] = 0

    out = (Friction_force, Friction_torque, static_check_prev, Friction_torque/(0.5 * dia_pipe_equiv))

    return out