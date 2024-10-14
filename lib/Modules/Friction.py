import numpy as np


def Friction(z, v, theta, omega, Forcing_F, Forcing_T, static_check_prev, constants, fric_mod):
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
    
    radius = 0.5 * dia_pipe_equiv      # Radius in inches
    epsilon = 1e-6                     # Small threshold for zero velocity detection
    Friction_limit = mu_static * Normal_force
    Force_dynamic = mu_dynamic * Normal_force
    v_tan = radius * omega
    resultant_vel = np.abs(np.sqrt(v**2 + v_tan**2))

    mu_effective = mu_dynamic + (mu_static - mu_dynamic) * np.exp(-resultant_vel / clc.v_cs)
    coloumb_r = np.where(mu_effective * Normal_force < 1e-10, 0, mu_effective * Normal_force)

    # Calculating forces
    # Fd_a = -KA_top_s.dot(z) - ca_array * v + Forcing_F
    # Fd_t = (-KT_top_s.dot(theta)/(0.0254 * 0.5 * dia_pipe_equiv)        # Converting from Torque to Force_tangential
    #         - ct_array * omega * (0.0254 * 0.5 * dia_pipe_equiv)        # Converting to Force_tangential
    #         + Forcing_T/(0.0254 * 0.5 * dia_pipe_equiv)                 # Converting from Torque to Force_tangential
    # )

    Fd_a = -KA_top_s.dot(z) - ca_array * v + Forcing_F
    Fd_t = -KT_top_s.dot(theta)/radius - ct_array*omega*radius + Forcing_T/radius
    Fd_resultant = np.sqrt(Fd_a**2 + Fd_t**2)

    comp_a = np.divide(v, resultant_vel, out=np.zeros_like(v), where=(resultant_vel!=0))
    comp_t = np.divide(v_tan, resultant_vel, out=np.zeros_like(v_tan), where=(resultant_vel!=0))
    
    if fric_mod.lower() == 'stribeck':
        # Update static checks
        static_check_prev[resultant_vel < 1e-06] = 0
        static_check_prev[Fd_resultant > Friction_limit] = 1    

        # Force and torque calculations
        Friction_force = np.where(static_check_prev == 0, Fd_a, coloumb_r * comp_a)
        Friction_torque = np.where(static_check_prev == 0, Fd_t * radius, coloumb_r * comp_t * radius)

        # # Handle motor element
        # if motor_elem != "N":
        #     Friction_force[motor_elem + 1] = 0
        #     Friction_torque[motor_elem + 1] = 0

        out = (Friction_force, Friction_torque, static_check_prev, Friction_torque/radius)
    
    elif fric_mod.lower() == 'coulomb':
        # Update static check based on velocities
        static_check_prev[resultant_vel < epsilon] = 0
        static_check_prev[(resultant_vel >= epsilon) & (resultant_vel < clc.v_cs)] = 1
        static_check_prev[resultant_vel >= clc.v_cs] = 1

        # Combined calculation of friction force
        Friction_force = np.where(
            static_check_prev == 0,
            np.clip(Fd_a, -Friction_limit, Friction_limit),                         # Full static 
            np.where(
                static_check_prev == 1,
                np.minimum(np.abs(Fd_a), Friction_limit) * np.sign(v),              # Slow move in static
                Force_dynamic * comp_a                                              # Full dynamic
            )
        )
        Friction_torque = np.where(
            static_check_prev == 0,
            np.clip(Fd_t, -Friction_limit, Friction_limit) * radius,                # Full static 
            np.where(
                static_check_prev == 1,
                np.minimum(np.abs(Fd_t), Friction_limit)*np.sign(v_tan)*radius,     # Slow move in static
                Force_dynamic * comp_t * radius                                     # Full dynamic
            )
        )

        # # Handle motor element
        # if motor_elem != "N":
        #     Friction_force[motor_elem + 1] = 0
        #     Friction_torque[motor_elem + 1] = 0

        out = Friction_force, Friction_torque, static_check_prev, Friction_torque/radius
    else:
        raise ValueError("Unknown friction type specified.")
    return out