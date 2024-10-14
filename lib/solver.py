from lib.topdrive import topdrive
from lib.Modules.Coulomb import Coulomb_Friction
from lib.Modules.Stribeck import Stribeck_friction
from lib.Modules.Friction import Friction
from lib.Modules.Bitrock import bit_rock as bit_rock
import numpy as np

def Main_Func(t, x, constants, p, fric_mod):

    """
    This function has the lumped spring mass damper dynamic equations in axial and torsional directions
    The second order ode is broken down to two first order ODE's to solve the problem
    This is the driver function to ODE solver

    Parameters
    ----------
    t : `float`,
        time in seconds

    x : `np.array`,
        global dof vector : [z,v,theta,omega]*noe

    p : `dict`,
        dictionary containing important constants

    Returns
    -------
    dx : `np.array`
        global dof_dot vector : [z_dot, v_dot, theta_dot, omega_dot]*noe

    """
    clc = constants
    MK_top_s = clc.global_mass_inv_ka_matrix
    m_inv_array = 1/clc.global_mass_array
    mca_array = clc.global_ca_array/clc.global_mass_array
    JK_top_s = clc.global_inertia_inv_kt_matrix
    J_inv_array = 1/clc.J
    Jct_array = clc.global_ct_array/clc.J
    dia_pipe_equiv = clc.DIA_EQ
    
    dx = np.zeros(4 * clc.noe)
    Forcing_F = np.zeros(clc.noe)
    Forcing_T = np.zeros(clc.noe)
    _, _, z_top_drive, theta_top_drive = topdrive(t, clc)
    Forcing_F[0] = clc.ka[0] * z_top_drive
    Forcing_T[0] = clc.kt[0] * theta_top_drive

    # Checkpoint for debugging
    checkpoint = float(np.round(t, 2))
    
    mm_omega = np.zeros(clc.noe)
    # if motor_elem != "N":
    #     motor_speed = 2 * np.pi * p[MOTOR_RPG] * p[MOTOR_FLOW_RATE]*(1/60)
    #     mm_omega[motor_elem::] = motor_speed

    z = x[0::4]
    v = x[1::4]
    theta = x[2::4] - mm_omega * t
    omega = x[3::4] - mm_omega
    doc = bit_rock(z[-1], theta[-1], p, clc)

    Forcing_F[-1] = -clc.K_WOB * doc  # - c_bit_axial * v[-1] * abs(np.sign(doc))
    Forcing_T[-1] = -clc.K_TQ * doc

    Friction_force, Friction_torque, p['STATIC_CHECK_PREV'], new_fric_force = Friction(
        z, v, theta, omega, Forcing_F, Forcing_T, p['STATIC_CHECK_PREV'], clc, fric_mod
    )
    
    # if fric_mod.lower() == "stribeck":
    #     Friction_force, Friction_torque, p['STATIC_CHECK_PREV'], new_fric_force = Stribeck_friction(
    #         z, v, theta, omega, Forcing_F, Forcing_T, p['STATIC_CHECK_PREV'], clc
    #     )
    # elif fric_mod.lower() == "coulomb":
    #     Friction_force, Friction_torque, p['STATIC_CHECK_PREV'], new_fric_force = Coulomb_Friction(
    #         z, v, theta, omega, Forcing_F, Forcing_T, p['STATIC_CHECK_PREV'], clc
    #     )

    # store weight , torque , depth of cut and solution time
    # TODO: Dense outputs needn't require storing solution time

    p['new_force'].append(new_fric_force[-1])
    p['DOWNHOLE_WEIGHT'].append(Forcing_F[-1])
    p['DOWNHOLE_TORQUE'].append(Forcing_T[-1])
    p['DOC'].append(doc)
    p['SOLUTION_TIME'].append(t)
    p['FRICTION_FORCE_STORE'].append(Friction_force)
    p['FRICTION_TORQUE_STORE'].append(Friction_torque)
    p['STATIC_CHECK_PREV_STORE'].append(p['STATIC_CHECK_PREV'])

    # Governing Equations
    dx[0::4] = np.array(v)
    dx[1::4] = (
        -MK_top_s.dot(z) + m_inv_array * (Forcing_F - Friction_force) - mca_array * v
    )
    dx[2::4] = np.array(omega + mm_omega)
    dx[3::4] = (
        -JK_top_s.dot(theta)
        + (J_inv_array) * (Forcing_T - Friction_torque)
        - (Jct_array * omega) * ((0.5 * dia_pipe_equiv)**2)
    )
    return dx
