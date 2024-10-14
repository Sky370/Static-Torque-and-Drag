# from scipy.integrate import cumtrapz
import numpy as np
from .init_xls import *

def topdrive(t, constants):
    # Time intervals definition
    a1, a2, a3, a4, a5, a6 = (
        constants.a1, constants.a2, constants.a3, 
        constants.a4, constants.a5, constants.a6
    )
    b1, b2, b3, b4, b5, b6 = (
        constants.b1, constants.b2, constants.b3, 
        constants.b4, constants.b5, constants.b6
    )

    # Axial and rotational velocities
    v1 = constants.v1
    v2 = constants.v2
    RPM1 = constants.rpm1 * (2*np.pi/60)
    RPM2 = constants.rpm2 * (2*np.pi/60)

    if 0 <= t < a1:
        ROP_topdrive = ((v1 / a1)) * t
        z_top_drive = ROP_topdrive * 0.5 * t
    elif a1 <= t < a2:
        ROP_topdrive = v1
        z_top_drive = ROP_topdrive * (t - 0.5 * a1)
    elif a2 <= t < a3:
        ROP_topdrive = ((v1) * (a3 - t)) / (a3 - a2)
        z_top_drive = v1 * (a2 - 0.5 * a1) + (0.5 * v1 * (t - a2) * (2 * a3 - a2 - t)) / (a3 - a2)
    elif a3 <= t < a4:
        ROP_topdrive = ((v2) * (a3 - t)) / (a4 - a3)
        z_top_drive = 0.5 * (a3 + a2 - a1) * v1 - (0.5 * (t - a3) ** 2) * (v2 / (a4 - a3))
    elif a4 <= t < a5:
        ROP_topdrive = -(v2)
        z_top_drive = 0.5 * (a3 + a2 - a1) * v1 - 0.5 * (a4 - a3) * v2 - (t - a4) * v2
    elif a5 <= t < a6:
        ROP_topdrive = ((v2) * (t - a6)) / (a6 - a5)
        z_top_drive = 0.5 * (a3 + a2 - a1) * v1 - 0.5 * (a6 - a3 + a5 - a4) * v2
    
    # RPM and theta_top_drive

    if 0 <= t < b1:
        RPM_topdrive = ((RPM1 / b1)) * t
        theta_top_drive = RPM_topdrive * 0.5 * t
    elif b1 <= t < b2:
        RPM_topdrive = RPM1
        theta_top_drive = RPM_topdrive * (t - 0.5 * b1)
    elif b2 <= t < b3:
        RPM_topdrive = ((RPM1) * (b3 - t)) / (b3 - b2)
        theta_top_drive = RPM1 * (b2 - 0.5 * b1) + (0.5 * RPM1 * (t - b2) * (2 * b3 - b2 - t)) / (b3 - b2)
    elif b3 <= t < b4:
        RPM_topdrive = ((RPM2) * (b3 - t)) / (b4 - b3)
        theta_top_drive = 0.5 * (b3 + b2 - b1) * RPM1 - (0.5 * (t - b3) ** 2) * (RPM2 / (b4 - b3))
    elif b4 <= t < b5:
        RPM_topdrive = -(RPM2)
        theta_top_drive = 0.5 * (b3 + b2 - b1) * RPM1 - 0.5 * (b4 - b3) * RPM2 - (t - b4) * RPM2
    elif b5 <= t < b6:
        RPM_topdrive = ((RPM2) * (t - b6)) / (b6 - b5)
        theta_top_drive = 0.5 * (b3 + b2 - b1) * RPM1 - 0.5 * (b6 - b3 + b5 - b4) * RPM2

    z_heave = 0
    ROP_heave = 0
    # if p[HEAVE_STATE] == "Y":
    #     if t <= p[HEAVE_DELAY]:
    #         z_heave = 0
    #         ROP_heave = 0
    #     else:
    #         z_heave = p[HEAVE_AMP] * np.sin(p[HEAVE_OMEGA] * (t - p[HEAVE_DELAY]))
    #         ROP_heave = (
    #             p[HEAVE_AMP]
    #             * p[HEAVE_OMEGA]
    #             * np.cos(p[HEAVE_OMEGA] * (t - p[HEAVE_DELAY]))
    #         )

    out = (ROP_topdrive + ROP_heave, RPM_topdrive, z_top_drive + z_heave, theta_top_drive)

    return out


# def topdrive(t, constants):
#     ''' 
#     This function takes a time array as input and returns velocity, displacement, 
#     and rotational displacement over time.
    
#     Parameters
#     ----------
#     t : np.array,
#         Time in seconds
    
#     Returns
#     -------
#     velocity : np.array,
#         Velocity in ft/s
#     rotational_velocity : np.array,
#         Rotational Velocity in rad/s
#     displacement : np.array,
#         Axial displacement (area under velocity curve)
#     rotational_disp : np.array,
#         Rotational displacement (area under rotational velocity curve)
#     '''

#     # Time intervals definition
#     a1, a2, a3, a4, a5, a6 = (
#         constants.a1, constants.a2, constants.a3, 
#         constants.a4, constants.a5, constants.a6
#     )
#     b1, b2, b3, b4, b5, b6 = (
#         constants.b1, constants.b2, constants.b3, 
#         constants.b4, constants.b5, constants.b6
#     )

#     # Velocity definition
#     v1 = constants.v1   # Axial velocity
#     v2 = constants.v2   # Rotational velocity

#     velocity = np.zeros_like(t)
#     rotational_velocity = np.zeros_like(t)

#     # Axial velocity calculation
#     velocity[(0 <= t) & (t < a1)] = v1 * (t[(0 <= t) & (t < a1)] / a1)
#     velocity[(a1 <= t) & (t < a2)] = v1
#     velocity[(a2 <= t) & (t < a3)] = v1 * (a3 - t[(a2 <= t) & (t < a3)]) / (a3 - a2)
#     velocity[(a3 <= t) & (t < a4)] = -v2 * (t[(a3 <= t) & (t < a4)] - a3) / (a4 - a3)
#     velocity[(a4 <= t) & (t < a5)] = -v2
#     velocity[(a5 <= t) & (t <= a6)] = v2 * (t[(a5 <= t) & (t <= a6)] - a6) / (a6 - a5)

#     # Rotational velocity calculation (similar to axial)
#     rotational_velocity[(0 <= t) & (t < b1)] = v2 * (t[(0 <= t) & (t < b1)] / b1)
#     rotational_velocity[(b1 <= t) & (t < b2)] = v2
#     rotational_velocity[(b2 <= t) & (t < b3)] = v2 * (b3 - t[(b2 <= t) & (t < b3)]) / (b3 - b2)
#     rotational_velocity[(b3 <= t) & (t < b4)] = -v2 * (t[(b3 <= t) & (t < b4)] - b3) / (b4 - b3)
#     rotational_velocity[(b4 <= t) & (t < b5)] = -v2
#     rotational_velocity[(b5 <= t) & (t <= b6)] = v2 * (t[(b5 <= t) & (t <= b6)] - b6) / (b6 - b5)

#     # Calculate displacement and rotational displacement using cumtrapz
#     Cum_disp = cumtrapz(velocity, t, initial=0)  
#     displacements = np.diff(Cum_disp, prepend=0)        # Axial displacement
#     rot_disp = cumtrapz(rotational_velocity, t, initial=0)  
#     rotational_disp = np.diff(rot_disp, prepend=0)      # Rotational displacement

#     return velocity, rotational_velocity, displacements, rotational_disp
