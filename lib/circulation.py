from .constants import *
from .init_xls import GPM2ms

def circl(length, flowrate):
    clc = Calculations(depth=length)
    P_drop_outer, P_drop_inner, v_o, v_i = pres_calc(
        Q=GPM2ms(flowrate),
        rho=clc.mud_density, 
        K=clc.visc_p, 
        tao=clc.tao_y,
        m=clc.m,
        clc=clc
    )
    C_tao_0 = (8 + 4*clc.global_eps + clc.global_eps**2)/8

    # Static_part = (A_i-A_h)*clc.mud_density*clc.g
    Annular_part = (1+C_tao_0)*P_drop_outer*clc.A_h
    Inertial_part = clc.mud_density*(v_o**2*clc.A_h - v_i**2*clc.A_i)*clc.DLS*clc.n_z

    lin_force = Inertial_part - Annular_part
    # fluid_force = np.flip(lin_force)
    # check = P_drop_outer*clc.A_h - P_drop_inner*clc.A_i
    return lin_force