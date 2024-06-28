from lib1.constants_opt import *

def circl(flowrate, length):
    clc = Calculations(lengths=length)
    P_drop_inner , P_drop_outer, v_i, v_o = p_drop(
        Q=flowrate,
        rho=clc.mud_density_ppg, 
        mu_p=clc.visc_p, 
        tao=clc.tao_y,
        D_o=clc.global_od_array, 
        D_i=clc.global_id_array,
        D_w=clc.HOLE_ARRAY,
        del_L=clc.global_length_array
    )
    Area_o = np.pi*(clc.global_od_array**2)/4
    Area_i = np.pi*(clc.global_id_array**2)/4
    C_tao_0 = (8 + 4*clc.global_eps + clc.global_eps**2)/8

    Static_part = (12/231)*(Area_i-Area_o)*clc.mud_density_ppg*clc.gravity
    Annular_part = (1+C_tao_0)*P_drop_outer*Area_o
    Inertial_part = (12/231)*clc.mud_density_ppg*(v_o**2*Area_o - v_i**2*Area_i)*clc.DLS*clc.n_z

    lin_force = Inertial_part - Annular_part
    fluid_force = np.flip(lin_force)

    return fluid_force