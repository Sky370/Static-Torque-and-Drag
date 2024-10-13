from .constants import *
from scipy.optimize import fsolve

def pres_calc(rho, mu_p, tao, Q, D_o, D_i, D_w):
    # Conversion to SI
    D_i_new, D_o_new, D_w_new = D_i*0.0254, D_o*0.0254, D_w*0.0254  # in to m
    Q_new = Q / (264.172 * 60)                                      # GPM to m3/s
    rho_new = rho * 119.8264273167                                  # lbm/ft^3 to kg/m^3
    tao_yield = tao * 0.4788                                        # psi to Pa
    mu_p_new = mu_p / 1000                                          # cp to Pa.s
    m = 0.65

    # Calculation part
    D_hy = D_w_new - D_o_new       # Hydraulic diameter
    v_in = Q_new / Calculations.A_i
    v_an = Q_new / Calculations.A_h

    def shear_stress_iteration(v, D, is_annular=True):
        # Iterative calculation for shear stress
        tao_initial = tao_yield + (12 if is_annular else 8) * v / D * mu_p_new**m
        tao_new = np.copy(tao_initial)
        tolerance = np.ones_like(tao_initial)

        while np.any(tolerance > 1e-4):
            x = tao_yield / tao_initial
            if is_annular:
                C_c = (1-x)*((m*x/(1+m))+1)
                D_e = 3*m/(2*m+1)*C_c*D
            else:
                C_c = (1-x)*(2*(m*x)**2/((1+2*m)*(1+m))+2*m*x/(1+2*m)+1)
                D_e = 4*m/(3*m+1)*C_c*D

            shear_rate = (12 if is_annular else 8) * v / D_e
            tao_new = tao_yield + mu_p_new * shear_rate**m
            tolerance = np.abs(tao_new - tao_initial)
            tao_initial = np.copy(tao_new)

        return tao_new

    # Iteration for shear stress in annulus and inner sections
    tao_new_an = shear_stress_iteration(v_an, D_hy, is_annular=True)
    tao_new_in = shear_stress_iteration(v_in, D_i_new, is_annular=False)

    # Reynolds number and friction factor calculations
    def calculate_friction_factor(N_RE, N, is_annular=True):
        N_RE_Crit = 3250 - 1150 * N
        laminar_mask = N_RE < N_RE_Crit
        fric_f = np.zeros_like(N_RE)
        fric_f[laminar_mask] = 24 / N_RE[laminar_mask] if is_annular else 16 / N_RE[laminar_mask]

        turbulent_mask = ~laminar_mask
        if np.any(turbulent_mask):
            def equation(f_f):
                return (1 / np.sqrt(f_f)) - (4 / (N[turbulent_mask] ** 0.75)) * \
                    np.log(N_RE[turbulent_mask] * (f_f ** (1 - N[turbulent_mask] / 2))) + (0.4 / (N[turbulent_mask] ** 1.2))

            initial_guess = np.full_like(N_RE[turbulent_mask], 0.01)
            fric_f[turbulent_mask] = fsolve(equation, initial_guess)

        return fric_f

    N_RE_an = 12 * rho_new * v_an**2 / tao_new_an
    N_RE_in = 8 * rho_new * v_in**2 / tao_new_in
    N_an = np.log(tao_new_an) / np.log(12 * v_an / D_hy)
    N_in = np.log(tao_new_in) / np.log(8 * v_in / D_i_new)

    fric_f_an = calculate_friction_factor(N_RE_an, N_an, is_annular=True)
    fric_f_in = calculate_friction_factor(N_RE_in, N_in, is_annular=False)

    # Pressure gradient calculation
    dPdL_an = 2 * fric_f_an * rho_new * v_an**2 / D_hy
    dPdL_in = 2 * fric_f_in * rho_new * v_in**2 / D_i_new

    return dPdL_an / 22620.40367, dPdL_in / 22620.40367, v_in * 3.28084, v_an * 3.28084

def circl(flowrate, length):
    clc = Calculations(lengths=length)
    P_drop_outer, P_drop_inner, v_i, v_o = pres_calc(
        Q=flowrate,
        rho=clc.mud_density_ppg, 
        mu_p=clc.visc_p, 
        tao=clc.tao_y,
        D_o=clc.global_od_array,
        D_i=clc.global_id_array, 
        D_w=clc.global_hole_array
    )
    # Area_surf = np.pi*clc.global_od_array*clc.global_length_array
    C_tao_0 = (8 + 4*clc.global_eps + clc.global_eps**2)/8

    Static_part = (12/231)*(clc.A_i-clc.A_o)*clc.mud_density_ppg*clc.gravity
    Annular_part = (1+C_tao_0)*P_drop_outer*clc.A_h
    # Annular_part = (1+C_tao_0)*P_drop_outer*Area_surf
    Inertial_part = (12/231)*clc.mud_density_ppg*(v_o**2*clc.A_h - v_i**2*clc.A_i)*clc.DLS

    lin_force = Inertial_part - Annular_part
    # fluid_force = np.flip(lin_force)
    fluid_force = lin_force

    return fluid_force