import pandas as pd 
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

THIS_FOLDER = os.path.dirname(os.path.abspath("__file__"))
outputFolderPath = os.path.join(THIS_FOLDER, 'Output')
input_excel_path = os.path.join(THIS_FOLDER, 'Input/NewData.xlsx')

# Import the data
sheet_names = ["BHA", "ADVANCED", "SURVEY"]
df_dict = pd.read_excel(input_excel_path, sheet_name=sheet_names)

df_BHA = df_dict["BHA"]
df_ADV = df_dict["ADVANCED"]
df_SRV = df_dict["SURVEY"]
df_TOP = df_dict["TOP_DRIVE"]

def nearestLength(a, b):
    """
    Takes two numbers a and b and tries to break a into 'n' integer parts with
    length close to b
    Parameters
    ----------
    a : `float`
        time in seconds
    b : `float`
        time in seconds

    Returns
    -------
    num , length : `tuple` of `float`

    """
    if a <= b:
        return 1, a
    ceil_var = 1.0 * np.ceil(a / b)
    floor_var = 1.0 * np.floor(a / b)
    dict_var = {}
    dict_var[abs(b - (a / ceil_var))] = a / ceil_var
    dict_var[abs(b - (a / floor_var))] = a / floor_var
    val = min(abs(b - (a / ceil_var)), abs(b - (a / floor_var)))
    dict_var[val]
    return a / dict_var[val], dict_var[val]

def survey_mod(df, MD):
    x = np.array(df["MD"].values)
    y = np.array(df["INC"].values)
    z = np.array(df["AZI"].values)
    y_interp = interp1d(x, y, kind='linear', fill_value="extrapolate")
    z_interp = interp1d(x, z, kind='linear', fill_value="extrapolate")
    theta_inclination = y_interp(MD)
    theta_azimuth = z_interp(MD)
    # Inclination Angle in rad, not in 'deg'
    return theta_inclination, theta_azimuth

def pres_calc(rho, mu_p, tao, Q, D_o, D_i, D_w):
    # Conversion to SI
    D_i_new, D_o_new, D_w_new = D_i*0.0254, D_o*0.0254, D_w*0.0254  # in to m
    Q_new = Q / (264.172 * 60)                                      # GPM to m3/s
    rho_new = rho * 119.8264273167                                  # lbm/ft^3 to kg/m^3
    tao_yield = tao * 0.4788                                        # psi to Pa
    mu_p_new = mu_p / 1000                                          # cp to Pa.s
    m = 0.65

    # Calculation part
    D_hy = D_w_new - D_o_new                            # Hydraulic diameter
    Area_an = np.pi / 4 * (D_w_new**2 - D_o_new**2)     # Annular area
    Area_in = np.pi / 4 * D_i_new**2                    # Inner area
    v_in = Q_new / Area_in
    v_an = Q_new / Area_an

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