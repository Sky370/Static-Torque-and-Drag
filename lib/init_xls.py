import pandas as pd 
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

THIS_FOLDER = os.path.dirname(os.path.abspath("__file__"))
outputFolderPath = os.path.join(THIS_FOLDER, 'Output')
input_excel_path = os.path.join(THIS_FOLDER, 'Input/NewData.xlsx')

# Import the data
sheet_names = ["PUMP", "BHA", "ADVANCED", "SURVEY", "Borehole_Properties"]
df_dict = pd.read_excel(input_excel_path, sheet_name=sheet_names)

df_PUMP = df_dict["PUMP"]
df_BHA = df_dict["BHA"]
df_ADV = df_dict["ADVANCED"]
df_SRV = df_dict["SURVEY"]
df_WELL = df_dict["Borehole_Properties"]

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

def survey_mod_SI(df, MD):
    x = ft2m(np.array(df["MD"].values))
    y = np.array(df["INC"].values)
    z = np.array(df["AZI"].values)
    dls = np.array(df["DLS"].values)
    y_interp = interp1d(x, y, kind='linear', fill_value="extrapolate")
    z_interp = interp1d(x, z, kind='linear', fill_value="extrapolate")
    k_interp = interp1d(x, dls, kind='linear', fill_value="extrapolate")
    theta_inclination = y_interp(MD)
    theta_azimuth = z_interp(MD)
    DLS = k_interp(MD) # deg/m
    # Inclination Angle in rad, not in 'deg'
    # Normal_force = bf * mass * g * np.sin(np.deg2rad((theta_inclination[:-1] + theta_inclination[1:])/2))
    return theta_inclination, theta_azimuth, DLS

def survey_mod_IMPERIAL(df, MD):
    x = (np.array(df["MD"].values))
    y = np.array(df["INC"].values)
    z = np.array(df["AZI"].values)
    dls = (np.array(df["DLS"].values))
    y_interp = interp1d(x, y, kind='linear', fill_value="extrapolate")
    z_interp = interp1d(x, z, kind='linear', fill_value="extrapolate")
    k_interp = interp1d(x, dls, kind='linear', fill_value="extrapolate")
    theta_inclination = y_interp(MD)
    theta_azimuth = z_interp(MD)
    DLS = k_interp(MD) # deg/m
    # Inclination Angle in rad, not in 'deg'
    # Normal_force = bf * mass * g * np.sin(np.deg2rad((theta_inclination[:-1] + theta_inclination[1:])/2))
    return theta_inclination, theta_azimuth, DLS

# New Pressure Loss Calculation (more accurate)
def pres_calc(rho, m, K, tao, Q, clc):
    v_in = Q / clc.A_i
    v_an = Q / clc.A_h

    def shear_stress_iteration(v, D, is_annular=True):
        # Iterative calculation for shear stress
        tao_initial = tao + K * ((12 if is_annular else 8) * v / D)**m
        tao_new = np.copy(tao_initial)
        tolerance = np.ones_like(tao_initial)

        while np.any(tolerance > 1e-4):
            x = tao / tao_initial
            if is_annular:
                C_c = (1-x)*((m*x/(1+m))+1)
                D_e = 3*m/(2*m+1)*C_c*D
            else:
                C_c = (1-x)*(2*(m*x)**2/((1+2*m)*(1+m))+2*m*x/(1+2*m)+1)
                D_e = 4*m/(3*m+1)*C_c*D

            shear_rate = (12 if is_annular else 8) * v / D_e
            tao_new = tao + K * shear_rate**m
            tolerance = np.abs(tao_new - tao_initial)
            tao_initial = np.copy(tao_new)

        return tao_new

    # Iteration for shear stress in annulus and inner sections
    tao_new_an = shear_stress_iteration(v_an, clc.D_h, is_annular=True)
    tao_new_in = shear_stress_iteration(v_in, clc.global_id_array, is_annular=False)

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

    N_RE_an = 12 * rho * v_an**2 / tao_new_an
    N_RE_in = 8 * rho * v_in**2 / tao_new_in
    N_an = np.log(tao_new_an) / np.log(12 * v_an / clc.D_h)
    N_in = np.log(tao_new_in) / np.log(8 * v_in / clc.global_id_array)

    fric_f_an = calculate_friction_factor(N_RE_an, N_an, is_annular=True)
    fric_f_in = calculate_friction_factor(N_RE_in, N_in, is_annular=False)

    # Pressure gradient calculation
    dPdL_an = 2 * fric_f_an * rho * v_an**2 / clc.D_h
    dPdL_in = 2 * fric_f_in * rho * v_in**2 / clc.global_id_array

    return dPdL_an, dPdL_in, v_an, v_in

ft2m = lambda ft: ft * 0.3048
in2m = lambda inch: inch * 0.0254
ft_min2ms = lambda f: f * 0.3048 / 60
ft_hr2ms = lambda f: f * 0.3048 / 3600
ppg2kgm = lambda rho: rho * 1.1983e+2
ppg2lbft3 = lambda ppg: ppg*7.4805e+0
psi2Pa = lambda psi: psi * 6894.76
psi2psf = lambda psi: psi * 144
pa2psi = lambda psi: psi / 6894.76
psf2Pa = lambda psf: psf * 47.880208
cp2Pas = lambda viscosity: viscosity/1000
cp2lbfft2 = lambda cp: cp*2.0885e-5
lbfft2Pas = lambda viscosity: 47.88*viscosity  # (lbf/ft^2).s to Pa.s
lbf100ft2Pa = lambda lbf: 0.4788*lbf 
GPM2ms = lambda flowrate: flowrate * 6.309e-5
lbs2kg = lambda lbs: lbs * 4.536e-1
lbf2N = lambda lbf: lbf * 4.4482e+0
rpm2rad_s = lambda rpm: rpm * 2*np.pi/60
# Unit conversion [metric-imperial]
m2ft = lambda m: m / 0.3048
m2in = lambda m: m / 0.0254
N2lbf = lambda lbf: lbf / 4.44822
Nm2lbfft = lambda Nm: Nm * 0.73756215