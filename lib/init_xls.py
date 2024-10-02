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

def p_drop(rho, mu_p, tao, Q, D_o, D_i, D_w, del_L):
    
    # Annular section
    P_f_an = []
    Area_an = 2.448*(D_w**2-D_o**2)
    v_an = Q/Area_an
    mu_a_an = mu_p + 5*tao*(D_w-D_o)/v_an
    N_re_an = 757*rho*v_an*(D_w-D_o)/mu_a_an
    
    for i in range(len(D_o)):
        if N_re_an[i] > 2100:
            P_f_an.append((rho**0.75*v_an[i]**1.75*mu_p**0.25)/(1396*(D_w[i]-D_o[i])**1.25)) # Turbulent flow case
        else:
            P_f_an.append(((mu_p*v_an[i])/(1000*(D_w[i]-D_o[i])**2) + tao/(200*(D_w[i]-D_o[i])))) # Laminar flow case

    # Inner section
    P_f_in = []
    Area_in = 2.448*D_i**2
    v_in = Q/Area_in
    mu_a_in = mu_p + 6.66*tao*D_i/v_in
    N_re_in = 928*rho*v_in*D_i/mu_a_in
    for i in range(len(D_i)):
        if N_re_in[i] > 2100:
            P_f_in.append((rho**0.75*v_an[i]**1.75*mu_p**0.25)/(1800*D_i[i]**1.25)) # Turbulent flow case
        else:
            P_f_in.append(((mu_p*v_an[i])/(1500*(D_i[i])**2) + tao/(225*D_i[i]))) # Laminar flow case     
    
    return P_f_in, P_f_an, v_in, v_an

def pressure_gradient(rho, mu_p, tao, Q, D_o, D_i, D_w):
    # Conversion to SI
    D_i_new = D_i * 0.0254 # in to m
    D_o_new = D_o * 0.0254 # in to m
    D_w_new = D_w * 0.0254 # in to m
    Q_new = Q / (264.172*60)  # GPM to m3/s
    rho_new = rho* 119.8264273167  # lbm/ft^3 to kg/m^3
    tao_new = tao * 0.4788  # psi to Pa
    mu_p_new = mu_p / 1000  # cp to Pa.s
    m = 1
    
    # Calculation part
    D_hy = D_w_new - D_o_new  # Hydraulic diameter
    Area_an = np.pi / 4 * (D_w_new**2 - D_o_new**2)
    Area_in = np.pi / 4 * D_i_new**2
    v_in = Q_new/Area_in
    v_an = Q_new/Area_an
    Re_a = (12**(1-m) * D_hy**m * v_an**(2-m) * rho_new) / (mu_p_new * ((2*m+1)/(3*m))**m * (1+((2*m+1)/(m+1)) * ((D_hy/(4*v_an)) * (m/(2*m+1)))**m * (tao_new/mu_p_new)))
    # Re_a = rho*v_an*D_hy / (mu_p*(1+(D_hy/(8*v_an))*(tao/mu_p)))

    a = np.where(Re_a < 2100, 24, (np.log10(m) + 3.93) / 50)
    b = np.where(Re_a < 2100, 1, (1.75 - np.log10(m)) / 7)

    fric_an = a / (Re_a**b)
    P_f = 2 * fric_an * rho_new * v_an**2 / D_hy / 22620.60367  # Converting back to Field units (psi/ft)

    return P_f, v_in*3.28084, v_an*3.28084

def pres_evren(rho, mu_p, tao, Q, D_o, D_i, D_w):
    # Conversion to SI
    D_i_new = D_i * 0.0254 # in to m
    D_o_new = D_o * 0.0254 # in to m
    D_w_new = D_w * 0.0254 # in to m
    Q_new = Q / (264.172*60)  # GPM to m3/s
    rho_new = rho* 119.8264273167  # lbm/ft^3 to kg/m^3
    tao_yield = tao * 0.4788  # psi to Pa
    mu_p_new = mu_p / 1000  # cp to Pa.s
    m = 1

    # Calculation part
    D_hy = D_w_new - D_o_new  # Hydraulic diameter
    Area_an = np.pi / 4 * (D_w_new**2 - D_o_new**2)
    Area_in = np.pi / 4 * D_i_new**2
    v_in = Q_new/Area_in
    v_an = Q_new/Area_an

    # Iteration for shear stress
    tao_initial = tao_yield + 12*v_an/D_hy*mu_p_new**m
    tao_new = np.copy(tao_initial)
    tolerance = np.ones_like(tao_initial)  # Array for tolerance

    # Vectorized iteration
    while np.any(tolerance > 1e-4):
        x = tao_yield/tao_initial
        C_c = (1-x)*((m*x/(1+m))+1)
        D_e = 3*m/(2*m+1)*C_c*D_hy
        shear_rate = 12*v_an/D_e
        tao_new = tao_yield + mu_p_new*(shear_rate)**m
        tolerance = np.abs(tao_new - tao_initial)
        tao_initial = np.copy(tao_new)

    N_RE = 12 * rho_new * v_an**2 / tao_new
    N = np.log(tao_new)/np.log(12*v_an/D_hy)
    N_RE_Crit = 3250 - 1150*N

    # Friction factor calculation
    fric_f = np.zeros_like(N_RE)  
    laminar_mask = N_RE < N_RE_Crit
    fric_f[laminar_mask] = 24 / N_RE[laminar_mask]

    # For turbulent flows
    turbulent_mask = ~laminar_mask
    if np.any(turbulent_mask):
        def equation(f_f):
            return (1 / np.sqrt(f_f)) - (4 / (N[turbulent_mask] ** 0.75)) * np.log(N_RE[turbulent_mask] * (f_f ** (1 - N[turbulent_mask]/2))) + (0.4 / (N[turbulent_mask] ** 1.2))

        initial_guess = np.full_like(N_RE[turbulent_mask], 0.01)
        fric_f[turbulent_mask] = fsolve(equation, initial_guess)

    # Pressure gradient calculation
    dPdL = 2 * fric_f * rho_new * v_an**2 / D_hy
    
    return dPdL/22620.40367, v_in*3.28084, v_an*3.28084