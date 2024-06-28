import pandas as pd 
import os
import numpy as np
from scipy.interpolate import interp1d

THIS_FOLDER = os.path.dirname(os.path.abspath("__file__"))
outputFolderPath = os.path.join(THIS_FOLDER, 'Output')
input_excel_path = os.path.join(THIS_FOLDER, 'Input/NewData.xlsx')

# Import the data
df_dict = {}
for sheet in ["BHA", "SURVEY"]:
    df_dict[sheet] = pd.read_excel(
        io=input_excel_path,
        sheet_name=sheet,
    )


df_BHA = df_dict["BHA"]
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
    Area_in = 2.448*D_i
    v_in = Q/Area_in
    mu_a_in = mu_p + 6.66*tao*D_i/v_in
    N_re_in = 928*rho*v_in*D_i/mu_a_in
    for i in range(len(D_i)):
        if N_re_in[i] > 2100:
            P_f_in.append((rho**0.75*v_an[i]**1.75*mu_p**0.25)/(1800*D_i[i]**1.25)) # Turbulent flow case
        else:
            P_f_in.append(((mu_p*v_an[i])/(1500*(D_i[i])**2) + tao/(225*D_i[i]))) # Laminar flow case     
    
    return P_f_in, P_f_an, v_in, v_an