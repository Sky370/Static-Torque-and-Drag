from scipy.interpolate import interp1d
import pandas as pd 
import numpy as np
import os

THIS_FOLDER = os.path.dirname(os.path.abspath("__file__"))
outputFolderPath = os.path.join(THIS_FOLDER, 'Output')
input_excel_path = os.path.join(THIS_FOLDER, 'Input/NewData.xlsx')

# Import the data
sheet_names = ["BHA", "ADVANCED", "SURVEY", "TOP_DRIVE", "Borehole_Properties", "steady_state_inputs"]
df_dict = pd.read_excel(input_excel_path, sheet_name=sheet_names)

df_BHA = df_dict["BHA"]
df_ADV = df_dict["ADVANCED"]
df_SRV = df_dict["SURVEY"]
df_TOP = df_dict["TOP_DRIVE"]
df_WELL = df_dict["Borehole_Properties"]
df_SS = df_dict["steady_state_inputs"]

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

def survey_mod(df, MD, bf, mass, g):
    x = np.array(df["MD"].values)
    y = np.array(df["INC"].values)
    z = np.array(df["AZI"].values)
    y_interp = interp1d(x, y, kind='linear', fill_value="extrapolate")
    z_interp = interp1d(x, z, kind='linear', fill_value="extrapolate")
    theta_inclination = y_interp(MD)
    theta_azimuth = z_interp(MD)
    # Inclination Angle in rad, not in 'deg'
    Normal_force = bf * mass * g * (np.sin((theta_inclination[:-1] + theta_inclination[1:])/2))
    return theta_inclination, theta_azimuth, Normal_force