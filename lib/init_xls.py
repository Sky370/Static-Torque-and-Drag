from scipy.interpolate import interp1d
import numpy as np
import pandas as  pd
from pathlib import Path

base_folder = Path(__file__).resolve().parent.parent
input_excel_path = base_folder / 'Input/NewData.xlsx'
output_folder_path = base_folder / 'Output'

# Import the data
sheet_names = ["BHA", "ADVANCED", "SURVEY", "TOP_DRIVE"]
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

def get_area(s, MD, areas):
    """
    Returns the cross-sectional area for a given depth s using measured depths (MD) and precomputed areas.
    
    Parameters:
    - s (float): The depth at which to find the area.
    - MD (list or array-like): The boundaries of the measured depths for each segment.
    - areas (list or array-like): The precomputed area for each segment defined by MD.
    
    Returns:
    - float: The cross-sectional area at depth s or None if out of range.
    """
    # Ensure MD and areas are numpy arrays for efficient operation
    MD = np.asarray(MD)
    areas = np.asarray(areas)
    
    # Use searchsorted to find the right interval for s
    idx = np.searchsorted(MD, s, side='left')
    
    # Check if the index is within the valid range
    if idx < 0 or idx > len(MD):
        return None
    elif idx == 0:
        return areas[0]
    return areas[idx-1]

def pressure_gradient(rho, mu_p, tao, Q, D_o, D_i, D_w):
    # Conversion to SI
    D_i_new = D_i * 0.0254 # in to m
    D_o_new = D_o * 0.0254 # in to m
    D_w_new = D_w * 0.0254 # in to m
    Q_new = Q / (264.172*60)  # Gallons per minute to cubic meters per second
    rho_new = rho* 119.83  # lbm/ft^3 to kg/m^3
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
