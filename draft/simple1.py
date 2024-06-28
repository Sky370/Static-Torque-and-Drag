import numpy as np 
import pandas as pd 
from scipy.interpolate import interp1d
import os 

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
outputFolderPath = os.path.join(THIS_FOLDER, 'Output')
input_excel_path = os.path.join(THIS_FOLDER, 'Input/Copy of RigData.xlsx')

# Import the data
df_dict = {}
for sheet in ["BHA", "SURVEY"]:
    df_dict[sheet] = pd.read_excel(
        io=input_excel_path,
        sheet_name=sheet,
    )


df_BHA = df_dict["BHA"]
df_SRV = df_dict["SURVEY"]



# Constants and parameters
gravitational_constant = 9.81  # m/s²
mud_density_ppg = 13.3  # ppg
mud_density = mud_density_ppg * 119.8264  # Convert ppg to kg/m³
total_length = 17000*0.3048  # m
friction_factor = 0.15
bf = (1-mud_density_ppg/85.9)
elem_length = float(df_BHA[df_BHA["BHA Type"] == "DP"]["Length (ft)"])*0.3048




# ********************************************************************************
# *****************************    Calculation    ********************************
# ********************************************************************************





# Collars:

COLLAR_LEN = 0.3048*df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["Length (ft)"]
COLLAR_MASS = 0.4535924*df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["Mass (lbs)"]
COLLAR_OD = np.array(0.0254*df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["OD (in)"])
COLLAR_ID = np.array(0.0254*df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["ID (in)"])



# Drill pipes:

dp_length = total_length - COLLAR_LEN.sum()
N_DP, L_DP = nearestLength(dp_length, elem_length)
N_DP = int(N_DP)
MASS_DP = float(0.4535924*df_BHA[df_BHA["BHA Type"] == "DP"]["Mass (lbs)"]*N_DP)
DP_OD = np.ones(N_DP) * float(0.0254*df_BHA[df_BHA["BHA Type"] == "DP"]["OD (in)"])
DP_ID = np.ones(N_DP) * float(0.0254*df_BHA[df_BHA["BHA Type"] == "DP"]["ID (in)"])
L_DP_ARRAY = np.ones(N_DP) * L_DP
DP_MASS_ARRAY = np.ones(N_DP) * (MASS_DP/N_DP)   


global_mass_array = np.concatenate([DP_MASS_ARRAY, COLLAR_MASS])
global_length_array = np.concatenate([L_DP_ARRAY, COLLAR_LEN])
global_od_array = np.concatenate([DP_OD, COLLAR_OD])
global_id_array = np.concatenate([DP_ID, COLLAR_ID])

# Vertical case:

# Calculating weight of the drillstring in air and buoyant force
theta_new, Normal_Force = survey_mod(df_SRV, global_length_array, global_mass_array, bf)
total_weight_air = global_mass_array*gravitational_constant*np.cos(theta_new)
axial_force = total_weight_air*bf


hookload_tripping_in = axial_force - Normal_Force*friction_factor
hookload_tripping_out = axial_force + Normal_Force*friction_factor






