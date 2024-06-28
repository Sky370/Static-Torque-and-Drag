# import math

# def calculate_hookload(mud_density_ppg, segment_mass, OD_in, ID_in, total_length, segment_length, inclinations, friction_factor):
#     """
#     Calculate the hookload for tripping in and out processes in a J-type wellbore.

#     Parameters:
#     mud_density_ppg (float): Density of the mud in pounds per gallon.
#     segment_mass (float): Mass of each drillstring segment in pounds.
#     OD_in (float): Outer diameter of the drillstring in inches.
#     ID_in (float): Inner diameter of the drillstring in inches.
#     total_length (float): Total measured depth of the wellbore in feet.
#     segment_length (float): Length of each drillstring segment in feet.
#     inclinations (list of floats): List of inclination angles in degrees for each wellbore segment.
#     friction_factor (float): Friction factor for the wellbore.

#     Returns:
#     float: Total frictional force.
#     """
#     gravitational_constant = 32.2  # ft/s²
#     mud_density = mud_density_ppg * 0.051948051948052  # Convert ppg to lbm/ft³
#     OD_ft = OD_in / 12  # Convert inches to feet
#     ID_ft = ID_in / 12  # Convert inches to feet
#     num_segments = total_length / segment_length
#     total_weight_air = num_segments * segment_mass * gravitational_constant
#     volume_drillstring = (math.pi / 4) * (OD_ft**2 - ID_ft**2) * total_length
#     buoyant_force = volume_drillstring * mud_density * gravitational_constant
#     total_weight_fluid = total_weight_air - buoyant_force
#     total_frictional_force = 0
#     segment_lengths = [100] * (len(inclinations) - 1)  # Assuming uniform spacing for simplicity

#     for i in range(len(segment_lengths)):
#         normal_force = total_weight_fluid * math.cos(math.radians(inclinations[i]))
#         frictional_force_segment = normal_force * friction_factor * segment_lengths[i]
#         total_frictional_force += frictional_force_segment

#     return total_frictional_force

# # Example Usage
# mud_density_ppg = 13
# segment_mass = 800
# OD_in = 7
# ID_in = 3
# total_length = 1000
# segment_length = 30
# inclinations = [0, 0, 0, 0, 15, 30, 45, 60, 80, 90, 90]
# friction_factor = 0.15

# total_frictional_force = calculate_hookload(mud_density_ppg, segment_mass, OD_in, ID_in, total_length, segment_length, inclinations, friction_factor)
# print("Total Frictional Force:", total_frictional_force)
 

 # Re-importing math module and redefining all necessary parameters and variables
import numpy as np 
import pandas as pd 
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

# Data for 13000 ft
sv13 = df_SRV[:592]
num_pipes13 = int((sv13.iloc[-1,0] - 99.3)/31.5)
df_BHA = pd.concat([pd.concat([df_BHA.iloc[0:1]]*num_pipes13, ignore_index=True),df_BHA], ignore_index=True)


# Constants and parameters
gravitational_constant = 32.2  # ft/s²
mud_density_ppg = 13.3  # ppg
mud_density = mud_density_ppg*7.48
# mud_density = mud_density_ppg * 0.051948051948052  # Convert ppg to lbm/ft³
segment_mass = 850  # lbs
total_length = 13000  # ft
segment_length = 31.5  # ft
OD_ft = df_BHA["OD (in)"]/12 # Convert in to ft
ID_ft = df_BHA["ID (in)"]/12 # Convert in to ft
A_i = np.pi*ID_ft**2/4
A_o = np.pi*OD_ft**2/4
D_h = 8.25/12 - OD_ft
A_h = np.pi*D_h**2/4
friction_factor = 0.15
num_segments = total_length / segment_length
inclinations = sv13.iloc[:,1]
segment_lengths = [22] * (len(inclinations) - 1)  # Assuming uniform spacing

# Calculating weight of the drillstring in air and buoyant force
total_weight_air = num_segments * segment_mass
volume_drillstring = (np.pi / 4) * (OD_ft**2 - ID_ft**2) * segment_length
buoyant_force = volume_drillstring * mud_density * gravitational_constant
total_weight_fluid = total_weight_air - buoyant_force.sum()

# Calculating total frictional force
bend_contact_length = 3  # ft, length of the contact at the bend
increased_friction_factor = 2 * friction_factor  # Assuming increased friction at the bend
total_frictional_force = 0
total_fr_f = 0

for i in range(len(segment_lengths)):
    bf = (1-mud_density_ppg/65.4)
    new_norm_f = bf*segment_mass*gravitational_constant*np.sin(np.radians(inclinations[i]))
    fr_force = new_norm_f*friction_factor*segment_lengths[i]
    normal_force = total_weight_fluid * np.cos(np.radians(inclinations[i]))
    frictional_force_segment = normal_force * friction_factor * segment_lengths[i]

    total_frictional_force += frictional_force_segment
    total_fr_f += fr_force



hookload_tripping_in = total_weight_air - buoyant_force - total_frictional_force
hookload_tripping_out = total_weight_air - buoyant_force + total_frictional_force

hookload_tripping_in, hookload_tripping_out

