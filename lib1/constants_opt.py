from .init_xls import *
import warnings
warnings.filterwarnings('ignore')

class Calculations:
    def __init__(self, lengths):
        # Field units
        self.gravity = 32.1740  # ft/s²
        self.total_length = lengths # ft
        self.elem_length = df_BHA[df_BHA["BHA Type"] == "DP"]["Length (ft)"].iloc[0] # ft
        self.mud_density_ppg = float(df_ADV[df_ADV["Parameter"] == "Mud Density"]["Value"]) # ppg
        self.visc_p = float(df_ADV[df_ADV["Parameter"] == "Plastic Viscosity"]["Value"]) # cP
        self.tao_y = float(df_ADV[df_ADV["Parameter"] == "Yield Point"]["Value"]) # lbf/100ft2
        self.ff = float(df_ADV[df_ADV["Parameter"] == "Friction Factor"]["Value"])
        self.bf = float((1-self.mud_density_ppg/65.4))

        # Collars:
        self.COLLAR_LEN = df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["Length (ft)"]
        self.COLLAR_MASS = df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["Mass (lbs)"]
        self.COLLAR_OD = np.array(df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["OD (in)"])
        self.COLLAR_ID = np.array(df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["ID (in)"])

        # Drill pipes:
        self.dp_length1 = 4000
        self.dp_length2 = self.total_length - self.COLLAR_LEN.sum() - self.dp_length1
        self.N_DP1, self.L_DP1 = nearestLength(self.dp_length1, self.elem_length)
        self.N_DP2, self.L_DP2 = nearestLength(self.dp_length2, self.elem_length)
        self.N_DP1 = round(self.N_DP1)
        self.N_DP2 = round(self.N_DP2)
        self.MASS_DP1 = df_BHA[df_BHA["BHA Type"] == "DP"]["Mass (lbs)"].iloc[0]*self.N_DP1
        self.MASS_DP2 = df_BHA[df_BHA["BHA Type"] == "DP"]["Mass (lbs)"].iloc[1]*self.N_DP2
        self.DP_OD1 = np.ones(self.N_DP1)*df_BHA[df_BHA["BHA Type"] == "DP"]["OD (in)"].iloc[0]
        self.DP_OD2 = np.ones(self.N_DP2)*df_BHA[df_BHA["BHA Type"] == "DP"]["OD (in)"].iloc[1]
        self.DP_ID1 = np.ones(self.N_DP1)*df_BHA[df_BHA["BHA Type"] == "DP"]["ID (in)"].iloc[0]
        self.DP_ID2 = np.ones(self.N_DP2)*df_BHA[df_BHA["BHA Type"] == "DP"]["ID (in)"].iloc[1]
        self.L_DP_ARRAY1 = np.ones(self.N_DP1)*self.L_DP1
        self.L_DP_ARRAY2 = np.ones(self.N_DP2)*self.L_DP2
        self.DP_MASS_ARRAY1 = np.ones(self.N_DP1)*(self.MASS_DP1/self.N_DP1)
        self.DP_MASS_ARRAY2 = np.ones(self.N_DP2)*(self.MASS_DP2/self.N_DP2)

        # Bottom hole
        self.HOLE_OD = df_ADV[df_ADV["Parameter"] == "Hole Diameter"]["Value"].iloc[0]
        self.HOLE_DEPTH = df_ADV[df_ADV["Parameter"] == "Hole Depth"]["Value"].iloc[0]
        self.HOLE_LENGTH = self.HOLE_DEPTH - self.total_length
        self.N_HOLE, self.L_HOLE = nearestLength(self.HOLE_LENGTH, self.elem_length)
        self.N_HOLE = round(self.N_HOLE)
        self.HOLE_ARRAY = np.ones(self.N_DP1 + self.N_DP2 + len(self.COLLAR_ID))*self.HOLE_OD
        
        if self.HOLE_LENGTH > 0:
            self.OP_HOLE_ARRAY = np.ones(self.N_HOLE) * self.HOLE_OD
            self.global_hole_array = np.concatenate([self.HOLE_ARRAY, self.OP_HOLE_ARRAY])
        elif self.HOLE_DEPTH < self.bit_depth:
            raise ValueError("Hole depth cannot be smaller than bit depth.")
        else:
            self.global_hole_array = self.HOLE_ARRAY

        # Concatenation
        self.global_mass_array = np.concatenate([self.DP_MASS_ARRAY1, self.DP_MASS_ARRAY2, self.COLLAR_MASS])
        self.global_length_array = np.concatenate([self.L_DP_ARRAY1, self.L_DP_ARRAY2, self.COLLAR_LEN])
        self.global_od_array = np.concatenate([self.DP_OD1, self.DP_OD2, self.COLLAR_OD])
        self.global_id_array = np.concatenate([self.DP_ID1, self.DP_ID2, self.COLLAR_ID])
        self.global_eps = (self.HOLE_ARRAY/self.global_od_array) - 1

        # Build and Turn rates calculation
        self.bw_pipe = self.bf*self.global_mass_array/self.global_length_array
        self.MD = np.insert(np.cumsum(self.global_length_array), 0, 0)
        self.inc, self.azi = survey_mod(df_SRV, self.MD)
        self.inc_rad = np.deg2rad(self.inc)
        self.azi_rad = np.deg2rad(self.azi)

        # Minimum Curvature Method optimized further with np.diff
        cos_delta_azi = np.cos(np.diff(self.azi_rad))
        self.md_diff = np.diff(self.MD)

        self.betta = np.arccos(np.sin(self.inc_rad[:-1])*np.sin(self.inc_rad[1:])*cos_delta_azi + np.cos(self.inc_rad[:-1])*np.cos(self.inc_rad[1:]))
        self.betta = np.nan_to_num(self.betta, nan=0.0)
        self.RF = np.where(self.betta == 0, self.md_diff/2, self.md_diff/self.betta*np.tan(self.betta/2))
        self.DLS = self.betta/self.md_diff

        # Displacement calculations optimized with numpy
        self.del_x = self.RF*(np.sin(self.inc_rad[:-1])*np.cos(self.azi_rad[:-1]) + np.sin(self.inc_rad[1:])*np.cos(self.azi_rad[1:]))
        self.del_y = self.RF*(np.sin(self.inc_rad[:-1])*np.sin(self.azi_rad[:-1]) + np.sin(self.inc_rad[1:])*np.sin(self.azi_rad[1:]))
        self.del_z = self.RF*(np.cos(self.inc_rad[:-1]) + np.cos(self.inc_rad[1:]))
        self.MD_s = (self.MD[:-1] + self.MD[1:])/2

        # # B, T, and inc_ave calculations
        # self.B = np.diff(self.inc)/self.md_diff
        # self.T = np.diff(self.azi)/self.md_diff
        # inc_ave = (self.inc_rad[:-1] + self.inc_rad[1:]) / 2

        # # Tension, Normal force, and Side force calculations optimized
        # self.t_z = np.cos(inc_ave)
        # self.n_z = -1 / np.degrees(self.DLS) * np.sin(inc_ave) * self.B
        # self.b_z = 1 / np.degrees(self.DLS) * np.sin(inc_ave)**2 * self.T
        
        # Tension, Normal force, and Side force calculations optimized
        
        self.t_z_min_val = np.cos(self.inc_rad[:-1])*np.cos(self.DLS*(self.MD_s - self.MD[:-1])) + \
            ((np.cos(self.inc_rad[1:]) - np.cos(self.inc_rad[:-1])*np.cos(self.betta))/np.sin(self.betta))*np.sin(self.DLS*(self.MD_s - self.MD[:-1]))
        self.t_z = np.where(np.isnan(self.t_z_min_val), np.cos((self.inc_rad[1:] + self.inc_rad[:-1]) / 2), self.t_z_min_val)
        self.n_z = -np.cos(self.inc_rad[:-1])*np.sin(self.DLS*(self.MD_s - self.MD[:-1])) + \
            ((np.cos(self.inc_rad[1:]) - np.cos(self.inc_rad[:-1])*np.cos(self.betta))/np.sin(self.betta))*np.cos(self.DLS*(self.MD_s - self.MD[:-1]))
        self.n_z = np.nan_to_num(self.n_z)
        self.b_z = np.sin(self.inc_rad[:-1])*np.sin(self.inc_rad[1:])*np.sin(self.azi_rad[1:] - self.azi_rad[:-1])/np.sin(self.betta)
        self.b_z = np.nan_to_num(self.b_z)