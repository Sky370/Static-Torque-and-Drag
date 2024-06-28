from .init_xls import *
import warnings
warnings.filterwarnings('ignore')

class Calculations:
    def __init__(self, lengths):
        # Field units
        self.gravity = 32.1740  # ft/s²
        self.total_length = lengths # ft
        self.elem_length = df_BHA[df_BHA["BHA Type"] == "DP"]["Length (ft)"].iloc[0] # ft
        self.mud_density_ppg = 13  # ppg
        self.visc_p = 28 # cP
        self.tao_y = 10 # lbf/100ft2
        self.bf = (1-self.mud_density_ppg/65.4)

        # Collars:
        self.COLLAR_LEN = df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["Length (ft)"]
        self.COLLAR_MASS = df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["Mass (lbs)"]
        self.COLLAR_OD = np.array(df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["OD (in)"])
        self.COLLAR_ID = np.array(df_BHA[df_BHA["BHA Type"] == "Collar and BHA"]["ID (in)"])

        # Drill pipes:
        self.dp_length = self.total_length - self.COLLAR_LEN.sum()
        self.N_DP, self.L_DP = nearestLength(self.dp_length, self.elem_length)
        self.N_DP = round(self.N_DP)
        self.MASS_DP = df_BHA[df_BHA["BHA Type"] == "DP"]["Mass (lbs)"].iloc[0]*self.N_DP
        self.DP_OD = np.ones(self.N_DP)*df_BHA[df_BHA["BHA Type"] == "DP"]["OD (in)"].iloc[0]
        self.DP_ID = np.ones(self.N_DP)*df_BHA[df_BHA["BHA Type"] == "DP"]["ID (in)"].iloc[0]
        self.L_DP_ARRAY = np.ones(self.N_DP)*self.L_DP
        self.DP_MASS_ARRAY = np.ones(self.N_DP)*(self.MASS_DP/self.N_DP)

        # Concatenation
        self.global_mass_array = np.concatenate([self.DP_MASS_ARRAY, self.COLLAR_MASS])
        self.global_length_array = np.concatenate([self.L_DP_ARRAY, self.COLLAR_LEN])
        self.global_od_array = np.concatenate([self.DP_OD, self.COLLAR_OD])
        self.global_id_array = np.concatenate([self.DP_ID, self.COLLAR_ID])
        self.global_hole_array = np.ones(len(self.global_id_array))*8.75
        self.global_eps = (self.global_hole_array/self.global_od_array) - 1

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