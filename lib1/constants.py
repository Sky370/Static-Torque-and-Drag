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
        self.dp_length1 = 4000
        self.dp_length2 = self.total_length - self.COLLAR_LEN.sum() - self.dp_length1
        self.N_DP1, self.L_DP1 = nearestLength(self.dp_length1, self.elem_length)
        self.N_DP2, self.L_DP2 = nearestLength(self.dp_length2, self.elem_length)
        self.N_DP1 = round(self.N_DP1)
        self.N_DP2 = round(self.N_DP2)
        self.MASS_DP1 = df_BHA[df_BHA["BHA Type"] == "DP"]["Mass (lbs)"].iloc[0]*self.N_DP1
        self.MASS_DP2 = df_BHA[df_BHA["BHA Type"] == "DP"]["Mass (lbs)"].iloc[1]*self.N_DP2
        self.DP_OD1 = np.ones(self.N_DP1) * df_BHA[df_BHA["BHA Type"] == "DP"]["OD (in)"].iloc[0]
        self.DP_OD2 = np.ones(self.N_DP2) * df_BHA[df_BHA["BHA Type"] == "DP"]["OD (in)"].iloc[1]
        self.DP_ID1 = np.ones(self.N_DP1) * df_BHA[df_BHA["BHA Type"] == "DP"]["ID (in)"].iloc[0]
        self.DP_ID2 = np.ones(self.N_DP2) * df_BHA[df_BHA["BHA Type"] == "DP"]["ID (in)"].iloc[1]
        self.L_DP_ARRAY1 = np.ones(self.N_DP1) * self.L_DP1
        self.L_DP_ARRAY2 = np.ones(self.N_DP2) * self.L_DP2
        self.DP_MASS_ARRAY1 = np.ones(self.N_DP1) * (self.MASS_DP1/self.N_DP1)
        self.DP_MASS_ARRAY2 = np.ones(self.N_DP2) * (self.MASS_DP2/self.N_DP2)

        # Concatenation
        self.global_mass_array = np.concatenate([self.DP_MASS_ARRAY1, self.DP_MASS_ARRAY2, self.COLLAR_MASS])
        self.global_length_array = np.concatenate([self.L_DP_ARRAY1, self.L_DP_ARRAY2, self.COLLAR_LEN])
        self.global_od_array = np.concatenate([self.DP_OD1, self.DP_OD2, self.COLLAR_OD])
        self.global_id_array = np.concatenate([self.DP_ID1, self.DP_ID2, self.COLLAR_ID])
        self.global_hole_array = np.ones(len(self.global_id_array))*8.75
        self.global_eps = (self.global_hole_array/self.global_od_array) - 1
        
        # Build and Turn rates calculation
        self.bw_pipe = self.bf*self.global_mass_array/self.global_length_array
        self.MD = np.cumsum(self.global_length_array)
        self.MD = np.insert(self.MD, 0, 0)
        self.inc, self.azi = survey_mod(df_SRV, self.MD)
        self.inc_rad = np.deg2rad(self.inc)
        self.azi_rad = np.deg2rad(self.azi)

        # Minimum Curvature Method
        self.betta = [np.arccos(np.sin(self.inc_rad[i])*np.sin(self.inc_rad[i+1])*np.cos(self.azi_rad[i+1]-self.azi_rad[i]) \
                                + np.cos(self.inc_rad[i])*np.cos(self.inc_rad[i+1])) for i in range(len(self.MD)-1)]
        self.betta = np.nan_to_num(self.betta)
        self.RF_val = [(self.MD[i+1]-self.MD[i])/self.betta[i]*np.tan(self.betta[i]/2) for i in range(len(self.betta))]
        self.RF_val = np.nan_to_num(self.RF_val)
        self.RF = [(self.MD[i+1] - self.MD[i])/2 if val == 0 else val for i, val in enumerate(self.RF_val)]
        self.DLS = [np.degrees(self.betta)[i]/(self.MD[i+1] - self.MD[i]) for i in range(len(self.betta))]

        self.del_x = [self.RF[i]*(np.sin(self.inc_rad[i])*np.cos(self.azi_rad[i]) + np.sin(self.inc_rad[i+1])*np.cos(self.azi_rad[i+1])) for i in range(len(self.RF))]
        self.del_y = [self.RF[i]*(np.sin(self.inc_rad[i])*np.sin(self.azi_rad[i]) + np.sin(self.inc_rad[i+1])*np.sin(self.azi_rad[i+1])) for i in range(len(self.RF))]
        self.del_z = [self.RF[i]*(np.cos(self.inc_rad[i])+np.cos(self.inc_rad[i+1])) for i in range(len(self.RF))]

        self.MD_s = [(self.MD[i]+self.MD[i+1])/2 for i in range(len(self.MD)-1)]

        self.t_z_min_val = [np.cos(self.inc_rad[i])*np.cos(np.deg2rad(self.DLS[i]*(self.MD_s[i]-self.MD[i]))) \
                        +((np.cos(self.inc_rad[i+1])-np.cos(self.inc_rad[i])*np.cos(self.betta[i])) \
                          /np.sin(self.betta[i]))*np.sin(np.deg2rad(self.DLS[i]*(self.MD_s[i]-self.MD[i]))) for i in range(len(self.MD)-1)]
        self.t_z_min = [np.cos((self.inc_rad[i+1] + self.inc_rad[i])/2) if np.isnan(val) else val for i, val in enumerate(self.t_z_min_val)]
        self.n_z_min = [-np.cos(self.inc_rad[i])*np.sin(np.deg2rad(self.DLS[i]*(self.MD_s[i]-self.MD[i]))) \
                        +((np.cos(self.inc_rad[i+1])-np.cos(self.inc_rad[i])*np.cos(self.betta[i])) \
                          /np.sin(self.betta[i]))*np.cos(np.deg2rad(self.DLS[i]*(self.MD_s[i]-self.MD[i]))) for i in range(len(self.MD)-1)]
        self.n_z_min = np.nan_to_num(self.n_z_min)
        self.b_z_min = [np.sin(self.inc_rad[i])*np.sin(self.inc_rad[i+1])*np.sin(self.azi_rad[i+1]-self.azi_rad[i])/np.sin(self.betta[i]) for i in range(len(self.MD)-1)]
        self.b_z_min = np.nan_to_num(self.b_z_min)
