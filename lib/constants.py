from .init_xls import *
import warnings
warnings.filterwarnings('ignore')

class Calculations:
    def __init__(self, depth):
        # SI units
        self.g = ft2m(32.1740)  # ft/s² to m/s²
        self.bit_depth = ft2m(depth) # ft to m
        self.elem_length = ft2m(float(df_ADV[df_ADV["Parameter"] == "Element Length"]["Value"]))        # ft to m
        self.mud_density = ppg2kgm(float(df_WELL[df_WELL["Parameter"] == "Mud Density"]["Value"]))  # ppg to kg/m3
        self.pipe_density = ppg2kgm(float(df_WELL[df_WELL["Parameter"] == "Steel Density"]["Value"]))   # ppg to kg/m3
        self.visc_p = cp2Pas(float(df_WELL[df_WELL["Parameter"] == "Plastic Viscosity"]["Value"]))      # cP to Pa.s 
        self.tao_y = (float(df_WELL[df_WELL["Parameter"] == "Yield Point"]["Value"]))                   # Pa
        self.m = (float(df_WELL[df_WELL["Parameter"] == "m"]["Value"]))                                 
        self.bf = float((1-self.mud_density/self.pipe_density))
        self.mu = float(df_WELL[df_WELL["Parameter"] == "Friction Factor"]["Value"])
        self.Q = GPM2ms(float(df_PUMP[df_PUMP["Parameter"] == "Flow Rate"]["Value"]))                   # GPM  to m3/s

        # Drill pipes:
        self.dp_length = self.bit_depth - ft2m(df_BHA["Total Length (ft)"].iloc[0])
        self.N_DP, self.L_DP = nearestLength(self.dp_length, self.elem_length)
        self.N_DP = round(self.N_DP)        
        self.DP_OD = np.ones(self.N_DP) * df_BHA["OD (in)"].iloc[0]  
        self.DP_ID = np.ones(self.N_DP) * df_BHA["ID (in)"].iloc[0]  
        self.DP_MASS_ARRAY = np.ones(self.N_DP) * df_BHA["Mass (lbs)"].iloc[0]  
        self.L_DP_ARRAY = np.ones(self.N_DP) * self.L_DP
        self.DP_TYPES = np.array(["DP"]*self.N_DP)
        self.TJ_OD = df_BHA["OD Tool Joint (in)"].iloc[0]   # Tool Joint OD, in

        # BHA:
        self.BHA_TYPES = np.array(np.repeat(df_BHA["BHA Type"], df_BHA["Number of Items"]))
        self.BHA_LEN = np.array(np.repeat(df_BHA["Length (ft)"]/df_BHA["Number of Items"] , df_BHA["Number of Items"]))
        self.BHA_MASS = np.array(np.repeat(df_BHA["Mass (lbs)"]/df_BHA["Number of Items"] , df_BHA["Number of Items"]))
        self.BHA_OD = np.array(np.repeat(df_BHA["OD (in)"], df_BHA["Number of Items"]))
        self.BHA_ID = np.array(np.repeat(df_BHA["ID (in)"], df_BHA["Number of Items"])) 

        # Bottom hole
        self.HOLE_OD = in2m(df_ADV[df_ADV["Parameter"] == "Hole Diameter"]["Value"].iloc[0])    # inch to m
        self.HOLE_DEPTH = ft2m(df_ADV[df_ADV["Parameter"] == "Hole Depth"]["Value"].iloc[0])    # ft to m
        self.noe = self.N_DP + len(self.BHA_OD)
        self.HOLE_ARRAY = np.ones(self.noe)*self.HOLE_OD                                        
        self.HOLE_LENGTH = self.HOLE_DEPTH - self.bit_depth
        self.N_HOLE, self.L_HOLE = nearestLength(self.HOLE_LENGTH, self.elem_length)
        self.N_HOLE = round(self.N_HOLE)
        
        if self.HOLE_LENGTH > 0:
            self.OP_HOLE_ARRAY = np.ones(self.N_HOLE) * self.HOLE_OD
            self.global_hole_array = np.concatenate([self.HOLE_ARRAY, self.OP_HOLE_ARRAY])
        elif self.HOLE_DEPTH < self.bit_depth:
            raise ValueError("Hole depth cannot be smaller than bit depth.")
        else:
            self.global_hole_array = self.HOLE_ARRAY
        
        # Concatenation
        self.global_mass_array = lbs2kg(np.concatenate([self.DP_MASS_ARRAY, self.BHA_MASS]))     # Final conversion from lbs to kg
        self.global_length_array = np.concatenate([self.L_DP_ARRAY, ft2m(self.BHA_LEN)])         # ft to m
        self.global_od_array = in2m(np.concatenate([self.DP_OD, self.BHA_OD]))                   # inch to m
        self.global_id_array = in2m(np.concatenate([self.DP_ID, self.BHA_ID]))                   # inch to m
        self.global_types = np.concatenate([self.DP_TYPES, self.BHA_TYPES])
        self.global_eps = (self.HOLE_ARRAY/self.global_od_array) - 1

        # Area, Volume and Density calculation
        self.D_h = self.global_hole_array[:self.noe] - self.global_od_array             # Hydraulic diameter
        self.A_i = np.pi/4*(self.global_id_array)**2                                    # Inner area of the pipe
        self.A_o = np.pi/4*(self.global_od_array)**2                                    # Outer area of the pipe
        self.A_cross = np.pi/4*(self.global_od_array**2-self.global_id_array**2)        # Cross-sectional area of the pipe
        self.A_h = np.pi/4*(self.HOLE_ARRAY**2-self.global_od_array**2)                 # Annular flow area between the wellbore and the pipe
        self.Vol = self.A_cross*self.global_length_array                                # Volume of the pipe
        self.rho = self.global_mass_array/self.Vol                                      # Density of the pipe/bha
        
        # Build and Turn rates calculation
        self.bw_pipe = self.bf*self.global_mass_array/self.global_length_array * self.g
        self.MD = np.insert(np.cumsum(self.global_length_array), 0, 0)
        self.inc, self.azi, self.K = survey_mod_SI(df_SRV, self.MD)
        self.inc_rad = np.round(np.deg2rad(self.inc), 9)
        self.azi_rad = np.round(np.deg2rad(self.azi), 9)
        
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