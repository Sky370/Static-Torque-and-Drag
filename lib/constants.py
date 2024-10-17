from .init_xls import *
import scipy.sparse as sps
import warnings
warnings.filterwarnings('ignore')

class Calculations:
    def __init__(self, depth):
        # Field units
        self.g = ft2m(32.1740)  # ft/s² to m/s²
        self.bit_depth = ft2m(depth) # ft to m
        self.elem_length = ft2m(df_BHA[df_BHA["BHA Type"] == "DP"]["Length (ft)"].iloc[0])              # ft to m
        self.mud_density_ppg = ppg2kgm(float(df_WELL[df_WELL["Parameter"] == "Mud Density"]["Value"]))  # ppg to kg/m3
        self.pipe_density = ppg2kgm(float(df_WELL[df_WELL["Parameter"] == "Steel Density"]["Value"]))   # ppg to kg/m3
        self.visc_p = cp2Pas(float(df_WELL[df_WELL["Parameter"] == "Plastic Viscosity"]["Value"]))      # cP to Pa.s 
        self.tao_y = psi2Pa(float(df_WELL[df_WELL["Parameter"] == "Yield Point"]["Value"]))             # lbf/100ft2 to Pa
        self.bf = float((1-self.mud_density_ppg/self.pipe_density))
        self.mu_s = float(df_WELL[df_WELL["Parameter"] == "Static Friction Factor"]["Value"])           # Static Fric. 
        self.mu_d = float(df_WELL[df_WELL["Parameter"] == "Dynamic Friction Factor"]["Value"])          # Dynamic Fric.
        self.E = float(df_ADV[df_ADV["Parameter"] == "Young Modulus"]["Value"])                         # Young Modulus
        self.G = float(df_ADV[df_ADV["Parameter"] == "Bulk Modulus"]["Value"])                          # Shear Modulus
        self.ccs = float(df_ADV[df_ADV["Parameter"] == "CCS"]["Value"])                                 # ksi
        self.v_cs = float(df_WELL[df_WELL["Parameter"] == "Stribeck Critical Velocity"]["Value"])       # m/s
        self.CT_BOREHOLE = float(df_WELL[df_WELL["Parameter"] == "Torsional Drag Coefficient"]["Value"])    # N sec/m

        # Time intervals definition
        self.a1 = float(df_TOP[df_TOP["Parameter"] == "a1"]["Value"])
        self.a2 = float(df_TOP[df_TOP["Parameter"] == "a2"]["Value"])
        self.a3 = float(df_TOP[df_TOP["Parameter"] == "a3"]["Value"])
        self.a4 = float(df_TOP[df_TOP["Parameter"] == "a4"]["Value"])
        self.a5 = float(df_TOP[df_TOP["Parameter"] == "a5"]["Value"])
        self.a6 = float(df_TOP[df_TOP["Parameter"] == "a6"]["Value"])
        self.b1 = float(df_TOP[df_TOP["Parameter"] == "b1"]["Value"])
        self.b2 = float(df_TOP[df_TOP["Parameter"] == "b2"]["Value"])
        self.b3 = float(df_TOP[df_TOP["Parameter"] == "b3"]["Value"])
        self.b4 = float(df_TOP[df_TOP["Parameter"] == "b4"]["Value"])
        self.b5 = float(df_TOP[df_TOP["Parameter"] == "b5"]["Value"])
        self.b6 = float(df_TOP[df_TOP["Parameter"] == "b6"]["Value"])
        
        # Velocity definition
        self.v1 = ft_min2ms(float(df_TOP[df_TOP["Parameter"] == "Top Drive Axial Velocity Magnitude 1 (ft/min)"]["Value"]))     # ft/min to m/s
        self.v2 = ft_min2ms(float(df_TOP[df_TOP["Parameter"] == "Top Drive Axial Velocity Magnitude 2 (ft/min)"]["Value"]))     # ft/min to m/s
        self.rpm1 = float(df_TOP[df_TOP["Parameter"] == "Top Drive RPM Magnitude 1 (RPM)"]["Value"])                            # rev/min
        self.rpm2 = float(df_TOP[df_TOP["Parameter"] == "Top Drive RPM Magnitude 2 (RPM)"]["Value"])                            # rev/min

        # Time array
        self.time_val = float(df_ADV[df_ADV["Parameter"] == "Run Time"]["Value"])   # seconds

        # Extra
        # --------- Load steady state inputs -------- #
        WOB_SS = float(df_SS[df_SS["Parameter"] == "WOB initial"]["Value"])                 # lbs
        ROP_SS = float(df_SS[df_SS["Parameter"] == "ROP steady state"]["Value"]) / 3600     # m/hr to m/s
        RPM_SS = float(df_SS[df_SS["Parameter"] == "RPM steady state"]["Value"]) / 60       # rev/min to rev/s  

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
        self.HOLE_ARRAY = np.ones(self.noe)*self.HOLE_OD                                        # inch to m
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
        
        # Area calculation
        self.A_i = np.pi/4*(self.global_id_array)**2                                    # Inner area of the pipe
        self.A_o = np.pi/4*(self.global_od_array)**2                                    # Outer area of the pipe
        self.A_cross = np.pi/4*(self.global_od_array**2-self.global_id_array**2)        # Cross-sectional area of the pipe
        self.A_h = np.pi/4*(self.HOLE_ARRAY**2-self.global_od_array**2)                 # Annular flow area between the wellbore and the pipe

        # Prelim calculations:
        self.J = (0.5/4) * (self.global_od_array**2 + self.global_id_array** 2) * self.global_mass_array        # Moment of Inertia
        self.J_polar = np.pi / 32 * (self.global_od_array**4 - self.global_id_array** 4)                        # Polar moment of Inertia
        self.ka = self.E * self.A_cross / self.global_length_array                                              # Axial stiffness 
        self.kt = self.G * self.J_polar / self.global_length_array                                              # Torsional stiffness

        # Parallel stifness technique
        hwdp_mask = self.global_types == "HWDP"            # Mask for HWDP components
        collar_mask = self.global_types == "Collar"        # Mask for Collar components

        if np.any(hwdp_mask):
            hwdp_ka = 1 / np.sum(1 / self.ka[hwdp_mask])        # Combined axial stiffness
            hwdp_kt = 1 / np.sum(1 / self.kt[hwdp_mask])        # Combined torsional stiffness
            self.ka[hwdp_mask] = hwdp_ka * len(self.BHA_TYPES[self.BHA_TYPES == "HWDP"]) 
            self.kt[hwdp_mask] = hwdp_kt * len(self.BHA_TYPES[self.BHA_TYPES == "HWDP"]) 

        if np.any(collar_mask):
            collar_ka = 1 / np.sum(1 / self.ka[collar_mask])    # Combined axial stiffness
            collar_kt = 1 / np.sum(1 / self.kt[collar_mask])    # Combined torsional stiffness
            self.ka[collar_mask] = collar_ka * len(self.BHA_TYPES[self.BHA_TYPES == "Collar"]) 
            self.kt[collar_mask] = collar_kt * len(self.BHA_TYPES[self.BHA_TYPES == "Collar"]) 

        # Build and Turn rates calculation
        self.bw_pipe = self.bf*self.global_mass_array/self.global_length_array
        self.MD = np.insert(np.cumsum(self.global_length_array), 0, 0)
        self.inc, self.azi, self.Normal_force = survey_mod(df_SRV, m2ft(self.MD), self.bf, self.global_mass_array, self.g)
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

        self.DIA_EQ = (27*self.global_od_array + 3*in2m(self.TJ_OD)) / 30     # in to m
        AXIAL_VEL_MULTIPLIER = self.global_od_array**2 / (self.HOLE_OD**2 - self.global_od_array**2)    # Accounting for mud velocity drag effects along axial direction
        DOC_SS = ROP_SS / RPM_SS   # m/rev
        # units of CCS of formation in ksi
        # units of k_CCS are in '1/ksi'
        MU_ROCK = -0.349 * np.log(self.ccs) + 2.0436
        # mu_rock = -0.0201 * CCS + 1.5333
        # coefficient of friction for different rock-strength
        self.K_WOB = 0.8 * (self.ccs*0.5) * lbf2N(WOB_SS) / DOC_SS * (m2in(self.HOLE_OD) / 12.25)    # units of k_WOB are in (N-rev)/(m)
        self.K_TQ = MU_ROCK / 3 * (self.K_WOB / 0.8) * self.HOLE_OD                            # units of k_TQ are in (N-rev)

        self.CA_BOREHOLE = self.CT_BOREHOLE * (
        (ROP_SS * AXIAL_VEL_MULTIPLIER) / (RPM_SS) / (self.DIA_EQ * np.pi)
        )
        self.global_ct_array = np.where(self.global_length_array == 0, 0, self.CT_BOREHOLE / self.global_length_array)
        self.global_ca_array = np.where(self.global_length_array == 0, 0, self.CA_BOREHOLE / self.global_length_array)
        
        # Sparse Matrices       
        self.global_ka_matrix = sps.diags(
            [-self.ka[1:], self.ka + np.append(self.ka[1:], 0), -self.ka[1:]],
            offsets=[-1, 0, 1],
            format='csr'
        )
        self.global_kt_matrix = sps.diags(
            [-self.kt[1:], self.kt + np.append(self.kt[1:], 0), -self.kt[1:]],
            offsets=[-1, 0, 1],
            format='csr'
        )
        self.global_mass_inv_matrix = sps.diags(1 / self.global_mass_array, format='csr')
        self.global_inertia_inv_matrix = sps.diags(1 / self.J, format='csr')
        self.global_mass_inv_ka_matrix = self.global_mass_inv_matrix @ self.global_ka_matrix        # Sparse Matrix
        self.global_inertia_inv_kt_matrix = self.global_inertia_inv_matrix @ self.global_kt_matrix  # Sparse Matrix
        