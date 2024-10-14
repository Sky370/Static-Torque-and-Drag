from .init_xls import *
import scipy.sparse as sps
import warnings
warnings.filterwarnings('ignore')

class Calculations:
    def __init__(self, depth):
        # Field units
        self.gravity = 32.1740  # ft/s²
        self.bit_depth = depth # ft
        self.elem_length = df_BHA[df_BHA["BHA Type"] == "DP"]["Length (ft)"].iloc[0] # ft
        self.mud_density_ppg = float(df_WELL[df_WELL["Parameter"] == "Mud Density"]["Value"]) # ppg
        self.pipe_density = float(df_WELL[df_WELL["Parameter"] == "Steel Density"]["Value"]) # ppg
        self.visc_p = float(df_WELL[df_WELL["Parameter"] == "Plastic Viscosity"]["Value"]) # cP
        self.tao_y = float(df_WELL[df_WELL["Parameter"] == "Yield Point"]["Value"]) # lbf/100ft2
        self.bf = float((1-self.mud_density_ppg/self.pipe_density))
        self.mu_s = float(df_WELL[df_WELL["Parameter"] == "Static Friction Factor"]["Value"])     # Static Fric. 
        self.mu_d = float(df_WELL[df_WELL["Parameter"] == "Dynamic Friction Factor"]["Value"])    # Dynamic Fric.
        self.E = float(df_ADV[df_ADV["Parameter"] == "Young Modulus"]["Value"])                   # Young Modulus
        self.G = float(df_ADV[df_ADV["Parameter"] == "Bulk Modulus"]["Value"])                    # Torsional Modulus
        self.ccs = float(df_ADV[df_ADV["Parameter"] == "CCS"]["Value"])                           # Psi
        self.v_cs = float(df_WELL[df_WELL["Parameter"] == "Stribeck Critical Velocity"]["Value"]) # ft/s
        self.CT_BOREHOLE = float(df_WELL[df_WELL["Parameter"] == "Torsional Drag Coefficient"]["Value"])

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
        self.v1 = float(df_TOP[df_TOP["Parameter"] == "Top Drive Axial Velocity Magnitude 1 (ft/min)"]["Value"]) / 60
        self.v2 = float(df_TOP[df_TOP["Parameter"] == "Top Drive Axial Velocity Magnitude 2 (ft/min)"]["Value"]) / 60
        self.rpm1 = float(df_TOP[df_TOP["Parameter"] == "Top Drive RPM Magnitude 1 (RPM)"]["Value"]) / 60
        self.rpm2 = float(df_TOP[df_TOP["Parameter"] == "Top Drive RPM Magnitude 2 (RPM)"]["Value"]) / 60

        # Time array
        self.time_val = float(df_ADV[df_ADV["Parameter"] == "Run Time"]["Value"])

        # Extra
        # --------- Load steady state inputs -------- #
        WOB_SS = float(df_SS[df_SS["Parameter"] == "WOB initial"]["Value"])        # lbs
        ROP_SS = float(df_SS[df_SS["Parameter"] == "ROP steady state"]["Value"])   # ft/hr
        RPM_SS = float(df_SS[df_SS["Parameter"] == "RPM steady state"]["Value"])   

        # Drill pipes:
        self.dp_length = self.bit_depth - df_BHA["Total Length (ft)"].iloc[0]
        self.N_DP, self.L_DP = nearestLength(self.dp_length, self.elem_length)
        self.N_DP = round(self.N_DP)        
        self.DP_OD = np.ones(self.N_DP) * df_BHA["OD (in)"].iloc[0]  
        self.DP_ID = np.ones(self.N_DP) * df_BHA["ID (in)"].iloc[0]  
        self.DP_MASS_ARRAY = np.ones(self.N_DP) * df_BHA["Mass (lbs)"].iloc[0]  
        self.L_DP_ARRAY = np.ones(self.N_DP) * self.L_DP
        self.TJ_OD = df_BHA["OD Tool Joint (in)"].iloc[0]   # Tool Joint OD, in

        # Collars / BHA:
        self.COLLAR_LEN = np.array(np.repeat(df_BHA["Length (ft)"]/df_BHA["Number of Items"] , df_BHA["Number of Items"]))
        self.COLLAR_MASS = np.array(np.repeat(df_BHA["Mass (lbs)"]/df_BHA["Number of Items"] , df_BHA["Number of Items"]))
        self.COLLAR_OD = np.array(np.repeat(df_BHA["OD (in)"], df_BHA["Number of Items"]))
        self.COLLAR_ID = np.array(np.repeat(df_BHA["ID (in)"], df_BHA["Number of Items"])) 

        # Bottom hole
        self.HOLE_OD = df_ADV[df_ADV["Parameter"] == "Hole Diameter"]["Value"].iloc[0]
        self.HOLE_DEPTH = df_ADV[df_ADV["Parameter"] == "Hole Depth"]["Value"].iloc[0]
        self.noe = self.N_DP + len(self.COLLAR_OD)
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
        self.global_mass_array = np.concatenate([self.DP_MASS_ARRAY, self.COLLAR_MASS])
        self.global_length_array = np.concatenate([self.L_DP_ARRAY, self.COLLAR_LEN])
        self.global_od_array = np.concatenate([self.DP_OD, self.COLLAR_OD])
        self.global_id_array = np.concatenate([self.DP_ID, self.COLLAR_ID])
        self.global_eps = (self.HOLE_ARRAY/self.global_od_array) - 1
        
        # Area calculation
        self.A_i = np.pi/4*(self.global_id_array)**2                                    # Inner area of the pipe
        self.A_o = np.pi/4*(self.global_od_array)**2                                    # Outer area of the pipe
        self.A_cross = np.pi/4*(self.global_od_array**2-self.global_id_array**2)        # Cross-sectional area of the pipe
        self.A_h = np.pi/4*(self.HOLE_ARRAY**2-self.global_od_array**2)          # Annular flow area between the wellbore and the pipe

        # Build and Turn rates calculation
        self.bw_pipe = self.bf*self.global_mass_array/self.global_length_array
        self.MD = np.insert(np.cumsum(self.global_length_array), 0, 0)
        self.inc, self.azi, self.Normal_force = survey_mod(df_SRV, self.MD, self.bf, self.global_mass_array, self.gravity)
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

         # Prelim calculations:
        self.J = (0.5/4) * (self.global_od_array**2 + self.global_id_array** 2) * self.global_mass_array                        # Moment of Inertia
        self.ka = self.E * self.A_cross / self.global_length_array                                                              # Axial stiffness [imperial]
        self.kt = self.G * (np.pi / 32 * (self.global_od_array**4 - self.global_id_array** 4)) / self.global_length_array       # Torsional stiff. [imperial]

        self.DIA_EQ = (27*self.global_od_array + 3*self.TJ_OD) / 30     # in
        AXIAL_VEL_MULTIPLIER = self.global_od_array**2 / (self.HOLE_ARRAY**2 - self.global_od_array**2)    # Accounting for mud velocity drag effects along axial direction
        DOC_SS = ROP_SS / RPM_SS * 12 / 60  # in/rev
        # units of CCS of formation in ksi
        # units of k_CCS are in '1/ksi'
        MU_ROCK = -0.349 * np.log(self.ccs) + 2.0436
        # mu_rock = -0.0201 * CCS + 1.5333
        # coefficient of friction for different rock-strength
        self.K_WOB = 0.8 * (self.ccs*0.5) * WOB_SS / DOC_SS * (self.HOLE_ARRAY[0] / 12.25)      # units of k_WOB are in (lbf-rev)/(in)
        self.K_TQ = MU_ROCK / 3 * (self.K_WOB / 0.8) * (self.HOLE_ARRAY[0])                     # units of k_TQ are in (lbf-rev)

        self.CA_BOREHOLE = self.CT_BOREHOLE * (12/60) * (ROP_SS * AXIAL_VEL_MULTIPLIER / (RPM_SS * self.DIA_EQ * np.pi))
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
        # self.global_mass_inv_ca_matrix = sps.diags(1 / self.global_mass_array * self.global_ca_array, format='csr')
        # self.global_inertia_inv_ct_matrix = sps.diags(1 / self.J * self.global_ct_array, format='csr')
        