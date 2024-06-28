MD_US = [6750, 7050]
inc = [30.9, 43.2]
azi = [78.2, 100.6]
bw_pipe = [16]
friction_factor = 0.3
import numpy as np
import pandas as pd

inc_rad = np.deg2rad(inc)
azi_rad = np.deg2rad(azi)

# Minimum Curvature Method
betta = [np.arccos(np.sin(inc_rad[i])*np.sin(inc_rad[i+1])*np.cos(azi_rad[i+1]-azi_rad[i]) + np.cos(inc_rad[i])*np.cos(inc_rad[i+1])) for i in range(len(MD_US)-1)]
RF = [(MD_US[i+1]-MD_US[i])/betta[i]*np.tan(betta[i]/2) for i in range(len(betta))]
RF = np.nan_to_num(RF)
DLS = [np.degrees(betta)[i]/(MD_US[i+1] - MD_US[i]) for i in range(len(betta))]

del_x = [RF[i]*(np.sin(inc_rad[i])*np.cos(azi_rad[i]) + np.sin(inc_rad[i+1])*np.cos(azi_rad[i+1])) for i in range(len(RF))]
del_y = [RF[i]*(np.sin(inc_rad[i])*np.sin(azi_rad[i]) + np.sin(inc_rad[i+1])*np.sin(azi_rad[i+1])) for i in range(len(RF))]
del_z = [RF[i]*(np.cos(inc_rad[i])+np.cos(inc_rad[i+1])) for i in range(len(RF))]
MD_s = [(MD_US[i]+MD_US[i+1])/2 for i in range(len(MD_US)-1)]
t_z_min = [np.cos(inc_rad[i])*np.cos(np.deg2rad(DLS[i]*(MD_s[i]-MD_US[i])))+((np.cos(inc_rad[i+1])-np.cos(inc_rad[i])*np.cos(betta[i]))/np.sin(betta[i]))*np.sin(np.deg2rad(DLS[i]*(MD_s[i]-MD_US[i]))) for i in range(len(MD_US)-1)]
inc_ave = np.rad2deg(np.arccos(t_z_min))
t_z_min = [np.cos(inc_rad[i])*np.cos(np.deg2rad(DLS[i]*(MD_s[i]-MD_US[i])))+((np.cos(inc_rad[i+1])-np.cos(inc_rad[i])*np.cos(betta[i]))/np.sin(betta[i]))*np.sin(np.deg2rad(DLS[i]*(MD_s[i]-MD_US[i]))) for i in range(len(MD_US)-1)]
n_z_min = [-np.cos(inc_rad[i])*np.sin(np.deg2rad(DLS[i]*(MD_s[i]-MD_US[i])))+((np.cos(inc_rad[i+1])-np.cos(inc_rad[i])*np.cos(betta[i]))/np.sin(betta[i]))*np.cos(np.deg2rad(DLS[i]*(MD_s[i]-MD_US[i]))) for i in range(len(MD_US)-1)]
b_z_min = [np.sin(inc_rad[i])*np.sin(inc_rad[i+1])*np.sin(azi_rad[i+1]-azi_rad[i])/np.sin(betta[i]) for i in range(len(MD_US)-1)]

def approx(F_b):
    final_appx = [F_b]
    for i in range(-1, -len(MD_US), -1):
        F_t = final_appx[-1] + bw_pipe[i]*del_z[i]
        tol = 1
        F_prev = 0
        while tol>=0.01:
            w_c = ((F_t*np.radians(DLS[i]) + bw_pipe[i]*n_z_min[i])**2 + (bw_pipe[i]*b_z_min[i])**2)**0.5
            F_t_apprx = final_appx[-1] + bw_pipe[i]*del_z[i] + friction_factor*w_c*(MD_US[i]-MD_US[i-1])
            F_t = (F_t_apprx + final_appx[-1])/2
            tol = abs((F_prev - F_t_apprx)/F_t_apprx)
            F_prev = F_t_apprx
        final_appx.append(F_t_apprx)
    return final_appx

teze = approx(25000)
teze