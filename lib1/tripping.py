from .constants import Calculations
import numpy as np

def trip_POOH(F_b, ff, length):
    clc = Calculations(lengths=length)
    final_appx = [F_b]
    for i in range(-1, -len(clc.MD), -1):
        F_t = final_appx[-1] + clc.bw_pipe[i]*clc.del_z[i]
        tol = 1
        F_prev = 0
        while tol>=0.01:
            w_c = ((F_t*np.radians(clc.DLS[i]) + clc.bw_pipe[i]*clc.n_z_min[i])**2 + (clc.bw_pipe[i]*clc.b_z_min[i])**2)**0.5
            F_t_apprx = final_appx[-1] + clc.bw_pipe[i]*clc.del_z[i] + ff*w_c*(clc.MD[i]-clc.MD[i-1])
            F_t = (F_t_apprx + final_appx[-1])/2
            tol = abs((F_prev - F_t_apprx)/F_t_apprx)
            F_prev = F_t_apprx
        final_appx.append(F_t_apprx)
    diff = np.diff(final_appx)
    return diff

def trip_RIH(F_b, ff, length):
    clc = Calculations(lengths=length)
    final_appx = [F_b]
    for i in range(-1, -len(clc.MD), -1):
        F_t = final_appx[-1] + clc.bw_pipe[i]*clc.del_z[i]
        tol = 1
        F_prev = 0
        while tol>=0.01:
            w_c = ((F_t*np.radians(clc.DLS[i]) + clc.bw_pipe[i]*clc.n_z_min[i])**2 + (clc.bw_pipe[i]*clc.b_z_min[i])**2)**0.5
            F_t_apprx = final_appx[-1] + clc.bw_pipe[i]*clc.del_z[i] - ff*w_c*(clc.MD[i]-clc.MD[i-1])
            F_t = (F_t_apprx + final_appx[-1])/2
            tol = abs((F_prev - F_t_apprx)/F_t_apprx)
            F_prev = F_t_apprx
        final_appx.append(F_t_apprx)
    diff = np.diff(final_appx)
    return diff
