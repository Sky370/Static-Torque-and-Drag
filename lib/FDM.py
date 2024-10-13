from lib.constants import Calculations
from lib.topdrive import vel_acc
from lib.circulation import circl
from scipy.integrate import solve_ivp
import numpy as np

clc = Calculations(depth=length)
initial_vel, initial_acc = vel_acc(clc.time_array, clc)

# Initialize arrays for simulation variables
displacement        = np.zeros((len(clc.time_array), len(clc.global_mass_array)))
velocity            = np.zeros((len(clc.time_array), len(clc.global_mass_array)))
axial_force         = np.zeros((len(clc.time_array), len(clc.global_mass_array)))
inner_pressure      = np.zeros((len(clc.time_array), len(clc.global_mass_array)))
annular_pressure    = np.zeros((len(clc.time_array), len(clc.global_mass_array)))
friction_factor     = np.zeros((len(clc.time_array), len(clc.global_mass_array)))

# Time steps, nums, etc.
dt = 0.005
time = np.arange(0, clc.time_val, dt)
num_steps = len(time)

# Initialize variables based on the given equation
rho_s_As = clc.global_mass_array / clc.global_length_array  # Mass per unit length for all segments
E_As = clc.E*clc.A_cross                                    # Elastic modulus * area of cross-section for all segments
delta_s = clc.global_length_array                           # Segment lengths for all elements

for step in range(1, num_steps):
    for i in range(1, clc.noe - 1):  # Skipping boundaries      

        # Calculate the right-hand side terms for axial forces
        F_axial = (E_As[i] * (displacement[step-1, i+1] + displacement[step-1, i]) / delta_s[i] -
                   E_As[i-1] * (displacement[step-1, i] - displacement[step-1, i-1]) / delta_s[i-1])

        # External forces (friction, pressure, etc.)
        F_external = external_force[i] + friction_factor[step, i] * normal_force[i] + (
            inner_pressure[step-1, i] * clc.A_i[i] - outer_pressure[step-1, i] * clc.A_o[i])

        # Calculate U_{i,j+1} based on the derived equation
        U_j_next = (2 * displacement[step-1, i] - displacement[step-2, i] +
                    (dt**2 / (rho_s_As[i] * delta_s[i] / 2 + rho_s_As[i-1] * delta_s[i-1] / 2)) *
                    (F_axial + F_external))

        # Update displacement and velocity
        displacement[step, i] = U_j_next
        velocity[step, i] = (displacement[step, i] - displacement[step-1, i]) / dt









for step, time in enumerate(clc.time_array, start=1):
    for elm in 
    velocity[step,:] = initial_vel[step]
    pipe_acceleration = initial_acc[step]
    # Update friction factors (vectorized)
    friction_factor[step, :] = np.where(np.abs(velocity) < delta,
                                        clc.mu_s*np.sign(pipe_velocity),
                                        clc.mu_d*np.sign(pipe_velocity)
    )

    # Pressure drop calculations (vectorized)
    dPdL_an, dPdL_in, v_in, v_an = pres_calc(
        rho=calculations.mud_density_ppg,
        mu_p=calculations.visc_p,
        tao=calculations.tao_y,
        Q=some_flow_rate,         # Provide appropriate flow rate
        D_o=global_od_array,
        D_i=global_id_array,
        D_w=global_hole_array
    )

    # Update pressures (vectorized)
    inner_pressure[step, :] = inner_pressure[step-1, :] - dPdL_in * segment_length
    annular_pressure[step, :] = annular_pressure[step-1, :] - dPdL_an * segment_length

    # Compute contact forces (vectorized)
    # For example, using your bw_pipe and DLS arrays
    contact_force = calculations.bw_pipe * calculations.DLS  # Adjust as needed

    # Update axial forces (vectorized)
    axial_force[step, :] = axial_force[step-1, :] + dt * (
        friction_factor[step, :] * contact_force / global_mass_array
    )

    # Update velocities and displacements (vectorized)
    velocity[step, :] = pipe_velocity  # Assuming rigid body motion
    displacement[step, :] = displacement[step-1, :] + dt * velocity[step, :]

    # Apply boundary conditions
    apply_boundary_conditions(
        displacement=displacement,
        velocity=velocity,
        axial_force=axial_force,
        U_hook=trapezoidal_hook_displacement,
        WOB=WOB,
        step=step,
        n_segments=len(global_mass_array)
    )
    fric_factor = something



def trip_operation(F_b, flowrate, length, operation):
    clc = Calculations(depth=length)
    fluid_force = circl(flowrate, length) if flowrate != 0 else np.zeros_like(clc.MD[:-1])

    def dFds(s, F, velocity, acceleration):
        cross_area = np.pi*(clc.global_od_array**2 - clc.global_id_array**2)/4
        
        # Interpolate the necessary values for the current s
        DLS_interp = np.interp(s, clc.MD[:-1], clc.DLS)
        n_z_interp = np.interp(s, clc.MD[:-1], clc.n_z)
        b_z_interp = np.interp(s, clc.MD[:-1], clc.b_z)
        t_z_interp = np.interp(s, clc.MD[:-1], clc.t_z)
        bw_pipe_interp = np.interp(s, clc.MD[:-1], clc.bw_pipe)
        fluid_interp = np.interp(s, clc.MD[:-1], fluid_force)

        # Calculate the contact force with new dynamic terms
        vel_term = 12/231*clc.pipe_density * get_area(s, clc.MD, cross_area) * velocity**2                       # Velocity term
        acc_term = 12/231*clc.pipe_density * get_area(s, clc.MD, cross_area) * acceleration                      # Acceleration term
        w_c = np.sqrt((F * DLS_interp + bw_pipe_interp * n_z_interp - vel_term)**2 + (bw_pipe_interp * b_z_interp)**2)

        if operation == 'POOH':
            return -(bw_pipe_interp*t_z_interp + clc.ff*w_c + fluid_interp + acc_term)
        elif operation == 'RIH':
            return -(bw_pipe_interp*t_z_interp - clc.ff*w_c + fluid_interp + acc_term)
        elif operation == 'ROB':
            return -(bw_pipe_interp*t_z_interp + fluid_interp)
        else:
            raise ValueError("Unknown operation type specified.")

    s_span = [clc.MD[-1], clc.MD[0]]  # Span from the end to the start
        
    if operation == 'POOH':
        times = np.arange(0, clc.a3, 0.5)           # Time domain
        vel, acc = vel_acc(times, clc)                   # Computing velocity and acceleration
        
    elif operation == 'RIH':
        times = np.arange(clc.a3, clc.a6, 0.5)      # Time domain
        vel, acc = vel_acc(times, clc)                   # Computing velocity and acceleration

    elif operation == 'ROB':
        times = np.arange(0, 30, 0.5)               # Time domain
        vel, acc = (
            np.zeros(len(times)), 
            np.zeros(len(times))
            )                                       # Computing velocity and acceleration
    else:
        raise ValueError("Unknown operation type specified.")
    
    sol = [solve_ivp(dFds, 
                    s_span, 
                    [F_b], 
                    args=(vel[i], acc[i]), 
                    method='RK45', 
                    t_eval=np.flip(clc.MD),
                    max_step=10,
                    rtol=0.001,
                    atol=0.001,
                    dense_output=True
                    ).y.flatten()[-1] for i,_ in enumerate(times)]
    
    return sol, times