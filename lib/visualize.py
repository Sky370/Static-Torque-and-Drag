import plotly.graph_objects as go
import numpy as np
from plotly.subplots import make_subplots
from lib.topdrive import topdrive as top_drive_input


def visualize_vel(time_arr, sol_arr, noe, p, constant, path, fric_mod):

    subplot_titles = [
        "Axial Velocities (m/hr)",
        "Weight (lbf)",
        "Depth of Cut (in)",
        "Rotational velocities (revs/min)",
        "Torque (lbf-ft)",
        "Bit depth vs hole depth (m)",
    ]
    fig = make_subplots(
        rows=2,
        cols=3,
        subplot_titles=subplot_titles,
    )
    store = {}
    store["ROP_topdrive"] = np.zeros(len(time_arr))
    store["RPM_topdrive"] = np.zeros(len(time_arr))
    store["z_topdrive"] = np.zeros(len(time_arr))
    store["theta_topdrive"] = np.zeros(len(time_arr))
    for i, t in enumerate(time_arr):
        qq1, qq2, qq3, qq4 = top_drive_input(t, constant)
        store["ROP_topdrive"][i] = qq1
        store["RPM_topdrive"][i] = qq2
        store["z_topdrive"][i] = qq3
        store["theta_topdrive"][i] = qq4
    bit_axial_vel = go.Scatter(
        x=time_arr,
        y=sol_arr[(noe - 1) * 4 + 1] * 3600,
        name="Bit Axial Velocity (m/hr)",
        mode="lines",
    )
    topdrive_axial_vel = go.Scatter(
        x=time_arr,
        y=store["ROP_topdrive"] * 3600,
        mode="lines",
        name="Topdrive Axial Velocity (m/hr)",
    )

    b_rpm = sol_arr[(noe - 1) * 4 + 3] * (30 / np.pi) 
    t_rpm = store["RPM_topdrive"] * (30 / np.pi)

    bit_rpm = go.Scatter(
        x=time_arr,
        y=b_rpm,
        name="Bit RPM",
        mode="lines",
    )
    topdrive_rpm = go.Scatter(
        x=time_arr,
        y=t_rpm,
        mode="lines",
        name="Topdrive RPM",
    )

    surf_weight = (
        constant.ka[0] * (store["z_topdrive"] - sol_arr[0 * 4 + 0]) * (1 / 4.45)
    )
    surf_torque = (
        constant.kt[0] * (store["theta_topdrive"] - sol_arr[0 * 4 + 2])
        * (1 / (4.45 * 12 * 0.0254))
    )

    down_torque = (
        -np.array(p['DOWNHOLE_TORQUE']) * (1 / (4.45 * 12 * 0.0254))
    )

    down_weight = (
        -np.array(p['DOWNHOLE_WEIGHT']) * (1 / 4.45)
    )
    downhole_weight = go.Scatter(
        x=p['SOLUTION_TIME'],
        y=down_weight,
        mode="lines",
        name="Downhole Weight (lbf)",
    )
    surface_weight = go.Scatter(
        x=time_arr, y=surf_weight, mode="lines", name="Surface Weight (lbf)"
    )
    downhole_torque = go.Scatter(
        x=p['SOLUTION_TIME'],
        y=down_torque,
        mode="lines",
        name="Downhole Torque (lbf-ft)",
    )
    surface_torque = go.Scatter(
        x=time_arr, y=surf_torque, mode="lines", name="Surface Torque (lbf-ft)"
    )
    depth_of_cut = go.Scatter(
        x=p['SOLUTION_TIME'],
        y=39.3701 * np.array(p['DOC']),
        mode="lines",
        name="Depth of Cut",
    )
    bit_axial_depth = go.Scatter(
        x=time_arr,
        y=sol_arr[-4] + constant.bit_depth,
        name="Bit Depth (m)",
        mode="lines",
    )

    hole_depth_store = go.Scatter(
        x=p['SOLUTION_TIME'],
        y=p['HOLE_DEPTH_PREV'],
        name="Hole Depth (m)",
        mode="lines",
    )

    fig.add_trace(topdrive_axial_vel, row=1, col=1)
    fig.add_trace(bit_axial_vel, row=1, col=1)
    fig.add_trace(bit_rpm, row=2, col=1)
    fig.add_trace(topdrive_rpm, row=2, col=1)
    fig.add_trace(downhole_torque, row=2, col=2)
    fig.add_trace(surface_torque, row=2, col=2)
    fig.add_trace(downhole_weight, row=1, col=2)
    fig.add_trace(surface_weight, row=1, col=2)
    fig.add_trace(depth_of_cut, row=1, col=3)
    fig.add_trace(bit_axial_depth, row=2, col=3)
    fig.add_trace(hole_depth_store, row=2, col=3)
    fig.update_xaxes(matches="x")

    keys = ["xaxis"]
    keys = keys + [f"xaxis{i}" for i in range(2, 7)]
    layout = {key: {"title": "time (sec)"} for key in keys}
    fig.update_layout(layout)
    keys = ["yaxis"]
    keys = keys + [f"yaxis{i}" for i in range(2, 7)]
    layout = {key: {"title": f"{subplot_titles[i]}"} for i, key in enumerate(keys)}

    fig.update_layout(layout)
    layout = {"yaxis6": {"autorange": "reversed"}}
    fig.update_layout(layout)
    fig.update_layout(title_text="Drilling Parameters")
    fig.write_html(f"{path}/Outputs-{fric_mod}.html")