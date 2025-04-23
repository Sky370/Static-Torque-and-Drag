from lib.plots import plot_hookload_torque
from lib.plots_multi import plot_profiles_vs_gpm

if __name__ == "__main__":
    fig = plot_hookload_torque(
        total_length=13000,
        wob=0.0,            # Weight on bit
        tob=0.0,            # Torque at bit
        GPM = (0),          # Flowrate
        # html_out="Plot.html",
        # image_out="Plot.png",
        # image_scale=3,
        show=True
    )

    fig = plot_profiles_vs_gpm(
        total_length=13000,
        wob=0.0,                   # Weight on bit
        tob=0.0,                   # Torque at bit
        gpm_list=[0, 250, 600],    # any number of GPMs
        operation="ROB",           # or "POOH","RIH"
        # html_out="Plot.html",
        # image_out="Plot.png",
        # image_scale=3,
        show=True
    )