import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from lib.tripping import trip_operation
from lib.init_xls import N2lbf, Nm2lbfft

def plot_profiles_vs_gpm(
    total_length: float,
    wob: float = 0.0,
    tob: float = 0.0,
    gpm_list: list[float] = None,
    operation: str = "ROB",
    *,
    html_out: str = None,
    image_out: str = None,
    image_scale: int = 1,
    show: bool = True
) -> go.Figure:

    if not gpm_list:
        raise ValueError("`gpm_list` must be a non-empty list of flow rates.")

    profiles = []
    for gpm in gpm_list:
        F, T = trip_operation(
            F_b=wob,
            T_b=tob,
            GPM=gpm,
            length=total_length,
            operation=operation
        )
        profiles.append((gpm, F, T))

    # shared MD axis (use the first profile length)
    md_axis = np.linspace(total_length, 0, len(profiles[0][1]))

    fig = make_subplots(
        rows=1, cols=1,
        subplot_titles=(operation,),
        vertical_spacing=0,
        horizontal_spacing=0
    )

    for gpm, F, T in profiles:
        fig.add_trace(
            go.Scatter(
                x=N2lbf(F),
                y=md_axis,
                name=f"{int(gpm)} GPM",
                line=dict(width=3, dash="dash")
            ),
            row=1, col=1
        )

    fig.update_xaxes(
        title_text="Axial Force [lbf]",
        title_font=dict(size=20, family="Arial"),
        showgrid=True, gridcolor="lightgray",
        tickfont=dict(size=14)
    )
    fig.update_yaxes(
        title_text="MD [ft]",
        title_font=dict(size=20, family="Arial"),
        autorange="reversed",
        showgrid=True, gridcolor="lightgray",
        tickfont=dict(size=14)
    )

    fig.update_layout(
        template="plotly_white",
        title_text=f"Hookload vs MD for {operation} at {', '.join(str(int(g)) for g in gpm_list)} GPM",
        title_font=dict(size=24, family="Arial"),
        # legend=dict(orientation="h", yanchor="bottom", y=1.02, font=dict(size=14)),
        legend=dict(traceorder='normal', font=dict(size=16)),
        margin=dict(l=80, r=40, t=80, b=60),
        width=800,
        height=600
    )
    if html_out:
        fig.write_html(html_out)
    if image_out:
        fig.write_image(image_out, scale=image_scale)
    if show:
        fig.show()

    return fig
