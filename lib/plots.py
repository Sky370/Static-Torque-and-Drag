import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from lib.tripping import trip_operation
from lib.init_xls import N2lbf, Nm2lbfft

def plot_hookload_torque(
    total_length: float,
    wob: float = 0.0,
    tob: float = 0.0,
    GPM: float = 0.0,
    *,
    html_out: str = None,
    image_out: str = None,
    image_scale: int = 1,
    show: bool = True
) -> go.Figure:

    hkld_pooh, _ = trip_operation(F_b=wob, T_b=tob, GPM=GPM, length=total_length, operation="POOH")
    hkld_rih,  _ = trip_operation(F_b=wob, T_b=tob, GPM=GPM, length=total_length, operation="RIH")
    hkld_rob,  _ = trip_operation(F_b=wob, T_b=tob, GPM=GPM, length=total_length, operation="ROB")
    _, torque    = trip_operation(F_b=wob, T_b=tob, GPM=GPM, length=total_length, operation="TQ")

    md_axis = np.linspace(total_length, 0, len(hkld_pooh))
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Hookload", "Torque"),
        horizontal_spacing=0.12
    )

    # Hookload traces
    fig.add_trace(go.Scatter(
        x=N2lbf(hkld_rih), y=md_axis,
        name="RIH", line=dict(color="navy", width=3)
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=N2lbf(hkld_rob), y=md_axis,
        name="ROB", line=dict(color="seagreen", width=3)
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=N2lbf(hkld_pooh), y=md_axis,
        name="POOH", line=dict(color="violet", width=3)
    ), row=1, col=1)

    # Torque trace
    fig.add_trace(go.Scatter(
        x=Nm2lbfft(torque), y=md_axis,
        name="Torque", line=dict(color="red", width=3)
    ), row=1, col=2)

    for col in (1, 2):
        fig.update_xaxes(
            title_text=("Axial Force [lbf]" if col==1 else "Torque [lbf·ft]"),
            title_font=dict(size=22, family="Arial"),
            row=1, col=col,
            showgrid=True, gridcolor="lightgray",
            tickfont=dict(size=18),
        )
        fig.update_yaxes(
            title_text="MD [ft]",
            title_font=dict(size=22, family="Arial"),
            row=1, col=col,
            autorange="reversed",
            showgrid=True, gridcolor="lightgray",
            tickfont=dict(size=18),
        )
    fig.update_annotations(font_size=22, font_family="Arial")

    fig.update_layout(
        template="plotly_white",
        title_text="Hookload & Torque vs. Measured Depth",
        title_font=dict(size=30, family="Arial"),
        legend=dict(
            x=0.5, y=1.05, orientation="h",
            xanchor="center", font=dict(size=18)
        ),
        margin=dict(l=80, r=40, t=100, b=60),
        width=1200, height=600
    )
    
    # Export if requested
    if html_out:
        fig.write_html(html_out)
    if image_out:
        fig.write_image(image_out, scale=image_scale)

    if show:
        fig.show()

    return fig
