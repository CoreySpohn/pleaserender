import copy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
import xarray as xr
from exoverses.exovista import ExovistaSystem

from pleaserender.core import Figure, Plot, Scatter

# Create some simple data for the plots
t_values = np.linspace(0, 10, 100)
ds = xr.Dataset(coords={"time": t_values})
ds["x"] = ("time", np.sin(ds["time"].values))
ds["y"] = ("time", np.cos(ds["time"].values))

# Create Plot instances
plot1 = Plot(
    axis_keys={"x": "time", "y": "x"},
    animation_kwargs={"animation_style": "Cumulative"},
)
plot2 = Plot(
    axis_keys={"x": "x", "y": "y", "z": "time"},
    # animation_style="Trailing",
    ax_kwargs={
        "ylabel": "What the hell",
        "ylim": (-1, 1),
        "title": "'Test {animation_value:.2f}'",
    },
    animation_kwargs={
        "animation_style": "Trailing",
        "rotate": {"azim": (15, 75)},
    },
)
# plot3 = Scatter(
#     axis_keys={"x": "x", "y": "y"},
#     animation_style="Single point",
# )
# plot4 = Scatter(axis_keys={"x": "y", "y": "y"}, animation_style="Cumulative")


# Create a figure and add the plots
figure_kwargs = {"figsize": (20, 10), "layout": None}
gs_kwargs = {"height_ratios": [1, 2], "hspace": 1}
main_figure = Figure(
    nrows=2,
    ncols=1,
    fig_kwargs=figure_kwargs,
    gs_kwargs=gs_kwargs,
)
# new_figure = Figure(nrows=4, ncols=1)
# extra_figure = Figure(ncols=2)


# Add plots to the figures
main_figure.please_add_plot(plot1, row=0, col=0)
main_figure.please_add_plot(plot2, row=1, col=0)
# new_figure.please_add_plot(plot3, row=1, col=0)
# new_figure.please_add_plot(plot4, row=2, col=0, rowspan=2)
# extra_figure.please_add_plot(copy.deepcopy(plot1), row=0, col=0)
# extra_figure.please_add_plot(plot2, row=0, col=1)

# Add subfigures to the figures
# new_figure.please_add_subfigure(extra_figure, row=0, col=0)
# main_figure.please_add_subfigure(new_figure, row=0, col=1)

# Set the animation values and then render
font = {"size": 16}
plt.rc("font", **font)
plt.style.use("dark_background")
main_figure.please_add_dataset(ds)
main_figure.please_set_animation_values(t_values, "time")
# main_figure.please_preview(0)
# main_figure.please_render_images(Path("renders/main_figure/"))
render_settings = {"animation_duration": 3}
main_figure.please_render_video(
    Path("renders/test.mp4"), render_settings=render_settings
)
