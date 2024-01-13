import copy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
import xarray as xr

from pleaserender.core import Figure, Plot, Scatter

# Create some simple data for the plots
t_values = np.linspace(0, 10, 100)
ds = xr.Dataset(coords={"time": t_values})
ds["x"] = ("time", np.sin(ds["time"].values))
ds["y"] = ("time", np.cos(ds["time"].values))

# Create Plot instances
plot1 = Plot(
    axis_keys={"x": "time", "y": "x"},
    animation_style="Cumulative",
)
plot2 = Plot(
    axis_keys={"x": "time", "y": "y"},
    animation_style="Trailing",
    ax_kwargs={"ylabel": "What the hell", "ylim": (-1, 1), "title": "Test"},
)
# plot3 = Scatter(
#     axis_keys={"x": "x", "y": "y"},
#     animation_style="Single point",
# )
# plot4 = Scatter(axis_keys={"x": "y", "y": "y"}, animation_style="Cumulative")


# Create a figure and add the plots
main_figure = Figure(nrows=2, ncols=1)
# new_figure = Figure(nrows=4, ncols=1)
# extra_figure = Figure(ncols=2)


# Add plots to the figures
main_figure.please_add_plot(plot1, row=0, col=0)
main_figure.please_add_plot(plot2, row=1, col=0, sharex_plot=plot1)
# new_figure.please_add_plot(plot3, row=1, col=0)
# new_figure.please_add_plot(plot4, row=2, col=0, rowspan=2)
# extra_figure.please_add_plot(copy.deepcopy(plot1), row=0, col=0)
# extra_figure.please_add_plot(plot2, row=0, col=1)

# Add subfigures to the figures
# new_figure.please_add_subfigure(extra_figure, row=0, col=0)
# main_figure.please_add_subfigure(new_figure, row=0, col=1)

# Set the animation values and then render
main_figure.please_add_dataset(ds)
main_figure.please_set_animation_values(t_values, "time")
main_figure.please_render(Path("renders/main_figure/"))
# new_figure.please_set_animation_values(t_values, "time")
# new_figure.please_render(Path("renders/new_figure/"))
