import copy

import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
import xarray as xr

from pleaserender.core import Frame, Plot, Scatter

# Create some simple data for the plots
t_values = np.linspace(0, 10, 100)
# data1 = pd.DataFrame({"x": x_values, "y": np.sin(x_values)})
# data2 = pd.DataFrame({"x": x_values, "y": np.cos(x_values)})
# data1 = xr.DataArray(np.sin(x_values), coords={"x": x_values})
# data2 = xr.DataArray(np.cos(x_values), coords={"x": x_values})
ds = xr.Dataset(coords={"time": t_values})
ds["x"] = ("time", np.sin(ds["time"].values))
ds["y"] = ("time", np.cos(ds["time"].values))

# Create Plot instances
# plot1 = Scatter(axis_keys={"x": "time", "y": "x"}, animation_style="Cumulative")
# plot2 = Scatter(axis_keys={"x": "time", "y": "y"}, animation_style="Trailing")
# plot3 = Scatter(axis_keys={"x": "x", "y": "y"}, animation_style="Single point")
# plot4 = Scatter(axis_keys={"x": "y", "y": "y"}, animation_style="Cumulative")
plot1 = Plot(
    axis_keys={"x": "time", "y": "x"},
    animation_style="Cumulative",
)
plot2 = Plot(axis_keys={"x": "time", "y": "y"}, animation_style="Trailing")
plot3 = Scatter(
    axis_keys={"x": "x", "y": "y"},
    animation_style="Single point",
)
plot4 = Scatter(axis_keys={"x": "y", "y": "y"}, animation_style="Cumulative")

# Create a Frame and add the plots
main_frame = Frame(nrows=1, ncols=2)
new_frame = Frame(nrows=4, ncols=1)
extra_frame = Frame(ncols=2)

# Add plots to the frames
main_frame.please_add_plot(plot1, row=0, col=0)
new_frame.please_add_plot(plot3, row=1, col=0)
new_frame.please_add_plot(plot4, row=2, col=0, rowspan=2)
extra_frame.please_add_plot(copy.deepcopy(plot1), row=0, col=0)
extra_frame.please_add_plot(plot2, row=0, col=1)

# Add subframes to the frames
new_frame.please_add_subframe(extra_frame, row=0, col=0)
main_frame.please_add_subframe(new_frame, row=0, col=1)

# Set the animation values and then render
main_frame.please_add_dataset(ds)
main_frame.please_set_animation_values(t_values, "time")
main_frame.please_render()
