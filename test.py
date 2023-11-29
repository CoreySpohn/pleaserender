import copy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pleaserender.core import Frame, Scatter

# Create some simple data for the plots
x_values = np.linspace(0, 10, 100)
data1 = pd.DataFrame({"x": x_values, "y": np.sin(x_values)})
data2 = pd.DataFrame({"x": x_values, "y": np.cos(x_values)})

# Create Plot instances
plot1 = Scatter(data1, animation_style="Cumulative")
plot2 = Scatter(data2, animation_style="Trailing")
plot3 = Scatter(data2, animation_style="Single point")
plot4 = Scatter(data2, animation_style="Cumulative")

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
main_frame.please_set_animation_values(x_values, "x")
main_frame.please_render()
