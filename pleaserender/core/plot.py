import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pleaserender.core.plot_object import PlotObject


class Plot:
    def __init__(
        self, plot_obj, axis_keys=None, animation_style="Cumulative", plot_kwargs=None
    ):
        self.plot_obj = PlotObject(plot_obj)
        self.ax = None
        self.animation_style = animation_style
        default_plot_kwargs = {
            "color": "k",
        }
        self.plot_kwargs = default_plot_kwargs
        if plot_kwargs is not None:
            self.plot_kwargs.update(plot_kwargs)

        default_axis_keys = {"x": "x", "y": "y", "z": None}
        self.axis_keys = default_axis_keys
        if axis_keys is not None:
            self.axis_keys.update(axis_keys)

    def set_title(self, title):
        self.ax.set_title(title)

    def set_xlabel(self, label):
        self.ax.set_xlabel(label)

    def set_ylabel(self, label):
        self.ax.set_ylabel(label)

    def draw_plot(self):
        """
        This method should be overridden by subclasses to create specific plot types.
        """
        raise NotImplementedError("draw_plot method must be implemented by subclasses")

    def show(self):
        plt.show()

    def save(self, filename, dpi=300):
        self.fig.savefig(filename, dpi=dpi)

    def close(self):
        plt.close(self.fig)

    def verify_data(self, animation_values, animation_key):
        """
        Method to check if all the data required for the plot exists.
        Args:
            animation_values (numpy.ndarray):
                List of values for the animation
            animation_key (str):
                Key for the animation values in the plot object's dataframe
        """
        all_values_present = (
            self.plot_obj.data[animation_key].isin(animation_values).all()
        )
        return all_values_present

    def generate_data(self, animation_values, animation_key):
        """
        Method to generate all the data required. Not implemented in the
        core classes.
        """
        raise NotImplementedError(
            "generate_data method must be implemented by plot/ classes."
        )
        # if self.plot_obj.is_df:
        #     # Make sure all the data is generated
        #     all_values_present = self.plot_obj.data.isin(animation_values).all()
        #     if not all_values_present:
        #         # Cannot generate data if only a dataframe is provided
        #         raise ValueError("Cannot generate data from a dataframe")
        # else:
        #     # Generate data
        #     data = getattr(
        #         self.plot_obj, self.plot_obj.data_gen_method(animation_values)
        #     )()
        #     # Insert data into the plot object's dataframe
        #     self.plot_obj.data = pd.DataFrame(data)

    def create_axes_config(self, data):
        default_axes_config = {}
        xdata = data[self.axis_keys.get("x")]
        default_axes_config["x"] = {
            "label": self.axis_keys.get("x"),
            "limits": (xdata.min(), xdata.max()),
        }
        ydata = data[self.axis_keys.get("y")]
        default_axes_config["y"] = {
            "label": self.axis_keys.get("y"),
            "limits": (ydata.min(), ydata.max()),
        }
        if self.axis_keys.get("z") is not None:
            zdata = data[self.axis_keys.get("z")]
            default_axes_config["z"] = {
                "label": self.axis_keys.get("z"),
                "limits": (zdata.min(), zdata.max()),
            }
        self.axes_config = default_axes_config

    def apply_axes_config(self):
        for axis, settings in self.axes_config.items():
            if axis == "x":
                self.ax.set_xlabel(settings["label"])
                self.ax.set_xlim(settings["limits"])
            elif axis == "y":
                self.ax.set_ylabel(settings["label"])
                self.ax.set_ylim(settings["limits"])
            elif axis == "z" and hasattr(self.ax, "set_zlabel"):
                # Only apply z-axis settings if it's a 3D plot
                self.ax.set_zlabel(settings["label"])
                self.ax.set_zlim(settings["limits"])
