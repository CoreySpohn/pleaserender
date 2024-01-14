import matplotlib.pyplot as plt
import numpy as np

from pleaserender.core.util import filter_data


class Plot:
    def __init__(
        self,
        axis_keys=None,
        # animation_style="Cumulative",
        plot_kwargs=None,
        ax_kwargs=None,
        animation_kwargs=None,
    ):
        self.ax = None
        # self.animation_style = animation_style
        self.plot_method = "plot"

        default_axis_keys = {"x": "x", "y": "y", "z": None}
        self.axis_keys = default_axis_keys
        if axis_keys is not None:
            self.axis_keys.update(axis_keys)
        self.given_axis_keys = [
            key for key, val in self.axis_keys.items() if val is not None
        ]
        self.is_3d = len(self.given_axis_keys) == 3
        self.projection = "3d" if self.is_3d else None

        # Set the default plot_kwargs, then update with any user-specified
        # kwargs. The plot_kwargs are used in the plot call.
        default_plot_kwargs = {
            "color": "k",
        }
        self.plot_kwargs = default_plot_kwargs
        if plot_kwargs is not None:
            self.plot_kwargs.update(plot_kwargs)

        # Set the default ax_kwargs, then update with any user-specified
        # kwargs. The ax_kwargs are used to set the axes properties with the
        # ax.set method. Instead of using all the ax.set_xlabel, ax.set_title,
        # ax.set_yticks, etc. functions we can pass in a dictionary that has
        # {'xlabel': 'my_xlabel', 'xlim': (0, 10), 'title': "My title"}.
        default_ax_kwargs = {
            "xlabel": self.axis_keys.get("x"),
            "ylabel": self.axis_keys.get("y"),
        }
        self.ax_kwargs = default_ax_kwargs
        if ax_kwargs is not None:
            self.ax_kwargs.update(ax_kwargs)

        # animation_kwargs are properties that specify something about the animation.
        # To rotate in azim, elevation, or roll provide a dictionary such as
        # {'rotate': {'azim': (15, 75)}} to rotate from 15 to 75 degrees in azim.
        default_animation_kwargs = {
            "animation_style": "Cumulative",
            "rotate": None,
            "elev": 30,
            "azim": 45,
            "roll": 0,
        }
        self.animation_kwargs = default_animation_kwargs
        if animation_kwargs is not None:
            self.animation_kwargs.update(animation_kwargs)

    def draw_plot(self, dataset, animation_value, animation_key):
        plot_method = getattr(self.ax, self.plot_method)
        # if self.axis_keys["z"] is None:
        data = filter_data(
            dataset,
            animation_value,
            animation_key,
            self.animation_kwargs["animation_style"],
        )
        # Allows for 2 and 3 dimensional data with the same call
        separated_data = [
            data[self.axis_keys[axis_key]] for axis_key in self.given_axis_keys
        ]
        plot_method(*separated_data, **self.plot_kwargs)

        # Evaulate any f-strings provided in ax_kwargs
        for key, val in self.ax_kwargs.items():
            if "animation_value" in val:
                self.ax_kwargs[key] = eval(f"f{val}")

        # Set all the axes properties based on ax_kwargs
        self.ax.set(**self.ax_kwargs)

        if self.is_3d:
            frame_view = {
                "elev": self.animation_kwargs["elev"],
                "azim": self.animation_kwargs["azim"],
                "roll": self.animation_kwargs["roll"],
            }
            if self.animation_kwargs["rotate"] is not None:
                for key, val in self.animation_kwargs["rotate"].items():
                    current_index = np.argmax(
                        dataset[animation_key].values == animation_value
                    )
                    frame_view[key] = np.linspace(*val, dataset[animation_key].size)[
                        current_index
                    ]
            self.ax.view_init(**frame_view)

    def verify_data(self, dataset, animation_values, animation_key):
        """
        Method to check if all the data required for the plot exists.
        Args:
            dataset (xr.Dataset):
                Dataset containing the data and coordinates
            animation_values (numpy.ndarray):
                List of values for the animation
            animation_key (str):
                Key for the animation values in the plot object's dataframe
        """
        animation_data = dataset[animation_key].values

        # Check if each element in animation_values is in animation_data
        all_values_present = np.isin(animation_values, animation_data).all()
        return all_values_present

    def generate_data(self, animation_values, animation_key):
        """
        Method to generate all the data required. Not implemented in the
        core classes.
        """
        raise NotImplementedError(
            "generate_data method must be implemented by plot/ classes."
        )

    def create_axes_config(self, data):
        """
        Method that fills ax_kwargs with defaults for things not already
        specified.
        Currently sets the axis limits.
        """
        # Set the axis limits if they are not specified
        necessary_axes = ["x", "y"]
        if self.axis_keys.get("z") is not None:
            necessary_axes.append("z")

        for axis in necessary_axes:
            if self.ax_kwargs.get(f"{axis}lim") is None:
                # Set the axis limits to be offset from the min and max values
                # of the data if they are not specified
                min_value = data[self.axis_keys.get(axis)].min()
                max_value = data[self.axis_keys.get(axis)].max()
                offset = 0.025 * np.abs(max_value - min_value)
                self.ax_kwargs[f"{axis}lim"] = (min_value - offset, max_value + offset)
