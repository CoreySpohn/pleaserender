import matplotlib.pyplot as plt
import numpy as np

from pleaserender.core.util import filter_data


class Plot:
    def __init__(
        self,
        axis_keys=None,
        animation_style="Cumulative",
        plot_kwargs=None,
        ax_kwargs=None,
    ):
        self.ax = None
        self.animation_style = animation_style
        self.plot_method = "plot"

        default_axis_keys = {"x": "x", "y": "y", "z": None}
        self.axis_keys = default_axis_keys
        if axis_keys is not None:
            self.axis_keys.update(axis_keys)

        default_plot_kwargs = {
            "color": "k",
        }
        self.plot_kwargs = default_plot_kwargs
        if plot_kwargs is not None:
            self.plot_kwargs.update(plot_kwargs)

        default_ax_kwargs = {
            "xlabel": self.axis_keys.get("x"),
            "ylabel": self.axis_keys.get("y"),
        }
        self.ax_kwargs = default_ax_kwargs
        if ax_kwargs is not None:
            self.ax_kwargs.update(ax_kwargs)

    def set_title(self, title):
        self.ax.set_title(title)

    def set_xlabel(self, label):
        self.ax.set_xlabel(label)

    def set_ylabel(self, label):
        self.ax.set_ylabel(label)

    def draw_plot(self, dataset, animation_value, animation_key):
        plot_method = getattr(self.ax, self.plot_method)
        if self.axis_keys["z"] is None:
            data = filter_data(
                dataset, animation_value, animation_key, self.animation_style
            )
            plot_method(
                data[self.axis_keys["x"]], data[self.axis_keys["y"]], **self.plot_kwargs
            )
        else:
            self.ax.scatter(
                dataset[self.axis_keys["x"]],
                dataset[self.axis_keys["y"]],
                dataset[self.axis_keys["z"]],
                **self.plot_kwargs,
            )
        for key, val in self.ax_kwargs.items():
            if "animation_value" in val:
                self.ax_kwargs[key] = eval(f"f{val}")
        self.ax.set(**self.ax_kwargs)
        # self.apply_axes_config()

    def show(self):
        plt.show()

    def save(self, filename, dpi=300):
        self.fig.savefig(filename, dpi=dpi)

    def close(self):
        plt.close(self.fig)

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
