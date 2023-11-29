from mpl_toolkits.mplot3d import Axes3D

from pleaserender.core.plot import Plot


class Scatter(Plot):
    def __init__(self, plot_obj, **kwargs):
        """
        A generic scatter plot
        Args:
            plot_obj (Python object):
                The primary object of concern for the scatter plot, such as a
                planetary system or star that stores position data, or a dataframe
            axis_keys (dictionary):
                The keys of
        """
        super().__init__(plot_obj, **kwargs)

        self.kwargs = kwargs

    def draw_plot(self, animation_value, animation_key):
        if self.axis_keys["z"] is None:
            data = self.plot_obj.filter_data(
                animation_value, animation_key, self.animation_style
            )
            self.ax.scatter(
                data[self.axis_keys["x"]], data[self.axis_keys["y"]], **self.plot_kwargs
            )
        else:
            self.ax.scatter(
                self.plot_obj[self.axis_keys["x"]],
                self.plot_obj[self.axis_keys["y"]],
                self.plot_obj[self.axis_keys["z"]],
                **self.plot_kwargs
            )
        self.apply_axes_config()
