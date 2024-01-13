from pleaserender.core.plot import Plot


class Scatter(Plot):
    def __init__(self, **kwargs):
        """
        A generic scatter plot
        Args:
            plot_obj (Python object):
                The primary object of concern for the scatter plot, such as a
                planetary system or star that stores position data, or a dataframe
            axis_keys (dictionary):
                The keys of
        """
        super().__init__(**kwargs)
        self.plot_method = "scatter"

        self.kwargs = kwargs
