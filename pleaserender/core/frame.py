import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


class Frame:
    def __init__(self, nrows=1, ncols=1):
        self.nrows = nrows
        self.ncols = ncols

        # List to store plots
        self.plots = []

        # List to store subframes
        self.subframes = []

        self.dataset = None

    def please_add_plot(self, plot, row, col, rowspan=1, colspan=1):
        plot.grid_position = (row, col, rowspan, colspan)
        self.plots.append(plot)

    def please_add_subframe(self, subframe, row, col, rowspan=1, colspan=1):
        subframe.grid_position = (row, col, rowspan, colspan)
        self.subframes.append(subframe)

    def please_set_animation_values(self, animation_values, animation_key):
        self.animation_values = animation_values
        self.animation_key = animation_key
        self.is_animated = len(animation_values) > 1

    def please_add_dataset(self, dataset):
        self.dataset = dataset
        for subframe in self.subframes:
            subframe.please_add_dataset(dataset)

    def please_render(self):
        if self.is_animated:
            # Check that all the data has been generated
            self.verify_frame_data(self.animation_values, self.animation_key)

        # Create figure (and all subfigures)
        self.fig = self.create_figure()

        # Now draw the figure
        for i, animation_value in enumerate(self.animation_values):
            self.render(animation_value)
            self.fig.savefig(f"output/{i:003}.png", dpi=300)
            self.clear()

    def render(self, animation_value):
        # Render subframes first
        for subframe in self.subframes:
            # Recursive call for subframes
            subframe.render(animation_value)

        # Render plots
        for plot in self.plots:
            plot.draw_plot(self.dataset, animation_value, self.animation_key)

    def clear(self):
        # clear subframes first
        for subframe in self.subframes:
            # Recursive call for subframes
            subframe.clear()

        # Render plots
        for plot in self.plots:
            plot.ax.clear()

    def verify_frame_data(self, animation_values, animation_key):
        for subframe in self.subframes:
            # Recursive call for subframes
            subframe.verify_frame_data(animation_values, animation_key)
        for plot in self.plots:
            necessary_data_exists = plot.verify_data(
                self.dataset, animation_values, animation_key
            )
            if not necessary_data_exists:
                plot.generate_data(animation_values, animation_key)

    def create_figure(self, parent_fig=None, parent_spec=None, animation_key=None):
        if parent_fig is None:
            # Main frame
            self.fig = plt.figure(constrained_layout=True)
        else:
            # Subframe
            row, col, rowspan, colspan = self.grid_position
            subfigspec = parent_spec[row : row + rowspan, col : col + colspan]
            self.fig = parent_fig.add_subfigure(subfigspec)
            self.animation_key = animation_key
        gridspec = GridSpec(self.nrows, self.ncols, figure=self.fig)

        # Assign axes to plots
        for plot in self.plots:
            row, col, rowspan, colspan = plot.grid_position
            plot.ax = self.fig.add_subplot(
                gridspec.new_subplotspec((row, col), rowspan=rowspan, colspan=colspan)
            )
            plot.create_axes_config(self.dataset)

        # Recursively create figures for each subframe
        for subframe in self.subframes:
            subframe.create_figure(
                parent_fig=self.fig,
                parent_spec=gridspec,
                animation_key=self.animation_key,
            )
        return self.fig
