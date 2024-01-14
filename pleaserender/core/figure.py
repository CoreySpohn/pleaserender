from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm


class Figure:
    def __init__(self, nrows=1, ncols=1):
        self.nrows = nrows
        self.ncols = ncols

        # List to store plots
        self.plots = []

        # List to store subfigures
        self.subfigures = []

        self.dataset = None

        self.shared_axes = {}

    def please_add_plot(
        self, plot, row, col, rowspan=1, colspan=1, sharex_plot=None, sharey_plot=None
    ):
        plot.grid_position = (row, col, rowspan, colspan)
        self.shared_axes[plot] = {
            "x": True if sharex_plot is not None else False,
            "y": True if sharey_plot is not None else False,
            "sharex_plot": sharex_plot,
            "sharey_plot": sharey_plot,
        }
        self.plots.append(plot)

    def please_add_subfigure(self, subfigure, row, col, rowspan=1, colspan=1):
        subfigure.grid_position = (row, col, rowspan, colspan)
        self.subfigures.append(subfigure)

    def please_set_animation_values(self, animation_values, animation_key):
        self.animation_values = animation_values
        self.animation_key = animation_key
        self.is_animated = len(animation_values) > 1

    def please_add_dataset(self, dataset):
        self.dataset = dataset
        for subfigure in self.subfigures:
            subfigure.please_add_dataset(dataset)

    def please_preview(self):
        # Create figure (and all subfigures)
        self.fig = self.create_figure()
        self.render(self.animation_values[-1])
        plt.show()

    def please_render(self, save_dir):
        self.render_setup()

        # Ensure output directory exists
        if not save_dir.exists():
            save_dir.mkdir(parents=True, exist_ok=True)

        # Now draw the figure
        for i, animation_value in enumerate(
            tqdm(self.animation_values, desc=f"Rendering into {save_dir}")
        ):
            self.render(animation_value)
            self.fig.savefig(Path(save_dir, f"{i:003}.png"), dpi=300)
            self.clear()

    def render_setup(self):
        if self.is_animated:
            # Check that all the data has been generated
            self.verify_figure_data(self.animation_values, self.animation_key)

        # Create figure (and all subfigures)
        self.fig = self.create_figure()

    def render(self, animation_value):
        # Render subfigures first
        for subfigure in self.subfigures:
            # Recursive call for subfigures
            subfigure.render(animation_value)

        # Render plots
        for plot in self.plots:
            plot.draw_plot(self.dataset, animation_value, self.animation_key)

    def clear(self):
        # Clear plots in the subfigures first
        for subfigure in self.subfigures:
            # Recursive call for subfigures
            subfigure.clear()

        # Clear plots in the current figure
        for plot in self.plots:
            plot.ax.clear()

    def verify_figure_data(self, animation_values, animation_key):
        for subfigure in self.subfigures:
            # Recursive call for subfigures
            subfigure.verify_figure_data(animation_values, animation_key)
        for plot in self.plots:
            necessary_data_exists = plot.verify_data(
                self.dataset, animation_values, animation_key
            )
            if not necessary_data_exists:
                plot.generate_data(animation_values, animation_key)

    def create_figure(self, parent_fig=None, parent_spec=None, animation_key=None):
        if parent_fig is None:
            # Main figure
            self.fig = plt.figure(constrained_layout=True)
        else:
            # Subfigure
            row, col, rowspan, colspan = self.grid_position
            subfigspec = parent_spec[row : row + rowspan, col : col + colspan]
            self.fig = parent_fig.add_subfigure(subfigspec)
            self.animation_key = animation_key
        gridspec = GridSpec(self.nrows, self.ncols, figure=self.fig)

        # Assign axes to plots
        for plot in self.plots:
            row, col, rowspan, colspan = plot.grid_position
            # if (self.shared_axes[plot]["x"]) or (self.shared_axes[plot]["y"]):
            plot.ax = self.fig.add_subplot(
                gridspec.new_subplotspec((row, col), rowspan=rowspan, colspan=colspan),
                sharex=self.shared_axes[plot]["sharex_plot"].ax
                if self.shared_axes[plot]["x"]
                else None,
                sharey=self.shared_axes[plot]["sharey_plot"].ax
                if self.shared_axes[plot]["y"]
                else None,
                projection=plot.projection,
            )
            plot.create_axes_config(self.dataset)

        # Recursively create figures for each subfigure
        for subfigure in self.subfigures:
            subfigure.create_figure(
                parent_fig=self.fig,
                parent_spec=gridspec,
                animation_key=self.animation_key,
            )
        return self.fig
