from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation
from matplotlib.gridspec import GridSpec
from tqdm import tqdm


class Figure:
    def __init__(self, nrows=1, ncols=1, fig_kwargs=None, gs_kwargs=None):
        self.nrows = nrows
        self.ncols = ncols

        # Default figure kwargs
        default_fig_kwargs = {"figsize": (5, 10), "layout": None}
        self.fig_kwargs = default_fig_kwargs
        if fig_kwargs is not None:
            self.fig_kwargs.update(fig_kwargs)

        # Default gridspec kwargs
        default_gs_kwargs = {"height_ratios": [1] * nrows, "width_ratios": [1] * ncols}
        self.gs_kwargs = default_gs_kwargs
        if gs_kwargs is not None:
            self.gs_kwargs.update(gs_kwargs)

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

    def please_preview(self, frame):
        # Create figure (and all subfigures)
        self.fig = self.create_figure()
        self.render(self.animation_values[frame])
        plt.show()

    def please_render_images(self, save_path, dpi=300, file_format="png"):
        self.render_setup(save_path)

        # Ensure output directory exists
        if not save_path.exists():
            save_path.mkdir(parents=True, exist_ok=True)

        # Now draw the figure
        for i, animation_value in enumerate(self.animation_values):
            self.render(animation_value)
            self.fig.savefig(Path(save_path, f"{i:003}.{file_format}"), dpi=dpi)
            self.clear()
        self.pbar.close()

    def please_render_video(self, save_path, render_settings=None):
        # Set up settings
        final_render_settings = {
            "animation_duration": 10,
            "codec": "h264",
            "bitrate": -1,
            "extra_args": ["-vcodec", "libx264", "-pix_fmt", "yuv420p"],
            "img_dpi": 300,
        }
        if render_settings is not None:
            final_render_settings.update(render_settings)
        final_render_settings["framerate"] = (
            len(self.animation_values) / final_render_settings["animation_duration"]
        )

        # Create writer class for the animation
        writer = FFMpegWriter(
            fps=final_render_settings["framerate"],
            codec=final_render_settings["codec"],
            bitrate=final_render_settings["bitrate"],
            extra_args=final_render_settings["extra_args"],
        )

        save_settings = {
            "dpi": final_render_settings["img_dpi"],
            "writer": writer,
        }

        # Create the figure (and all subfigures)
        self.render_setup(save_path)

        anim = FuncAnimation(self.fig, self.render, frames=self.animation_values)
        anim.save(save_path, **save_settings)

        self.pbar.close()

    def render_setup(self, save_path):
        if self.is_animated:
            # Check that all the data has been generated
            self.verify_figure_data(self.animation_values, self.animation_key)

        # Create figure (and all subfigures)
        self.fig = self.create_figure()

        # Create a progress bar
        self.pbar = tqdm(
            total=len(self.animation_values) + 1, desc=f"Rendering into {save_path}"
        )

    def render(self, animation_value):
        """
        Render the figure for the given animation value.
        """
        # Render subfigures first
        for subfigure in self.subfigures:
            # Recursive call for subfigures
            subfigure.render(animation_value)

        # Render plots
        for plot in self.plots:
            plot.draw_plot(self.dataset, animation_value, self.animation_key)
            plot.adjust_settings(self.dataset, animation_value, self.animation_key)

        self.pbar.update(1)

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
            self.fig = plt.figure(**self.fig_kwargs)
        else:
            # Subfigure
            row, col, rowspan, colspan = self.grid_position
            subfigspec = parent_spec[row : row + rowspan, col : col + colspan]
            self.fig = parent_fig.add_subfigure(subfigspec)
            self.animation_key = animation_key
        gridspec = GridSpec(self.nrows, self.ncols, figure=self.fig, **self.gs_kwargs)

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
