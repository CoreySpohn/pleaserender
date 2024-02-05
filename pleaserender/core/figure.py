import datetime
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from astropy.time import Time
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

        # self.dataset = None

        self.plot_calls = {}
        self.plots_by_type = {}

        self.shared_axes = {}

    def please_add_plot(
        self,
        plot,
        row=0,
        col=0,
        rowspan=1,
        colspan=1,
        sharex_plot=None,
        sharey_plot=None,
        shared_plot_data=None,
        dataset=None,
    ):
        plot_type = plot.__class__.__name__
        if plot_type not in self.plots_by_type:
            self.plots_by_type[plot_type] = []
        self.plots_by_type[plot_type].append(self)

        plot.grid_position = (row, col, rowspan, colspan)
        self.shared_axes[plot] = {
            "x": True if sharex_plot is not None else False,
            "y": True if sharey_plot is not None else False,
            "sharex_plot": sharex_plot,
            "sharey_plot": sharey_plot,
        }
        plot.data_provided = dataset is not None
        plot.shared_plot_data = shared_plot_data
        if plot.data_provided:
            plot.data = dataset
        self.plots.append(plot)

    def please_add_subfigure(self, subfigure, row, col, rowspan=1, colspan=1):
        subfigure.grid_position = (row, col, rowspan, colspan)
        self.subfigures.append(subfigure)

    def please_set_animation_values(self, animation_values, animation_key):
        if type(animation_values) is Time:
            animation_values = animation_values.datetime64
        self.animation_values = animation_values
        self.animation_key = animation_key
        self.is_animated = len(animation_values) > 1

    def please_preview(self, frame):
        # Create figure (and all subfigures)
        self.save_path = None
        self.render_setup()
        self.render(self.animation_values[frame])
        plt.show()

    def please_render_images(self, save_path, dpi=300, file_format="png"):
        self.save_path = save_path
        self.render_setup()

        # Ensure output directory exists
        if not save_path.exists():
            save_path.mkdir(parents=True, exist_ok=True)

        # Now draw the figure
        for i, animation_value in enumerate(self.animation_values):
            self.render(animation_value)
            self.fig.savefig(Path(save_path, f"{i:003}.{file_format}"), dpi=dpi)
        self.pbar.close()

    def please_render_video(self, save_path, render_settings=None):
        self.save_path = save_path

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
        self.render_setup()

        anim = FuncAnimation(self.fig, self.render, frames=self.animation_values)
        anim.save(save_path, **save_settings)

        self.pbar.close()

    def render_setup(self):
        # Check that all the data has been generated
        self.generate_data(self.animation_values, self.animation_key)

        # Create figure (and all subfigures)
        self.fig = self.create_figure()

        # Create a progress bar
        self.pbar = tqdm(
            total=len(self.animation_values) + 1,
            desc=f"Rendering into {self.save_path}",
        )

    def render(self, animation_value):
        """
        Render the figure for the given animation value.
        """
        # Clear previous frame
        self.clear()

        # Render subfigures first
        for subfigure in self.subfigures:
            # Recursive call for subfigures
            subfigure.render(animation_value)

        # Render plots
        for plot in self.plots:
            plot.draw_plot(animation_value, self.animation_key)
            plot.adjust_settings(animation_value, self.animation_key)

        if "title" not in self.fig_kwargs:
            if isinstance(animation_value, np.datetime64):
                title = f"{Time(animation_value).decimalyear:.2f}"
            elif isinstance(animation_value, u.Quantity):
                title = (
                    f"{self.animation_key.replace('_', ' ')} "
                    "{animation_value.value:.2f}({animation_value.unit})"
                )
            elif isinstance(animation_value, float):
                title = (
                    f"{self.animation_key.replace('_', ' ')} "
                    "{animation_value.value:.2f}"
                )
            elif isinstance(animation_value, np.int64):
                title = f"{self.animation_key.replace('_', ' ')} {animation_value}"
            else:
                breakpoint()
                raise NotImplementedError(
                    f"Type {type(animation_value)} not implemented yet"
                )
        else:
            title = self.fig_kwargs["title"]
        self.fig.suptitle(title)
        if hasattr(self, "pbar"):
            self.pbar.update(1)

    def clear(self):
        # Clear plots in the subfigures first
        for subfigure in self.subfigures:
            # Recursive call for subfigures
            subfigure.clear()

        # Clear plots in the current figure
        for plot in self.plots:
            plot.ax.clear()

    def generate_data(self, animation_values, animation_key):
        for subfigure in self.subfigures:
            subfigure.generate_data(animation_values, animation_key)
        for plot in self.plots:
            if plot.shared_plot_data is not None:
                # Get the data from the other plot
                if getattr(plot.shared_plot_data, "data") is None:
                    raise ValueError("The shared plot data has not been generated yet.")
                plot.data = plot.shared_plot_data.data
            else:
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
            plot.create_axes_config()

        # Recursively create figures for each subfigure
        for subfigure in self.subfigures:
            subfigure.create_figure(
                parent_fig=self.fig,
                parent_spec=gridspec,
                animation_key=self.animation_key,
            )
        return self.fig
