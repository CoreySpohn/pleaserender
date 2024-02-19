import copy
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
        self.primary_plots = []

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
        primary=True,
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
        if plot.data_provided:
            plot.data = dataset
        plot.shared_plot_data = shared_plot_data
        self.plots.append(plot)
        if primary:
            self.primary_plots.append(plot)

    def please_add_subfigure(self, subfigure, row, col, rowspan=1, colspan=1):
        subfigure.grid_position = (row, col, rowspan, colspan)
        self.subfigures.append(subfigure)

    def please_set_animation_levels(self, levels):
        self.levels = levels
        self.key_levels = {}
        for level, level_keys in levels.items():
            for level_key in level_keys:
                self.key_levels[level_key] = level
        # default_animation_info = {"method": "index", "initial": 0}
        # self.animation_info = {}
        # for key, value in animation_info.items():
        #     self.animation_info[key] = copy.copy(default_animation_info)
        #     if value is not None:
        #         self.animation_info[key].update(value)

        # # Animation order defines how we are going to iterate over the data
        # self.animation_order = []
        # for key in animation_info.keys():
        #     self.animation_order.append(key)

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
            "framerate": 30,
        }
        if render_settings is not None:
            final_render_settings.update(render_settings)
        # final_render_settings["framerate"] = (
        #     len(self.animation_values) / final_render_settings["animation_duration"]
        # )

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
        self.plots = self.collect_plots(self.plots)
        self.subfigures = self.collect_subfigures(self.subfigures)
        # times = Time(np.linspace(2000, 2001, 4), format="decimalyear")
        # levels = {0: ["time"], 1: ["frame"]}
        # keys = ["time", "frame"]

        # self.state = AnimationRenderState(levels)
        # state1 = PlotRenderState(keys, levels, {"time": times[0], "frame": 0})
        # state2 = PlotRenderState(keys, levels, {"time": times[0], "frame": 1})
        # state3 = PlotRenderState(keys, levels, {"time": times[1], "frame": 0})
        # state4 = PlotRenderState(["time"], levels, {"time": times[0]})

        with writer.saving(self.fig, save_path, 300):
            self.render(writer)
        # with writer.saving(self.fig, save_path, len(self.animation_values)):
        # anim = FuncAnimation(self.fig, self.render, frames=self.animation_values)
        # anim.save(save_path, **save_settings)

        # self.pbar.close()

    def collect_plots(self, plots):
        for subfigure in self.subfigures:
            _plots = subfigure.collect_plots(plots)
            plots.extend(_plots)
        return plots

    def collect_subfigures(self, subfigures):
        for subfigure in self.subfigures:
            _subfigures = subfigure.collect_subfigures(subfigures)
            subfigures.extend(_subfigures)
        return subfigures

    def render_setup(self):
        # Check that all the data has been generated
        self.generate_data()

        # Create figure (and all subfigures)
        self.fig = self.create_figure()

        # Create a progress bar
        # self.pbar = tqdm(
        #     total=len(self.animation_values) + 1,
        #     desc=f"Rendering into {self.save_path}",
        # )

        # Set the initial state
        # self.render_state = {"nframes": 0, "draw_data": {}}
        # for i, animation_key in enumerate(self.animation_order):
        #     self.render_state[animation_key] = {}
        #     self.render_state[animation_key]["order"] = i
        #     self.render_state[animation_key]["method"] = self.animation_info[
        #         animation_key
        #     ].get("method")
        #     initial = self.animation_info[animation_key]["initial"]
        #     self.render_state["draw_data"][animation_key] = initial
        #     # if self.render_state[animation_key]["method"] == "index":
        #     #     self.render_state[animation_key]["current_ind"] = initial
        #     # elif self.render_state[animation_key]["method"] == "value":
        #     #     self.render_state[animation_key]["current_val"] = initial
        #     # else:
        #     #     raise NotImplementedError("Only index and value supported")

        #     # Progress bar for this key
        #     if animation_key in self.generation_data.keys():
        #         # self.render_state[animation_key]["est_bar"] = False
        #         ndata = len(self.generation_data[animation_key])
        #         self.render_state[animation_key]["pbar"] = tqdm(
        #             total=ndata,
        #             desc=f"Rendering into {self.save_path}. {animation_key}",
        #         )
        # else:
        #     self.render_state[animation_key]["pbar"] = None
        #     self.render_state[animation_key]["est_bar"] = True
        #     self.render_state[animation_key]["est_max"] = 0

    def render(self, writer):
        # Setting up a loop to render the frames until all are finished
        states = [plot.state for plot in self.plots]
        self.finished = False
        # next_state = sorted(states)[0]
        # Context refers to the current state of the animation
        # for each level value
        while True:
            self.context = self.get_next_context()
            # print(self.context)
            if self.finished:
                return

            for plot in self.plots:
                plot.render(self.context)
            writer.grab_frame()
            # print("\n")

            # Render the current frame
            # next_state.plot.ax.clear()
            # next_state.plot.draw_plot()

            # min_level = min(next_state.levels.keys())
            # min_level_key = next_state.levels[min_level]
            # current_val = next_state.next_sel[min_level_key]
            # next_state.plot.adjust_settings(current_val, min_level_key)

            # state_order = sorted(states)
            # Check if any states are equal to the state just drawn
            # any_left_to_draw = any(
            #     [other_state == state_order[0] for other_state in state_order[1:]]
            # )
            # any_parents_left_to_draw = any(
            #     [
            #         state_order[0].is_parent(other_state)
            #         for other_state in state_order[1:]
            #     ]
            # )

            # Iterate the state we just drew
            # next_state.iterate_sel()
            # print(next_state)

            # Get the next state to draw
            # next_state = sorted(states)[0]
            # all_finished = all([state.finished for state in states])

    def get_next_context(self):
        all_key_vals, ignored_keys = self.generate_all_key_vals()
        if not hasattr(self, "context"):
            # Dynamically generate all_key_vals and initialize context if it's
            # the first call
            self.context = self.get_first_context(all_key_vals)
            if self.check_valid_context_for_any_plot(self.context):
                # First context is valid, return
                return self.context

        first_context = self.get_first_context(all_key_vals)
        # Save the original context to detect a full cycle
        while True:
            for _, keys in sorted(self.levels.items(), reverse=True):
                for key in keys:
                    if key in ignored_keys:
                        # Skip iteration for keys with no values
                        continue
                    current_index = all_key_vals[key].index(self.context[key])
                    if current_index + 1 < len(all_key_vals[key]):
                        key_val = all_key_vals[key][current_index + 1]
                        self.context[key] = key_val
                        print(f"Checking {key}:{key_val} from Figure")
                        if self.check_valid_context_for_any_plot(self.context):
                            # Valid context found
                            return self.context
                        # Break out of the inner loop to reset context in outer
                        # loop if needed
                        break
                    else:
                        self.context[key] = all_key_vals[key][0]
                else:
                    # Continue if the inner loop wasn't broken
                    continue
                # Inner loop was broken, break the outer
                break

            # Check if we have looped back to the original context without
            # finding a valid one
            if self.context == first_context:
                self.finished = True
                return
                # raise StopIteration("All possible contexts have been exhausted.")

    def generate_all_key_vals(self):
        all_key_vals = {}
        ignored_keys = []
        for _, level_keys in self.levels.items():
            for level_key in level_keys:
                _key_vals = set()
                for plot in self.plots:
                    if level_key in plot.required_keys:
                        _key_vals.update(plot.state.key_values[level_key])
                if len(_key_vals) == 0:
                    ignored_keys.append(level_key)
                else:
                    all_key_vals[level_key] = sorted(_key_vals)
        return all_key_vals, ignored_keys

    def get_first_context(self, all_key_vals):
        # Initialize context with the first values from the dynamically
        # generated all_key_vals
        context = {k: v[0] for k, v in all_key_vals.items()}
        return context

    def check_valid_context_for_any_plot(self, context):
        # Check if the context is valid for at least one plot
        valid = any(plot.check_valid_context(context) for plot in self.primary_plots)
        # print(f"Figure: {valid}")
        return valid

    def clear(self):
        # Clear plots in the subfigures first
        for subfigure in self.subfigures:
            # Recursive call for subfigures
            subfigure.clear()

        # Clear plots in the current figure
        for plot in self.plots:
            plot.ax.clear()

    def generate_data(self):
        for subfigure in self.subfigures:
            subfigure.generate_data()
        for plot in self.plots:
            if plot.shared_plot_data is not None:
                # Get the data from the other plot
                if getattr(plot.shared_plot_data, "data") is None:
                    raise ValueError("The shared plot data has not been generated yet.")
                plot.data = plot.shared_plot_data.data
            else:
                plot.generate_data()
            plot.get_required_keys()
            plot.create_render_state(self.levels)

    def create_figure(self, parent_fig=None, parent_spec=None, animation_key=None):
        if parent_fig is None:
            # Main figure
            self.fig = plt.figure(**self.fig_kwargs)
        else:
            # Subfigure
            row, col, rowspan, colspan = self.grid_position
            subfigspec = parent_spec[row : row + rowspan, col : col + colspan]
            self.fig = parent_fig.add_subfigure(subfigspec)
            self.primary_animation_key = animation_key
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
                animation_key=self.primary_animation_key,
            )
        return self.fig

    def create_title(self, animation_value, animation_key):
        if "title" not in self.fig_kwargs:
            if isinstance(animation_value, np.datetime64):
                title = f"{Time(animation_value).decimalyear:.2f}"
            elif isinstance(animation_value, u.Quantity):
                title = (
                    f"{animation_key.replace('_', ' ')} "
                    "{animation_value.value:.2f}({animation_value.unit})"
                )
            elif isinstance(animation_value, float):
                title = (
                    f"{animation_key.replace('_', ' ')} " "{animation_value.value:.2f}"
                )
            elif isinstance(animation_value, np.int64):
                title = f"{animation_key.replace('_', ' ')} {animation_value}"
            else:
                breakpoint()
                raise NotImplementedError(
                    f"Type {type(animation_value)} not implemented yet"
                )
        else:
            title = self.fig_kwargs["title"]
        return title
