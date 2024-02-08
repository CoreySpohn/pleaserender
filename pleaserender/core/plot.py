import copy
import inspect

import numpy as np

import pleaserender.util as util


class Plot:
    def __init__(
        self,
        draw_key,
        axis_keys=None,
        plot_kwargs=None,
        ax_kwargs=None,
        animation_kwargs=None,
    ):
        self.ax = None
        self.draw_key = draw_key
        self.required_keys = [draw_key]

        # self.animation_style = animation_style
        self.plot_method = "plot"

        default_axis_keys = {"x": "x", "y": "y", "z": None}
        self.axis_keys = default_axis_keys
        if axis_keys is not None:
            self.axis_keys.update(axis_keys)
        self.given_axis_keys = [
            key for key, val in self.axis_keys.items() if val is not None
        ]
        self.given_axis_names = [
            val for val in self.axis_keys.values() if val is not None
        ]
        self.is_3d = len(self.given_axis_keys) == 3
        self.projection = "3d" if self.is_3d else None

        # Set the default plot_kwargs, then update with any user-specified
        # kwargs. The plot_kwargs are used in the plot call.
        default_plot_kwargs = {}
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
        if self.is_3d:
            default_ax_kwargs["zlabel"] = self.axis_keys.get("z")
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

    def draw_plot(self, animation_value, animation_key, plot_kwargs=None):
        """
        Method to draw the plot. Nothing past the generic plotting in this base class
        """
        self.generic_plot(
            self.data, animation_value, animation_key, plot_kwargs=plot_kwargs
        )

    def generic_plot(
        self,
        data,
        animation_value,
        animation_key,
        axis_keys=None,
        plot_method=None,
        plot_kwargs=None,
        animation_kwargs=None,
    ):
        # Get the correct plotting method, I think this is an over-engineered way to
        # call either ax.plot or ax.scatter. Might work for other things, idk
        if plot_method is None:
            plot_method = self.plot_method
        plot_method = getattr(self.ax, plot_method)
        if animation_kwargs is None:
            animation_kwargs = self.animation_kwargs
        data = util.filter_data(
            data,
            animation_value,
            animation_key,
            animation_kwargs["animation_style"],
        )
        ax_keys = self.given_axis_keys
        if axis_keys is not None:
            ax_keys = axis_keys

        # Allows for 2 and 3 dimensional data with the same call
        separated_data = [data[self.axis_keys[axis_key]].values for axis_key in ax_keys]
        if plot_kwargs is None:
            plot_kwargs = self.plot_kwargs
        plot_method(*separated_data, **plot_kwargs)

    def adjust_settings(self, animation_value, animation_key):
        self.generic_settings(animation_value, animation_key)

    def generic_settings(self, animation_value, animation_key):
        # Evaulate any f-strings provided in ax_kwargs
        original_kwargs = {}
        for key, val in self.ax_kwargs.items():
            if type(val) is str:
                if "animation_value" in val:
                    original_kwargs[key] = val
                    self.ax_kwargs[key] = eval(f"f{val}")

        # Get the keys that are settable with ax.set
        ax_set_keys = inspect.signature(self.ax.set).parameters.keys()
        settable_kwargs = {
            key: val for key, val in self.ax_kwargs.items() if key in ax_set_keys
        }
        # Set all the axes properties based on ax_kwargs
        self.ax.set(**settable_kwargs)

        unsettable_kwargs = {
            key: val for key, val in self.ax_kwargs.items() if key not in ax_set_keys
        }
        self.handle_unsettable_kwargs(unsettable_kwargs)

        # replace original kwargs
        self.ax_kwargs.update(original_kwargs)

        if self.is_3d:
            n_vals = self.data[animation_key].size
            current_index = np.argmax(
                self.data[animation_key].values == animation_value
            )
            frame_view = self.calc_frame_view(current_index, n_vals)
            self.ax.view_init(**frame_view)

    def calc_frame_view(self, current_ind, n_vals):
        if self.is_3d:
            frame_view = {
                "elev": self.animation_kwargs["elev"],
                "azim": self.animation_kwargs["azim"],
                "roll": self.animation_kwargs["roll"],
            }
            if self.animation_kwargs["rotate"] is not None:
                for key, val in self.animation_kwargs["rotate"].items():
                    frame_view[key] = np.linspace(*val, n_vals)[current_ind]
        else:
            plane = [
                axis for axis in ["x", "y", "z"] if axis not in self.given_axis_keys
            ]
            frame_view = {"azim": 0, "elev": 0, "roll": 0}
            if plane == "x":
                pass
            elif plane == "y":
                frame_view["azim"] = 270
            else:
                frame_view["azim"] = 270
                frame_view["elev"] = 90
        return frame_view

    def verify_data(self, animation_values, animation_key):
        """
        Method to check if all the data required for the plot exists.
        Args:
            # dataset (xr.Dataset):
            #     Dataset containing the data and coordinates
            animation_values (numpy.ndarray):
                List of values for the animation
            animation_key (str):
                Key for the animation values in the plot object's dataframe
        """
        animation_data = self.data[animation_key].values

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

    def create_axes_config(self):
        """
        Method that fills ax_kwargs with defaults for things not already
        specified.
        Currently sets the axis limits.
        """
        # Set the axis limits if they are not specified
        self.handle_axes_limits_and_ticks(self.ax_kwargs.get("equal_lims"))

    def handle_axes_limits_and_ticks(self, data=None, equal=False):
        necessary_axes = ["x", "y"]
        if self.axis_keys.get("z") is not None:
            necessary_axes.append("z")

        using_set = [False] * len(necessary_axes)
        for i, ax_letter in enumerate(necessary_axes):
            if self.ax_kwargs.get(f"{ax_letter}lim") is not None:
                using_set[i] = True
        using_unset = [not val for val in using_set]
        if any(using_unset):
            lims_exist = self.ax_kwargs.get("lims") is not None
            if not lims_exist:
                self.ax_kwargs["lims"] = {}
        necessary_axes = np.array(necessary_axes)[using_unset]
        self.ax_lims_helper(necessary_axes, data=data, equal=equal)

    def ax_lims_helper(self, necessary_axes, data=None, equal=False):
        if data is None:
            data = self.data
        if equal:
            max_val = 0
            for ax_letter in necessary_axes:
                if ax_lims := self.ax_kwargs["lims"].get(ax_letter):
                    if len(ax_lims) == 2:
                        ax_lims = util.calculate_axis_limits_and_ticks(
                            *ax_lims, exact=True
                        )
                else:
                    ax_data = data[self.axis_keys.get(ax_letter)].values
                    max_val = max(max_val, np.abs(ax_data.max()))
                    max_val = max(max_val, np.abs(ax_data.min()))
                    ax_lims = util.calculate_axis_limits_and_ticks(-max_val, max_val)
            for ax_letter in necessary_axes:
                self.ax_kwargs["lims"][ax_letter] = ax_lims
        else:
            for ax_letter in necessary_axes:
                axlims = util.calculate_axis_limits_and_ticks(
                    data[self.axis_keys.get(ax_letter)].min(),
                    data[self.axis_keys.get(ax_letter)].max(),
                )

                self.ax_kwargs["lims"][ax_letter] = axlims

    def handle_unsettable_kwargs(self, kwargs):
        if "lims" in kwargs:
            # lims is a tuple of (val0, valf, dval, offset)
            # Check if the lims are already set
            for ax_letter in self.given_axis_keys:
                if self.ax_kwargs.get(f"{ax_letter}lim") is None:
                    axlims = kwargs["lims"][ax_letter]
                    use_minor = not self.is_3d
                    self.set_lims_and_ticks(*axlims, use_minor=use_minor)
                    # assert (
                    #     self.ax_kwargs.get(f"{ax_letter}lim") is None
                    # ), f"lims and {ax_letter}lim cannot both be set."

    def set_lims_and_ticks(self, val0, valf, dval, offset, use_minor=True):
        for ax_letter in self.given_axis_keys:
            set_lim = getattr(self.ax, f"set_{ax_letter}lim")
            set_ticks = getattr(self.ax, f"set_{ax_letter}ticks")
            set_lim([val0 - offset, valf + offset])
            set_ticks(np.arange(val0, valf + dval / 4, dval))
            if use_minor:
                set_ticks(np.arange(val0 + dval / 2, valf + dval / 4, dval), minor=True)

    def project_point(
        self,
        data,
        animation_value,
        animation_key,
        ax_letter,
        object,
        projection_kwargs=None,
    ):
        min_value = self.ax_kwargs["lims"][ax_letter][0]
        _args = (data, animation_value, animation_key)
        _axis_keys = ["x", "y", "z"]
        _axis_keys.remove(ax_letter)
        axis_keys = {key: key for key in _axis_keys}
        default_projection_kwargs = {"zs": min_value, "zdir": ax_letter, "alpha": 0.1}
        if projection_kwargs is not None:
            default_projection_kwargs.update(projection_kwargs)
        # Plot the planet's position with a marker
        position_kwargs = copy.deepcopy(object.plot_kwargs)
        position_kwargs.update(default_projection_kwargs)
        position_kwargs["facecolor"] = "grey"
        position_kwargs["edgecolor"] = "grey"
        self.generic_plot(
            *_args,
            axis_keys=axis_keys,
            plot_method="scatter",
            animation_kwargs={"animation_style": "Single point"},
            plot_kwargs=position_kwargs,
        )

    def project_trail(
        self,
        data,
        animation_value,
        animation_key,
        ax_letter,
        projection_kwargs=None,
    ):
        min_value = self.ax_kwargs["lims"][ax_letter][0]
        _args = (data, animation_value, animation_key)
        _axis_keys = ["x", "y", "z"]
        _axis_keys.remove(ax_letter)
        axis_keys = {key: key for key in _axis_keys}
        default_projection_kwargs = {"zs": min_value, "zdir": ax_letter, "alpha": 0.1}
        if projection_kwargs is not None:
            default_projection_kwargs.update(projection_kwargs)
        # Project the planet's trail
        projection_trail_kwargs = copy.deepcopy(default_projection_kwargs)
        projection_trail_kwargs["linestyle"] = "-"
        projection_trail_kwargs["color"] = "w"
        self.generic_plot(
            *_args,
            axis_keys=axis_keys,
            plot_method="plot",
            animation_kwargs={"animation_style": "Cumulative"},
            plot_kwargs=projection_trail_kwargs,
        )
