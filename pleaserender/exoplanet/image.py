import copy

import astropy.units as u
import coronagraphoto.util as cu
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto import Observation, Observations
from exoverses.util import misc
from lod_unit import lod, lod_eq
from matplotlib.colors import LogNorm, Normalize
from tqdm import tqdm

from pleaserender import util
from pleaserender.core import Plot


class Image(Plot):
    def __init__(self, observation, gen_data=None, imaging_params=None, **kwargs):
        default_imaging_params = {"object": "img", "plane": "coro", "unit": u.pix}
        self.imaging_params = default_imaging_params
        if imaging_params is not None:
            self.imaging_params.update(imaging_params)
        self.imsel = f"{self.imaging_params['object']}({self.imaging_params['plane']})"

        if kwargs.get("axis_keys") is None:
            kwargs["axis_keys"] = {
                "x": f"xpix({self.imaging_params['plane']})",
                "y": f"ypix({self.imaging_params['plane']})",
            }
        if not kwargs.get("ax_kwargs"):
            kwargs["ax_kwargs"] = {}
        if "auto_title" not in kwargs["ax_kwargs"]:
            kwargs["ax_kwargs"]["auto_title"] = True
        if "aspect" not in kwargs["ax_kwargs"]:
            kwargs["ax_kwargs"]["aspect"] = "equal"

        # Set the default plot_kwargs
        if not kwargs.get("plot_kwargs"):
            kwargs["plot_kwargs"] = {}
        if "cmap" not in kwargs["plot_kwargs"]:
            kwargs["plot_kwargs"]["cmap"] = "viridis"
        if "norm" not in kwargs["plot_kwargs"]:
            kwargs["plot_kwargs"]["norm"] = Normalize(vmin=0, vmax=4e6)
        if "origin" not in kwargs["plot_kwargs"]:
            kwargs["plot_kwargs"]["origin"] = "lower"

        super().__init__(**kwargs)

        self.gen_data = gen_data
        self.observation = observation
        self.system = observation.system
        self.coronagraph = observation.coronagraph
        self.observing_scenario = observation.scenario
        self.settings = observation.settings
        self.plot_method = "imshow"
        self.name = "image"

    def generate_data(self):
        assert self.gen_data is not None, (
            "The gen_data dictionary must be provided if"
            " the data is not already provided."
        )
        times = self.observation.scenario.start_time.reshape(1)
        wavelengths = self.observation.scenario.central_wavelength.reshape(1)
        has_start_time = "start_time" in self.gen_data
        has_wavelength = "central_wavelength" in self.gen_data
        assert (has_start_time and not has_wavelength) or (
            not has_start_time and has_wavelength
        ), (
            "The generation_data dictionary must contain either a 'time'"
            " or 'central_wavelength' key."
        )
        if "start_time" in self.gen_data:
            times = self.gen_data["start_time"]
        elif "central_wavelength" in self.gen_data:
            wavelengths = self.gen_data["central_wavelength"]
        all_obs = Observations(self.observation, times, wavelengths)
        self.data = all_obs.run()

    def get_required_keys(self):
        """
        Get the keys that are required to plot the data
        """
        all_possible_keys = [
            "start_time",
            "time",
            "central_wavelength",
            "spectral_wavelength(nm)",
        ]
        for key in all_possible_keys:
            if key in self.data.dims:
                self.required_keys.append(key)
        sum_time = (
            "time" in self.access_kwargs["sum_keys"]
            and "start_time" in self.required_keys
        )
        if "start_time" in self.required_keys and "time" in self.required_keys:
            self.skipna = True
        else:
            self.skipna = False

    def get_observation_object(self, animation_value, animation_key):
        obs_scen = copy.deepcopy(self.observing_scenario)
        if animation_key == "start_time":
            obs_scen.scenario["start_time"] = Time(animation_value)
        elif animation_key == "central_wavelength":
            obs_scen.scenario["central_wavelength"] = animation_value
        obs = Observation(
            self.coronagraph,
            self.system,
            obs_scen,
            self.settings,
            logging_level="WARNING",
        )
        return obs

    def draw_plot(self, plot_kwargs=None):
        # photons = self.get_frame_data()
        if plot_kwargs is None:
            plot_kwargs = self.plot_kwargs
        data = self.get_plot_data()
        extent = self.get_extent()
        self.ax.imshow(data, extent=extent, **plot_kwargs)

    def get_plot_data(self, context=None):
        """
        Method to get the data required to plot. Not implemented in the core classes.
        """
        sel = self.state.context_sel(context)
        plot_data = self.data.sel(**sel)

        # Sum over keys if necessary
        for key in self.access_kwargs["sum_keys"]:
            if key in plot_data.dims and len(plot_data[key].shape) > 0:
                plot_data = plot_data.sum(dim=key, skipna=self.skipna)
        data = plot_data[self.imsel].data
        return data

    def check_valid_context(self, fig_context):
        context = self.state.extract_plot_context(fig_context)
        if context is None:
            return False
        for key, val in context.items():
            if val not in self.data[key]:
                return False
        plot_data = self.get_plot_data(context)
        valid = np.any(~np.isnan(plot_data))
        # Check if the start time is after the time
        if "start_time" in context and "time" in context:
            exposure_start_time = Time(context["start_time"])
            exposure_end_time = (
                exposure_start_time + self.observation.scenario.exposure_time
            )
            valid_time = exposure_start_time <= context["time"] <= exposure_end_time
            valid = valid and valid_time
        return valid

    def add_dependent_context(self, context):
        if "start_time" in self.required_keys and "time" in context:
            obs_time = context["time"]
            if obs_time not in self.data["time"]:
                # If the observation time is not in the data, then the
                # data has no corresponding start time
                return context
            # Get the start times based on the given context
            start_times = self.data.sel(**context)["start_time"].values
            after_start = start_times <= obs_time

            # Get the end times based on the exposure time
            end_times = Time(start_times) + self.observation.scenario.exposure_time
            before_end = obs_time < end_times

            # Find the start time that is after the observation time
            start_time = Time(start_times[after_start & before_end]).datetime64
            assert len(start_time) != 0, "No valid start times found"
            assert len(start_time) < 2, (
                "Multiple valid start times found. Check that the exposure time is less"
                "than the time between observations"
            )
            context["start_time"] = start_time[0]
        return context

    def access_data(self, data):
        return data[self.imsel].data

    def ax_lims_helper(self, necessary_axes, data=None, equal=False):
        extent = self.get_extent()
        xrange = [extent[0], extent[1]]
        yrange = [extent[2], extent[3]]
        self.ax_kwargs["lims"]["x"] = util.calculate_axis_limits_and_ticks(*xrange)
        self.ax_kwargs["lims"]["y"] = util.calculate_axis_limits_and_ticks(*yrange)

    def get_extent(self):
        if self.imaging_params["unit"] == u.pix:
            if self.imaging_params["plane"] == "coro":
                extent = [0, self.coronagraph.npixels, 0, self.coronagraph.npixels]
            elif self.imaging_params["plane"] == "det":
                det_shape = self.observing_scenario.scenario["detector_shape"]
                extent = [0, det_shape[0], 0, det_shape[1]]
        else:
            xarr, yarr = cu.convert_pixels(
                self.imaging_params["unit"], self.observation, plane="coro"
            )
            xarr = xarr.value
            yarr = yarr.value
            extent = [xarr.min(), xarr.max(), yarr.min(), yarr.max()]
        return extent
