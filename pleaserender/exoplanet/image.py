import copy

import astropy.units as u
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto import Observation, Observations
from matplotlib.colors import LogNorm
from tqdm import tqdm

from pleaserender import util
from pleaserender.core import Plot


class Image(Plot):
    def __init__(
        self,
        observation,
        gen_data=None,
        imaging_params=None,
        **kwargs,
    ):
        default_imaging_params = {"object": "img", "plane": "coro"}
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

        super().__init__(**kwargs)

        self.gen_data = gen_data
        self.observation = observation
        self.system = observation.system
        self.coronagraph = observation.coronagraph
        self.observing_scenario = observation.observing_scenario
        self.plot_method = "imshow"
        self.name = "image"
        # self.valid_animation_keys = [
        #     "time",
        #     "central_wavelength",
        #     "frame",
        #     "spectral_wavelength(nm)",
        # ]

    # def generate_data_old(self, observations):
    #     obs_ds = Observations.run(self, observations=observations)
    #     return obs_ds

    def generate_data(self):
        assert self.gen_data is not None, (
            "The gen_data dictionary must be provided if"
            " the data is not already provided."
        )
        times = self.observation.time.reshape(1)
        wavelengths = self.observation.central_wavelength.reshape(1)
        has_time = "time" in self.gen_data
        has_wavelength = "central_wavelength" in self.gen_data
        assert (has_time and not has_wavelength) or (not has_time and has_wavelength), (
            "The generation_data dictionary must contain either a 'time'"
            " or 'central_wavelength' key."
        )
        if "time" in self.gen_data:
            times = self.gen_data["time"]
        elif "central_wavelength" in self.gen_data:
            wavelengths = self.gen_data["central_wavelength"]
        all_obs = Observations(self.observation, times, wavelengths)
        # self.observations = set(all_obs.create_observations())
        self.data = all_obs.run()

    def get_required_keys(self):
        """
        Get the keys that are required to plot the data
        """
        self.get_required_keys = []
        all_possible_keys = [
            "time",
            "central_wavelength",
            "frame",
            "spectral_wavelength(nm)",
        ]
        for key in all_possible_keys:
            if key in self.data.dims:
                self.required_keys.append(key)

    def compile_necessary_info(self, animation_values, animation_key, plot_calls):
        base_observation = Observation(
            self.coronagraph,
            self.system,
            self.observing_scenario,
            logging_level="WARNING",
        )
        if animation_key == "time":
            times = Time(animation_values)
            wavelengths = base_observation.central_wavelength.reshape(1)
        elif animation_key == "central_wavelength":
            times = base_observation.time.reshape(1)
            wavelengths = animation_values
        _observations = Observations(base_observation, times, wavelengths)
        observations = set(_observations.create_observations())
        if Image in plot_calls:
            plot_calls[Image]["args"][0].update(observations)
        else:
            plot_calls[Image] = {"args": [observations], "kwargs": {}}
        return plot_calls

    def get_observation_object(self, animation_value, animation_key):
        obs_scen = copy.deepcopy(self.observing_scenario)
        if animation_key == "time":
            obs_scen.scenario["time"] = Time(animation_value)
        elif animation_key == "central_wavelength":
            obs_scen.scenario["central_wavelength"] = animation_value
        obs = Observation(
            self.coronagraph, self.system, obs_scen, logging_level="WARNING"
        )
        return obs

    def draw_plot(self, plot_kwargs=None):
        photons = self.get_frame_data()
        if plot_kwargs is None:
            plot_kwargs = self.plot_kwargs
        self.ax.imshow(
            photons, origin="lower", cmap="viridis", norm=LogNorm(vmin=0, vmax=1e5)
        )

    def get_frame_data(self):
        # draw_data = self.state.next_frame_values
        # if "time" in draw_data.keys():
        #     animation_key = "time"
        #     animation_value = draw_data["time"]
        # elif "central_wavelength" in draw_data.keys():
        #     animation_key = "central_wavelength"
        #     animation_value = draw_data["central_wavelength"]
        # obs = self.get_observation_object(animation_value, animation_key)

        # necessary_dims = self.data.attrs["dist_attrs_for_images"]
        # sel_call = self.create_sel_call(obs, necessary_dims)
        # sel_call = self.add_extra_sel_call(sel_call, obs, draw_data)
        sel_call = copy.copy(self.state.next_sel)
        for key, val in sel_call.items():
            if key in self.cumulative_keys:
                sel_call[key] = slice(None, val)

        base_data = self.data.sel(**sel_call)
        photons = self.process_photons(base_data)
        # photons = self.data.sel(**sel_call)[self.imsel].data
        return photons

    def create_sel_call(self, obs, necessary_dims):
        sel_call = {}
        if self.imaging_params["plane"] == "coro":
            sel_call["xpix(coro)"] = np.arange(obs.coronagraph.npixels)
            sel_call["ypix(coro)"] = np.arange(obs.coronagraph.npixels)
        elif self.imaging_params["plane"] == "det":
            sel_call["xpix(det)"] = np.arange(obs.detector_shape[0])
            sel_call["ypix(det)"] = np.arange(obs.detector_shape[1])
        else:
            raise NotImplementedError("Only coro and det planes are implemented.")
        for dim in necessary_dims:
            val = getattr(obs, dim)
            if isinstance(val, u.Quantity):
                val = val.value
            elif isinstance(val, Time):
                val = val.datetime64
            sel_call[dim] = val
        return sel_call

    def process_photons(self, photons):
        return photons[self.imsel].data

    def add_extra_sel_call(self, sel_call, obs, animation_value, animation_key):
        """
        Add any additional selection calls to the sel_call dictionary. This method
        should be overridden by subclasses to add any additional selection calls
        """
        return sel_call

    def ax_lims_helper(self, necessary_axes, data=None, equal=False):
        if self.imaging_params["plane"] == "coro":
            xrange = [0, self.coronagraph.npixels]
            yrange = [0, self.coronagraph.npixels]
        elif self.imaging_params["plane"] == "det":
            det_shape = self.observing_scenario.scenario["detector_shape"]
            xrange = [0, det_shape[0]]
            yrange = [0, det_shape[1]]
        self.ax_kwargs["lims"]["x"] = util.calculate_axis_limits_and_ticks(
            *xrange, exact=True
        )
        self.ax_kwargs["lims"]["y"] = util.calculate_axis_limits_and_ticks(
            *yrange, exact=True
        )

    def auto_set_title(self):
        title = self.create_auto_title()
        self.ax.set_title(title[:-2])

    def create_auto_title(self):
        title = ""
        for key in self.required_keys:
            if key not in self.cumulative_keys:
                val = self.state.next_sel[key]
                title += util.create_val_str(key, val)
                title += ", "
        return title[:-2]
