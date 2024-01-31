import copy

import astropy.units as u
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto import Observation, Observations
from tqdm import tqdm

from pleaserender import util
from pleaserender.core import Plot


class Image(Plot):
    def __init__(
        self,
        system,
        coronagraph,
        observing_scenario,
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
        if not kwargs.get("animation_kwargs"):
            kwargs["animation_kwargs"] = {}
        if "auto_title" not in kwargs["animation_kwargs"]:
            kwargs["animation_kwargs"]["auto_title"] = True

        super().__init__(
            **kwargs,
        )

        self.system = system
        self.coronagraph = coronagraph
        self.observing_scenario = observing_scenario
        self.plot_method = "imshow"
        self.name = "image"

    def generate_data(self, observations):
        obs_ds = Observations.run(self, observations=observations)
        return obs_ds

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

    def draw_plot(self, dataset, animation_value, animation_key, plot_kwargs=None):
        # obs = self.get_observation_object(animation_value, animation_key)
        # necessary_dims = dataset.attrs["dist_attrs_for_images"]
        # sel_call = {}
        # title_id = []
        # if self.imaging_params["plane"] == "coro":
        #     sel_call["xpix(coro)"] = np.arange(obs.coronagraph.npixels)
        #     sel_call["ypix(coro)"] = np.arange(obs.coronagraph.npixels)
        # elif self.imaging_params["plane"] == "det":
        #     sel_call["xpix(det)"] = np.arange(obs.detector_shape[0])
        #     sel_call["ypix(det)"] = np.arange(obs.detector_shape[1])
        # for dim in necessary_dims:
        #     val = getattr(obs, dim)
        #     if isinstance(val, u.Quantity):
        #         val = val.value
        #     elif isinstance(val, Time):
        #         val = val.datetime64
        #     sel_call[dim] = val
        #     title_id.append(val)
        # title_id = tuple(title_id)
        # # f"{attr}={val_str}, "
        # obs_title = dataset.attrs["image_titles"][title_id]
        # photons = dataset.sel(**sel_call)[self.imsel].data
        photons = self.get_frame_data(dataset, animation_value, animation_key)
        if plot_kwargs is None:
            plot_kwargs = self.plot_kwargs
        self.ax.imshow(photons, origin="lower", cmap="viridis")
        # self.ax.set_title(obs_title)
        # photon_counts = dataset[self.coronagraph.name].sel(time=animation_value)

    def get_frame_data(self, dataset, animation_value, animation_key):
        obs = self.get_observation_object(animation_value, animation_key)
        necessary_dims = dataset.attrs["dist_attrs_for_images"]
        sel_call = {}
        # title_id = []
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
            # title_id.append(val)
        # title_id = tuple(title_id)
        # obs_title = dataset.attrs["image_titles"][title_id]
        photons = dataset.sel(**sel_call)[self.imsel].data
        return photons

    def ax_lims_helper(self, data, necessary_axes, equal=False):
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
