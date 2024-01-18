import copy

import astropy.units as u
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto.observation import Observation
from coronagraphoto.observing_scenario import ObservingScenario

from pleaserender.core import Plot

# from pleaserender.exoplanet_plots import Orbit


class Image(Plot):
    def __init__(
        self,
        system,
        coronagraph,
        observing_scenario,
        **kwargs,
    ):
        super().__init__(
            **kwargs,
        )
        self.axis_keys = {"x": "x (pix)", "y": "y (pix)"}
        self.system = system
        self.coronagraph = coronagraph
        self.observing_scenario = observing_scenario
        self.plot_method = "imshow"

        # self.observation = observation

    def verify_data(self, dataset, animation_values, animation_key):
        return False

    def generate_data(self, dataset, animation_values, animation_key):
        times = dataset["time"]
        # a
        pixels = np.arange(self.coronagraph.npixels)
        new_coords = {
            "x (pix)": pixels,
            "y (pix)": pixels,
            "coronagraph": self.coronagraph.name,
        }
        if "coronagraph" in dataset.coords:
            breakpoint()
        dataset = dataset.assign_coords(**new_coords)
        image_data = np.zeros(
            (len(times), self.coronagraph.npixels, self.coronagraph.npixels)
        )
        for i, time in enumerate(times):
            self.observing_scenario.scenario["time"] = Time(time)
            self.observation = Observation(
                self.coronagraph,
                self.system,
                self.observing_scenario,
            )
            self.observation.create_count_rates()
            self.observation.count_photons()
            image = self.observation.data
            image_data[i] = image["total"].data
        image_xr = xr.DataArray(
            image_data,
            coords=[times, pixels, pixels],
            dims=["time", "x (pix)", "y (pix)"],
        )
        dataset = dataset.assign(image=image_xr)
        return dataset

    def draw_plot(self, dataset, animation_value, animation_key):
        photon_counts = dataset.image.sel(
            time=animation_value, coronagraph=self.coronagraph.name
        )
        self.ax.imshow(photon_counts, origin="lower", cmap="viridis")

    def ax_lims_helper(self, data, necessary_axes, equal=False):
        for ax_letter in necessary_axes:
            vals = data[self.axis_keys[ax_letter]].data
            self.ax_kwargs["lims"][ax_letter] = (vals.min(), vals.max(), 25, 0)
