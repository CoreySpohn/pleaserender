import copy

import astropy.units as u
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto import Observation, ObservingScenario
from tqdm import tqdm

from pleaserender.core import Plot


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
        pixels = np.arange(self.coronagraph.npixels)
        if "coronagraph" in dataset.coords:
            # Add the coronagraph name to the dataset's coordinates
            # dataset = dataset.assign_coords(
            #     coronagraph=np.append(
            #         dataset["coronagraph"].values, self.coronagraph.name
            #     )
            # )
            pass
        else:
            new_coords = {
                "x (pix)": pixels,
                "y (pix)": pixels,
                "coronagraph": self.coronagraph.name,
            }
            dataset = dataset.assign_coords(**new_coords)
        image_data = np.zeros(
            (len(times), self.coronagraph.npixels, self.coronagraph.npixels)
        )
        for i, time in enumerate(
            tqdm(times, desc=f"Generating images for {self.coronagraph.name}")
        ):
            self.observing_scenario.scenario["time"] = Time(time)
            self.observation = Observation(
                self.coronagraph,
                self.system,
                self.observing_scenario,
                logging_level="WARNING",
            )
            self.observation.create_count_rates()
            self.observation.count_photons()
            image = self.observation.data["total"].data
            image_data[i] = image
        image_xr = xr.DataArray(
            image_data,
            coords=[times, pixels, pixels],
            dims=["time", "x (pix)", "y (pix)"],
        )
        # image_xr = xr.DataArray(
        #     image_data[np.newaxis, ...],
        #     coords={
        #         "coronagraph": [self.coronagraph.name],
        #         "times": times,
        #         "x (pix)": pixels,
        #         "y (pix)": pixels,
        #     },
        #     dims=["coronagraph", "time", "x (pix)", "y (pix)"],
        # )

        # if "image" in dataset.data_vars:
        #     new_image_data = xr.concat([dataset["image"], image_xr], dim="coronagraph")
        #     dataset["image"] = new_image_data
        #     breakpoint()
        #     # new_coronagraph_coord = np.append(dataset["coronagraph"].values, self.coronagraph.name)
        #     # dataset.assign({ "image": (("coronagraph", "time", "x (pix)", "y (pix)"), new_image_data), "coronagraph": ("coronagraph", new_coronagraph_coord) })
        # else:
        dataset[self.coronagraph.name] = image_xr
        return dataset

    def draw_plot(self, dataset, animation_value, animation_key):
        photon_counts = dataset[self.coronagraph.name].sel(time=animation_value)
        self.ax.imshow(photon_counts, origin="lower", cmap="viridis")

    def ax_lims_helper(self, data, necessary_axes, equal=False):
        for ax_letter in necessary_axes:
            vals = data[self.axis_keys[ax_letter]].data
            self.ax_kwargs["lims"][ax_letter] = (vals.min(), vals.max(), 25, 0)
