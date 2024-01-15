import copy

import numpy as np
from pleaseRender.core import Scatter


class Orbit(Scatter):
    def __init__(
        self,
        system,
        star_params=None,
        planet_params=None,
        disk_params=None,
        orbit_fit_params=None,
        **kwargs
    ):
        """
        A plot of the orbits of some set of objects in a planetary system

        Args:
            system (exoverses System):
                An exoverses system to be plotted
            star_params (dict):
                Parameters specifying how to plot the star
            planet_params (dict):
                Parameters specifying how to plot the planets
            disk_params(dict):
                Parameters specifying how to plot the disk
            orbit_fit_params(dict):
                Parameters specifying how to plot the orbit fits

        """
        super().__init__(**kwargs)
        # Create a copy of the system so we don't modify the original and we
        # allow for different kwargs for the objects in different plots
        self.system = copy.deepcopy(system)

        # The default is to plot all planets and nothing else
        default_star_params = {"plot": False}
        self.star_params = default_star_params
        if star_params is not None:
            self.star_params.update(star_params)

        default_planet_params = {
            "plot": True,
            "planets_to_plot": np.arange(len(system.planets)),
        }
        self.planet_params = default_planet_params
        if planet_params is not None:
            self.planet_params.update(planet_params)

        default_disk_params = {"plot": False}
        self.disk_params = default_disk_params
        if disk_params is not None:
            self.disk_params.update(disk_params)

        default_orbit_fit_params = {"plot": False}
        self.orbit_fit_params = default_orbit_fit_params
        if orbit_fit_params is not None:
            self.orbit_fit_params.update(orbit_fit_params)

    def verify_data(self, dataset, animation_values, animation_key):
        """
        Make sure we have x, y, z information for all objects
        Args:
            dataset (xr.Dataset):
                Dataset containing the data and coordinates
            animation_values (numpy.ndarray):
                List of values for the animation
            animation_key (str):
                Key for the animation values in the plot object's dataframe
        """
        # TODO Check for data of all objects that should be plotted
        animation_data = dataset[animation_key].values

        # Check if each element in animation_values is in animation_data
        all_values_present = np.isin(animation_values, animation_data).all()
        return all_values_present

    def draw_plot(self, dataset, animation_value, animation_key):
        if self.planet_params.get("plot"):
            for planet_ind in self.planet_params.planets_to_plot:
                self.draw_planet(planet_ind, dataset, animation_value, animation_key)

    def draw_planet(self, dataset, planet_ind, animation_value, animation_key):
        planet = self.system[planet_ind]

        planet_dataset = dataset.sel(planet=planet_ind)
        super().draw_plot(
            planet_dataset,
            animation_value,
            animation_key,
            plot_kwargs=planet.plot_kwargs,
        )

    def add_planet_kwargs(self, planet_ind, **kwargs):
        planet = self.system[planet_ind]
        default_kwargs = {"label": planet.name}
        if not hasattr(planet, "plot_kwargs"):
            planet.plot_kwargs = {}
        planet.plot_kwargs.update(kwargs)
