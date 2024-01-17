import copy

import astropy.units as u
import numpy as np
import scipy
import xarray as xr
from astropy.time import Time

from pleaserender.core import Scatter
from pleaserender.util import calc_observer_position


class Orbit(Scatter):
    def __init__(
        self,
        system,
        star_params=None,
        planet_params=None,
        disk_params=None,
        orbit_fit_params=None,
        **kwargs,
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
        self.unit = u.AU

        # Set the default kwargs
        if kwargs.get("axis_keys") is None:
            kwargs["axis_keys"] = {"x": "x", "y": "y", "z": "z"}

        default_animation_kwargs = {
            "animation_style": "Single point",
            "rotate": {"azim": (15, 75)},
        }
        default_animation_kwargs.update(
            kwargs.get(
                "animation_kwargs",
            )
        )
        kwargs["animation_kwargs"] = default_animation_kwargs
        # if kwargs.get("animation_kwargs") is None:
        #     kwargs["animation_kwargs"] = default_animation_kwargs

        # default_ax_kwargs = {}
        # if kwargs.get("ax_kwargs") is None:
        #     kwargs["ax_kwargs"] = default_ax_kwargs

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
            "add_trail": True,
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

        if self.planet_params.get("plot"):
            # Check that we have the planet coordinate
            if "planet" not in dataset.coords:
                return False
            for planet_ind in self.planet_params["planets_to_plot"]:
                # Check for data of all planets that should be plotted
                animation_data = dataset.sel(planet=planet_ind)[animation_key].values

                # Check if each element in animation_values is in animation_data
                all_values_present = np.isin(animation_values, animation_data).all()
                if not all_values_present:
                    return False
        return all_values_present

    def generate_data(self, dataset, animation_values, animation_key):
        times = dataset["time"]
        if self.planet_params.get("plot"):
            # Get the max a value of the planets to be plotted which is used
            # for marker size
            self.max_a = max(
                self.system.getpattr("a")[self.planet_params["planets_to_plot"]]
            )

            # Add a "planet" coordinate to the dataset
            dataset = dataset.assign_coords(
                planet=self.planet_params["planets_to_plot"]
            )
            planet_inds = self.planet_params["planets_to_plot"]
            x_data = np.zeros((len(times), len(planet_inds)))
            y_data = np.zeros((len(times), len(planet_inds)))
            z_data = np.zeros((len(times), len(planet_inds)))
            for i, planet_ind in enumerate(planet_inds):
                planet = self.system.planets[planet_ind]
                _x, _y, _z = planet.calc_vectors(Time(times))
                x_data[:, i], y_data[:, i], z_data[:, i] = (
                    _x.to(self.unit).value,
                    _y.to(self.unit).value,
                    _z.to(self.unit).value,
                )
                self.add_planet_kwargs(planet_ind)
            x_xr = xr.DataArray(
                x_data, coords=[("time", times.data), ("planet", planet_inds)]
            )
            y_xr = xr.DataArray(
                y_data, coords=[("time", times.data), ("planet", planet_inds)]
            )
            z_xr = xr.DataArray(
                z_data, coords=[("time", times.data), ("planet", planet_inds)]
            )

            # # Add the data to the dataset
            dataset = dataset.assign(x=x_xr, y=y_xr, z=z_xr)

        return dataset

    def draw_plot(self, dataset, animation_value, animation_key):
        n_vals = dataset[animation_key].size
        current_ind = np.argmax(dataset[animation_key].values == animation_value)
        current_view = self.calc_frame_view(dataset, current_ind, n_vals)
        if self.planet_params.get("plot"):
            # Start by getting all the planet plotting information
            self.planet_distances = (
                np.zeros(len(self.planet_params["planets_to_plot"])) * self.unit
            )
            for planet_ind in self.planet_params["planets_to_plot"]:
                planet = self.system.planets[planet_ind]
                planet_dataset = dataset.sel(planet=planet_ind)
                self.planet_marker_size(
                    planet, planet_dataset, current_ind, current_view
                )
            # Set the zorder of the planets
            zorders = 10 - np.argsort(self.planet_distances)
            string = ""
            for i, planet_ind in enumerate(self.planet_params["planets_to_plot"]):
                planet = self.system.planets[planet_ind]
                planet_zorder = zorders[i]
                planet.plot_kwargs["zorder"] = planet_zorder
                string += f"{planet_zorder:.2f}, "
            self.ax.set_title(string)

            # Finally, plot the planets
            for planet_ind in self.planet_params["planets_to_plot"]:
                self.draw_planet(planet_ind, dataset, animation_value, animation_key)

    def draw_planet(self, planet_ind, dataset, animation_value, animation_key):
        planet = self.system.planets[planet_ind]
        planet_dataset = dataset.sel(planet=planet_ind)

        # Plot the planet
        self.generic_plot(
            planet_dataset,
            animation_value,
            animation_key,
            plot_kwargs=planet.plot_kwargs,
        )
        if self.planet_params.get("add_trail"):
            # Add a dashed line behind the planet
            self.generic_plot(
                planet_dataset,
                animation_value,
                animation_key,
                plot_method="plot",
                animation_kwargs={"animation_style": "Cumulative"},
                plot_kwargs={"color": "w", "linestyle": "dashed", "alpha": 0.5},
            )

    def add_planet_kwargs(self, planet_ind, **kwargs):
        planet = self.system.planets[planet_ind]
        default_kwargs = {
            "s": 25,
            "color": "w",
            "edgecolor": "k",
        }
        if not hasattr(planet, "plot_kwargs"):
            planet.plot_kwargs = default_kwargs
        planet.plot_kwargs.update(kwargs)
        planet.base_size = self.planet_base_size(planet)

    def create_axes_config(self, data):
        """
        Method that fills ax_kwargs with defaults for things not already
        specified.
        Currently sets the axis limits.
        """
        # Set the axis limits if they are not specified
        necessary_axes = ["x", "y"]
        if self.axis_keys.get("z") is not None:
            necessary_axes.append("z")

        for axis in necessary_axes:
            max = 0
            # Set the axis limits to be offset from the min and max values
            # of the data if they are not specified
            min_value = data[self.axis_keys.get(axis)].min()
            max_value = data[self.axis_keys.get(axis)].max()
            offset = 0.025 * np.abs(max_value - min_value)
            if np.abs(min_value - offset) > max:
                max = np.abs(min_value - offset)
            if np.abs(max_value + offset) > max:
                max = np.abs(max_value + offset)
        for axis in necessary_axes:
            if self.ax_kwargs.get(f"{axis}lim") is None:
                self.ax_kwargs[f"{axis}lim"] = (-max, max)

    def planet_marker_size(
        self, planet, dataset, current_ind, current_view, factor=0.5
    ):
        """
        Make the planet marker smaller when the planet is behind the star in its orbit
        """
        obs_pos = calc_observer_position(
            current_view["elev"], current_view["azim"], self.max_a
        )
        planet_pos = (
            np.array(
                [
                    dataset["x"][current_ind],
                    dataset["y"][current_ind],
                    dataset["z"][current_ind],
                ]
            )
            * self.unit
        )
        # self.ax.set_title(
        #     f"obs_pos = {obs_pos}, planet_pos = [{dataset['x'][current_ind]:.2f},{dataset['y'][current_ind]:.2f},{dataset['z'][current_ind]:.2f}]"
        # )
        # b = np.linalg.norm(obs_pos)
        r_po = planet_pos - obs_pos
        r_ps = planet_pos
        angle = np.arccos(
            -np.dot(r_po, r_ps) / (np.linalg.norm(r_po) * np.linalg.norm(r_ps))
        )
        factor_to_use = max(min(20, factor * planet.base_size), 10)
        marker_size = planet.base_size + (
            (180 - angle.to(u.deg).value) / 180 * factor_to_use
        )
        planet.plot_kwargs["s"] = marker_size

        # Calculating the orthogonal
        r_os = obs_pos
        r_os_hat = r_os / np.linalg.norm(r_os)
        scalar_proj = np.dot(r_ps, r_os_hat)
        planet_orth_dist = np.linalg.norm(obs_pos) - scalar_proj
        self.planet_distances[int(dataset.planet)] = planet_orth_dist
        # self.ax.set_title(f"{planet.obs_dist:.2f}")
        # projected = np.dot(r_po, r_ps) / np.dot(r_ps, r_ps) * r_ps
        # planet.plot_kwargs["zorder"] = angle.to(u.deg).value

        # return marker_size

    def planet_base_size(self, planet):
        min_size = 20
        max_size = 40
        minv = 1 * u.M_earth
        maxv = 1 * u.M_jupiter
        mass_to_use = max(min(planet.mass, 1 * u.M_jupiter), 1 * u.M_earth)
        sizes = np.logspace(np.log10(min_size), np.log10(max_size), 10)
        masses = np.logspace(0, np.log10(maxv / minv).decompose().value, 10) * minv
        base_s = scipy.interpolate.interp1d(masses, sizes)(mass_to_use.to(u.M_earth))
        return base_s

    def adjust_settings(self, dataset, animation_value, animation_key):
        super().generic_settings(dataset, animation_value, animation_key)
        if self.is_3d:
            self.darken_3d_background()

    def darken_3d_background(self):
        color = (0.0, 0.0, 0.0, 1.0)
        self.ax.xaxis.set_pane_color(color)
        self.ax.yaxis.set_pane_color(color)
        self.ax.zaxis.set_pane_color(color)
