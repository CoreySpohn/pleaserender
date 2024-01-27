import copy

import astropy.units as u
import exoverses.util.misc as misc
import numpy as np
import scipy
import xarray as xr
from astropy.time import Time

import pleaserender.util as util
from pleaserender.core import Scatter


class Orbit(Scatter):
    def __init__(
        self,
        system,
        plane_2d=None,
        orbit_params=None,
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
            plane_2d (str):
                The plane to plot the 2d orbit in ('x', 'y', or 'z')
            star_params (dict):
                Parameters specifying how to plot the star
            planet_params (dict):
                Parameters specifying how to plot the planets
            disk_params(dict):
                Parameters specifying how to plot the disk
            orbit_fit_params(dict):
                Parameters specifying how to plot the orbit fits

        """
        # Set the default orbit params
        # coords can be barycentric or heliocentric
        # integration can be 'kepler' or 'nbody'
        # units can be distance, angular, or pixels
        # if units are angular, then a distance must be specified
        # if units are pixels, then a pixel scale and a distance must be
        # specified
        default_orbit_params = {
            "frame": "barycentric",
            "integration": "kepler",
            "convention": "exovista",
            "unit": u.AU,
            "distance": None,
            "pixel_scale": None,
        }
        self.orbit_params = default_orbit_params
        if orbit_params is not None:
            self.orbit_params.update(orbit_params)

        # Set the default kwargs
        if kwargs.get("axis_keys") is None:
            if plane_2d is None:
                kwargs["axis_keys"] = {"x": "x", "y": "y", "z": "z"}
            else:
                _axis_keys = ["x", "y", "z"]
                _axis_keys.remove(plane_2d)
                kwargs["axis_keys"] = {"x": _axis_keys[0], "y": _axis_keys[1]}

        default_animation_kwargs = {
            "animation_style": "Single point",
            "rotate": {"azim": (25, 65)},
        }
        if kwargs.get("animation_kwargs") is not None:
            default_animation_kwargs.update(
                kwargs.get(
                    "animation_kwargs",
                )
            )
        kwargs["animation_kwargs"] = default_animation_kwargs
        # if kwargs.get("animation_kwargs") is None:
        #     kwargs["animation_kwargs"] = default_animation_kwargs

        default_ax_kwargs = {"equal_lims": True}
        if kwargs.get("ax_kwargs") is not None:
            default_ax_kwargs.update(kwargs.get("ax_kwargs"))
        kwargs["ax_kwargs"] = default_ax_kwargs

        super().__init__(**kwargs)
        # Create a copy of the system so we don't modify the original and we
        # allow for different kwargs for the objects in different plots
        self.system = copy.deepcopy(system)

        # The default is to plot all planets and nothing else
        default_star_params = {"plot": False}
        self.star_params = default_star_params
        if star_params is not None:
            self.star_params.update(star_params)

        # project is used to project of the planet's position and trail
        # onto the x, y, and z axes. The default is to not project anything
        # but the axes can be specified as
        # project = {"point": ["x", "y", "z"], "trail": ["x", "y", "z"]}
        # in the planet_params
        default_project = {}
        default_planet_params = {
            "plot": True,
            "planets_to_plot": np.arange(len(system.planets)),
            "add_trail": True,
            "project": default_project,
            "planet_plot_kwargs": {"color": "w", "edgecolor": "k"},
        }
        self.planet_params = default_planet_params
        if planet_params is not None:
            if planet_params.get("planet_plot_kwargs") is not None:
                self.planet_params["planet_plot_kwargs"].update(
                    planet_params.get("planet_plot_kwargs")
                )
            self.planet_params.update(planet_params)

        default_disk_params = {"plot": False}
        self.disk_params = default_disk_params
        if disk_params is not None:
            self.disk_params.update(disk_params)

        default_orbit_fit_params = {"plot": False}
        self.orbit_fit_params = default_orbit_fit_params
        if orbit_fit_params is not None:
            self.orbit_fit_params.update(orbit_fit_params)

        # Get the max a value of the planets to be plotted which is used
        # for marker size
        self.max_a = max(
            self.system.getpattr("a")[self.planet_params["planets_to_plot"]]
        )
        for planet_ind in self.planet_params["planets_to_plot"]:
            self.add_planet_kwargs(planet_ind, self.planet_params["planet_plot_kwargs"])

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
        # Check that we have the x, y, z coordinates for all objects
        if self.planet_params.get("plot"):
            for planet_ind in self.planet_params["planets_to_plot"]:
                # Check for data of all planets that should be plotted
                planet_data = self.get_planet_da(planet_ind, dataset)
                for var in ["x", "y", "z"]:
                    if var not in planet_data.variables:
                        return False
                    if np.any(np.isnan(planet_data[var])):
                        return False
        dataset = self.convert_units(dataset)
        return True

    def generate_data(self, dataset, animation_values, animation_key):
        """
        Adds the x, y, z coordinates for the planets to the dataset
        """
        times = dataset["time"]
        _da = self.system.propagate(
            Time(times),
            method=self.orbit_params["integration"],
            frame=self.orbit_params["frame"],
            convention=self.orbit_params["convention"],
        )
        dataset = dataset.update(_da)
        dataset = self.convert_units(dataset)

        return dataset

    def draw_plot(self, dataset, animation_value, animation_key):
        n_vals = dataset[animation_key].size
        current_ind = np.argmax(dataset[animation_key].values == animation_value)
        current_view = self.calc_frame_view(current_ind, n_vals)
        if self.planet_params.get("plot"):
            # Start by getting all the planet plotting information
            self.planet_distances = (
                np.zeros(len(self.planet_params["planets_to_plot"])) * u.m
            )
            for i, planet_ind in enumerate(self.planet_params["planets_to_plot"]):
                planet = self.system.planets[planet_ind]
                planet_ds = self.get_planet_da(planet_ind, dataset)

                # Planet's info for marker size
                r_v, r_p, r_pv = util.calc_object_viewer_vectors(
                    planet_ds,
                    current_ind,
                    current_view,
                    self.max_a,
                )
                planet_viewer_angle = util.calc_object_viewer_angle(r_p, r_pv)
                self.planet_marker_size(planet, planet_viewer_angle, factor=0.5)

                # Used for z-ordering
                planet_viewer_orth_dist = util.calc_object_viewer_orth_dist(r_v, r_p)
                self.planet_distances[i] = planet_viewer_orth_dist

            # Set the zorder of the planets
            zorders = 10 - np.argsort(self.planet_distances)
            for i, planet_ind in enumerate(self.planet_params["planets_to_plot"]):
                planet = self.system.planets[planet_ind]
                planet_zorder = zorders[i]
                planet.plot_kwargs["zorder"] = planet_zorder

            # Finally, plot the planets
            for planet_ind in self.planet_params["planets_to_plot"]:
                self.draw_planet(planet_ind, dataset, animation_value, animation_key)

    def get_planet_da(self, planet_ind, dataset):
        planet_dataset = dataset.sel(
            object="planet", index=planet_ind, frame=self.orbit_params["frame"]
        )
        return planet_dataset

    def draw_planet(self, planet_ind, dataset, animation_value, animation_key):
        planet = self.system.planets[planet_ind]
        planet_dataset = self.get_planet_da(planet_ind, dataset)

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

        if self.planet_params.get("project"):
            point_projects = self.planet_params["project"].get("point")
            trail_projects = self.planet_params["project"].get("trail")
            if point_projects is not None:
                for ax_letter in point_projects:
                    self.project_trail(
                        planet_dataset, animation_value, animation_key, ax_letter
                    )
            if trail_projects is not None:
                for ax_letter in trail_projects:
                    self.project_trail(
                        planet_dataset, animation_value, animation_key, ax_letter
                    )

    def convert_units(self, ds):
        ds = misc.add_units(
            ds,
            self.orbit_params["unit"],
            distance=self.orbit_params.get("distance"),
            pixel_scale=self.orbit_params.get("pixel_scale"),
        )

        # Update the axis keys
        for key, val in self.axis_keys.items():
            if val in ["x", "y", "z"] and val is not None:
                new_name = f"{val}({self.orbit_params['unit']})"
                self.axis_keys[key] = new_name
                # Update the axis label
                if self.ax_kwargs[f"{key}label"] == val:
                    self.ax_kwargs[f"{key}label"] = new_name

        return ds

    def add_planet_kwargs(self, planet_ind, planet_plot_kwargs):
        planet = self.system.planets[planet_ind]
        if hasattr(planet, "plot_kwargs"):
            planet_plot_kwargs.update(planet.plot_kwargs)
        planet.plot_kwargs = planet_plot_kwargs
        planet.base_size = self.planet_base_size(planet)

    def planet_marker_size(self, planet, planet_viewer_angle, factor=0.5):
        """
        Make the planet marker smaller when the planet is behind the star in its orbit
        """
        factor_to_use = max(min(10, factor * planet.base_size), 20)
        marker_size = planet.base_size + (
            (180 - planet_viewer_angle.to(u.deg).value) / 180 * factor_to_use
        )
        planet.plot_kwargs["s"] = marker_size

    def planet_base_size(self, planet):
        min_size = 15
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
