import copy

import astropy.units as u
import exoverses.util.misc as misc
import numpy as np
import scipy

import pleaserender.util as util
from pleaserender.core import Scatter


class Orbit(Scatter):
    def __init__(
        self,
        system,
        gen_data=None,
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
        # propagation can be 'kepler' or 'nbody'
        # units can be distance, angular, or pixels
        # if units are angular, then a distance must be specified
        # if units are pixels, then a pixel scale and a distance must be
        # specified
        default_orbit_params = {
            "ref_frame": "bary",
            "propagation": "kepler",
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
            "title_key": "time",
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

        if (
            "auto_title" not in kwargs["ax_kwargs"]
            and "title" not in kwargs["ax_kwargs"]
        ):
            kwargs["ax_kwargs"]["auto_title"] = True

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

        self.gen_data = gen_data

        self.base_sel = {
            "ref_frame": self.orbit_params["ref_frame"],
            "prop": self.orbit_params["propagation"],
        }
        self.animation_key = "time"

    def get_required_keys(self):
        self.required_keys = ["time"]

    def draw_plot(self):
        animation_key = "time"
        animation_value = self.state.context["time"]
        n_vals = self.data[animation_key].size
        current_ind = np.argmax(self.data[animation_key].values == animation_value)
        current_view = self.calc_frame_view(current_ind, n_vals)
        if self.planet_params.get("plot"):
            # Start by getting all the planet plotting information
            self.planet_distances = (
                np.zeros(len(self.planet_params["planets_to_plot"])) * u.m
            )
            for i, planet_ind in enumerate(self.planet_params["planets_to_plot"]):
                planet = self.system.planets[planet_ind]
                planet_ds = self.get_planet_da(planet_ind)

                # Planet's info for marker size
                r_v, r_p, r_pv = util.calc_object_viewer_vectors(
                    planet_ds,
                    current_ind,
                    current_view,
                    self.max_a,
                )
                planet_viewer_angle = util.calc_object_viewer_angle(r_p, r_pv)
                self.planet_marker_size(planet_ind, planet_viewer_angle, factor=1)

                # Used for z-ordering
                planet_viewer_orth_dist = util.calc_object_viewer_orth_dist(r_v, r_p)
                self.planet_distances[i] = planet_viewer_orth_dist

            # Set the zorder of the planets
            zorders = 10 - np.argsort(self.planet_distances)
            for i, planet_ind in enumerate(self.planet_params["planets_to_plot"]):
                planet = self.system.planets[planet_ind]
                planet_zorder = zorders[i]
                planet.plot_kwargs["zorder"] = planet_zorder
                self.system.planets[planet_ind] = planet

            # Finally, plot the planets
            for planet_ind in self.planet_params["planets_to_plot"]:
                self.draw_planet(planet_ind)

    def get_planet_da(self, planet_ind):
        # dist_attrs = self.get_dist_attrs(self.system, lists=False)
        plan_sel = self.base_sel.copy()
        plan_sel.update({"object": "planet", "index": planet_ind})
        planet_dataset = self.data.sel(**plan_sel)
        return planet_dataset

    def get_system_da(self):
        system_dataset = self.data.sel(object="planet", **self.base_sel)
        return system_dataset

    def draw_planet(self, planet_ind):
        planet = self.system.planets[planet_ind]
        plan_sel = {"object": "planet", "index": planet_ind}

        # Plot the planet
        planet_context = self.state.context_sel(extra_sel=plan_sel)
        planet_dataset = self.data.sel(**planet_context)
        self.generic_plot(data=planet_dataset, plot_kwargs=planet.plot_kwargs)

        # Add trail
        if self.planet_params.get("add_trail"):
            # Add a dashed line behind the planet
            trail_sel = self.state.context_sel(
                extra_sel=plan_sel, extra_strategies={"time": "Cumulative"}
            )
            trail_dataset = self.data.sel(**trail_sel)
            trail_kwargs = {"color": "w", "linestyle": "dashed", "alpha": 0.5}
            self.generic_plot(
                data=trail_dataset, plot_method="plot", plot_kwargs=trail_kwargs
            )

        # Project onto the 3d axes
        if self.planet_params.get("project"):
            # Project planet positions as a trail
            trail_projects = self.planet_params["project"].get("trail")
            if trail_projects is not None:
                for ax_letter in trail_projects:
                    self.project("trail", ax_letter, extral_sel=plan_sel)

            # Project planet as a point
            point_projects = self.planet_params["project"].get("point")
            if point_projects is not None:
                for ax_letter in point_projects:
                    self.project("point", ax_letter, extral_sel=plan_sel)

    def project(self, type, ax_letter, extra_sel={}, plot_kwargs=None):
        min_value = self.ax_kwargs["lims"][ax_letter][0]
        default_plot_kwargs = {"zs": min_value, "zdir": ax_letter, "alpha": 0.1}
        if plot_kwargs is not None:
            default_plot_kwargs.update(plot_kwargs)

        plot_kwargs = copy.deepcopy(default_plot_kwargs)
        if type == "trail":
            plot_kwargs["linestyle"] = "-"
            plot_kwargs["color"] = "w"
            strategy = {"time": "Cumulative"}
            method = "plot"
        elif type == "point":
            plot_kwargs["facecolor"] = "grey"
            plot_kwargs["edgecolor"] = "grey"
            strategy = {"time": "Value"}
            method = "scatter"
        sel = self.state.context_sel(extra_sel=extra_sel, extra_strategies=strategy)
        dataset = self.data.sel(**sel)
        self.generic_plot(data=dataset, plot_method=method, plot_kwargs=plot_kwargs)

    def convert_units(self):
        self.data = misc.add_units(
            self.data,
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

    def add_planet_kwargs(self, planet_ind, planet_plot_kwargs):
        planet = self.system.planets[planet_ind]
        if hasattr(planet, "plot_kwargs"):
            planet_plot_kwargs.update(planet.plot_kwargs)
        planet.plot_kwargs = copy.copy(planet_plot_kwargs)
        planet.base_size = self.planet_base_size(planet)
        self.system.planets[planet_ind] = planet

    def planet_marker_size(self, planet_ind, planet_viewer_angle, factor=0.5):
        """
        Make the planet marker smaller when the planet is behind the star in its orbit
        """
        planet = self.system.planets[planet_ind]
        factor_to_use = max(min(35, factor * planet.base_size), 15)
        marker_size = planet.base_size + (
            (180 - planet_viewer_angle.to(u.deg).value) / 180 * factor_to_use
        )
        planet.plot_kwargs["s"] = marker_size
        self.system.planets[planet_ind] = planet

    def planet_base_size(self, planet):
        min_marker_size = 10
        max_marker_size = 40
        # minv = 1 * u.M_earth
        # Mercury mass
        min_mass = (3.3011e23 * u.kg).to(u.M_earth).value
        max_mass = (1 * u.M_jupiter).to(u.M_earth).value
        planet_mass = planet.mass.to(u.M_earth).value
        mass_to_use = max(min(planet_mass, max_mass), min_mass)
        sizes = np.logspace(np.log10(min_marker_size), np.log10(max_marker_size), 10)
        masses = np.logspace(0, np.log10(max_mass / min_mass), 10) * min_mass
        base_s = scipy.interpolate.interp1d(masses, sizes)(mass_to_use)
        return base_s

    def adjust_settings(self):
        super().generic_settings()
        if self.is_3d:
            self.darken_3d_background()

    def darken_3d_background(self):
        color = (0.0, 0.0, 0.0, 1.0)
        self.ax.xaxis.set_pane_color(color)
        self.ax.yaxis.set_pane_color(color)
        self.ax.zaxis.set_pane_color(color)

    def create_axes_config(self):
        # Set the axis limits if they are not specified
        data = self.get_system_da()
        self.handle_axes_limits_and_ticks(
            data=data, equal=self.ax_kwargs.get("equal_lims")
        )

    def check_valid_context(self, fig_context):
        valid = fig_context["time"] in self.data["time"]
        return valid

    def get_dist_attrs(self, system, lists=True):
        if lists:
            dist_attrs = {"name": [system.star.name], "origin": [system.origin]}
        else:
            dist_attrs = {"name": system.star.name, "origin": system.origin}
        return dist_attrs

    def generate_data(self):
        times = self.gen_data["time"]
        self.data = self.system.propagate(
            times,
            prop=self.orbit_params["propagation"],
            ref_frame=self.orbit_params["ref_frame"],
            convention=self.orbit_params["convention"],
            clean=True,
        )

    def process_data(self):
        self.convert_units()
