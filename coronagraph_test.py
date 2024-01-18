from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto import (coronagraph, observation, observations,
                            observing_scenario, render_engine)
from exoverses.base import Planet, Star, System
from exoverses.exovista import ExovistaSystem
from synphot import Observation, SourceSpectrum, SpectralElement
from synphot.models import (BlackBodyNorm1D, Box1D, Empirical1D, Gaussian1D,
                            GaussianFlux1D)

from pleaserender.core import Figure, Plot, Scatter
from pleaserender.exoplanet_plots import Image, Orbit

# Create some simple data for the plots
times = Time(np.linspace(2000, 2005, 100), format="decimalyear")
# times = Time(np.linspace(2000, 2010, 100), format="decimalyear")

# Input files
coronagraph_dir = Path("input/coronagraphs/LUVOIR-B-VC6_timeseries/")
scene = Path("input/scenes/999-HIP_-TYC_SUN-mv_4.83-L_1.00-d_10.00-Teff_5778.00.fits")
coronagraph_dir = Path("input/coronagraphs/LUVOIR-B-VC6_timeseries/")

# Create a system and coronagraph
coro = coronagraph.Coronagraph(coronagraph_dir)
system = ExovistaSystem(scene)

# GETTING EXOVISTA POSITION DATA
planet_inds = np.arange(len(system.planets))
x_data = np.zeros((len(times), len(planet_inds)))
y_data = np.zeros((len(times), len(planet_inds)))
z_data = np.zeros((len(times), len(planet_inds)))
times64 = times.datetime64
dataset = xr.Dataset(coords={"time": times64, "planet": planet_inds})
for i, planet_ind in enumerate(planet_inds):
    planet = system.planets[planet_ind]
    _x = np.interp(times.decimalyear, planet._t.decimalyear, planet._x.to(u.AU).value)
    _y = np.interp(times.decimalyear, planet._t.decimalyear, planet._y.to(u.AU).value)
    _z = np.interp(times.decimalyear, planet._t.decimalyear, planet._z.to(u.AU).value)
    # _x, _y, _z = planet._x, planet._y, planet._z
    x_data[:, i], y_data[:, i], z_data[:, i] = (_x, _y, _z)
x_xr = xr.DataArray(x_data, coords=[("time", times64), ("planet", planet_inds)])
y_xr = xr.DataArray(y_data, coords=[("time", times64), ("planet", planet_inds)])
z_xr = xr.DataArray(z_data, coords=[("time", times64), ("planet", planet_inds)])

# # Add the exovista data to the dataset
dataset = dataset.assign(_x=x_xr, _y=y_xr, _z=z_xr)


wavelength = 500 * u.nm
frac_bandwidth = 0.15
bandpass = SpectralElement(
    Gaussian1D,
    mean=wavelength,
    stddev=frac_bandwidth * wavelength / np.sqrt(2 * np.pi),
)

obs_scen = {
    "diameter": 8 * u.m,
    "wavelength": wavelength,
    "time": times[0],
    "exposure_time": 48 * u.hr,
    "frame_time": 1 * u.hr,
    "include_star": False,
    "include_planets": True,
    "include_disk": False,
    "bandpass": bandpass,
    "spectral_resolution": 100,
    "return_spectrum": False,
    "return_frames": False,
    "separate_sources": False,
    "wavelength_resolved_flux": False,
    "wavelength_resolved_transmission": False
    # "include_photon_noise": True,
}
observing_scenario = observing_scenario.ObservingScenario(obs_scen)


# Create Plot instances
animation_kwargs = {}
# animation_kwargs = {"rotate": {"azim": (1, 89)}, "elev": 0}
# animation_kwargs = {"rotate": {"azim": (44, 46)}}
# animation_kwargs = {"rotate": None, "azim": 270, "elev": 90}
# planet_params = {}
# planet_params_3d = {"planets_to_plot": [1], "project": ["x", "y", "z"]}
planet_params_3d = {"project": {"point": ["x", "y", "z"], "trail": ["x", "y", "z"]}}
# planet_params_2d = {"planets_to_plot": [1]}
planet_params_2d = {}

# system.planets[1].mass = 1 * u.Mjup
# system.star.midplane_I = 0 * u.deg
# system.star.midplanet_PA = 0 * u.deg
# system.planets[0].inc = 90 * u.deg
# system.planets[0].w = 0 * u.deg
# system.planets[0].W = 0 * u.deg
ax_kwargs = {}
plot3d = Orbit(
    system,
    planet_params=planet_params_3d,
    animation_kwargs=animation_kwargs,
    ax_kwargs=ax_kwargs,
)

plot2d = Orbit(
    system, planet_params=planet_params_2d, plane_2d="z", ax_kwargs={"aspect": "equal"}
)
plot_image = Image(system, coro, observing_scenario)

plot2d_exovista = Orbit(
    system,
    plane_2d="z",
    axis_keys={"x": "_x", "y": "_y"},
    ax_kwargs={"aspect": "equal"},
)
# Create a figure and add the plots
figure_kwargs = {"figsize": (15, 10), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, nrows=2, ncols=2)
main_figure.please_add_dataset(dataset)

# Add plots to the figures
main_figure.please_add_plot(plot3d)
main_figure.please_add_plot(plot_image, col=1)
main_figure.please_add_plot(plot2d, row=1)
main_figure.please_add_plot(plot2d_exovista, row=1, col=1)


# Set the animation values and then render
plt.style.use("dark_background")

main_figure.please_set_animation_values(times, "time")
render_settings = {"animation_duration": 5}
# main_figure.please_render_images(Path("renders/test"))
main_figure.please_render_video(
    Path("renders/test.mp4"), render_settings=render_settings
)
