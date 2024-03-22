from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto import coronagraph, observing_scenario
from exoverses.exovista import ExovistaSystem
from synphot import SpectralElement
from synphot.models import Gaussian1D

from pleaserender.core import Figure
from pleaserender.exoplanet import Orbit

plt.style.use("dark_background")

# Create some simple data for the plots
times = Time(np.linspace(2000, 2005, 100), format="decimalyear")

# Input files
coronagraph1_dir = Path("input/coronagraphs/LUVOIR-B-VC6_timeseries/")
coronagraph2_dir = Path(
    "input/coronagraphs/LUVOIR-A_APLC_18bw_medFPM_2021-05-07_Dyn10pm-nostaticabb"
)
scene = Path("input/scenes/999-HIP_-TYC_SUN-mv_4.83-L_1.00-d_10.00-Teff_5778.00.fits")

# Create a system and coronagraph
coro1 = coronagraph.Coronagraph(coronagraph1_dir)
coro2 = coronagraph.Coronagraph(coronagraph2_dir)
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
    "wavelength_resolved_transmission": False,
    # "include_photon_noise": True,
}
observing_scenario = observing_scenario.ObservingScenario(obs_scen)

1
# Define plots
planet_params_3d = {"project": {"point": ["x", "y", "z"], "trail": ["x", "y", "z"]}}
plot3d = Orbit(system, planet_params=planet_params_3d)
plot2d = Orbit(system, plane_2d="z", ax_kwargs={"aspect": "equal"})
# plot_image1 = Image(system, coro1, observing_scenario, ax_kwargs={"title": "Vector"})
# plot_image2 = Image(system, coro2, observing_scenario, ax_kwargs={"title": "APLC"})
# plot2d_exovista = Orbit(
#     system,
#     plane_2d="z",
#     axis_keys={"x": "_x", "y": "_y"},
#     ax_kwargs={"aspect": "equal"},
# )

# Create a figure and add the plots
figure_kwargs = {"figsize": (15, 10), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, nrows=2, ncols=2)
main_figure.please_add_dataset(dataset)
main_figure.please_set_animation_values(times, "time")

# Add plots to the figures
main_figure.please_add_plot(plot3d)
main_figure.please_add_plot(plot2d, col=1)
# main_figure.please_add_plot(plot_image1, row=1)
# main_figure.please_add_plot(plot_image2, row=1, col=1)

# Set the animation values and then render
render_settings = {"animation_duration": 5}
main_figure.please_render_video(
    Path("renders/coro_test.mp4"), render_settings=render_settings
)
