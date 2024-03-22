from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from coronagraphoto import coronagraph
from exoverses.exovista import ExovistaSystem

from pleaserender.core import Figure
from pleaserender.exoplanet import Orbit

plt.style.use("dark_background")

# Create some simple data for the plots
times = Time(np.linspace(2000, 2002, 100), format="decimalyear")

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
system_conv = ExovistaSystem(scene, convert=True)

# conv = system_conv.propagate(times)
plot2d = Orbit(
    system,
    gen_data={"time": times},
    plane_2d="z",
    orbit_params={
        "ref_frame": "bary-sky",
        "propagation": "nbody",
        "unit": u.AU,
    },
    ax_kwargs={"title": "ExoVista nbody"},
)
plot2d_conv = Orbit(
    system_conv,
    gen_data={"time": times},
    plane_2d="z",
    orbit_params={
        "ref_frame": "bary",
        "propagation": "kepler",
        "unit": u.AU,
    },
    ax_kwargs={"title": "EXOSIMS Keplerian"},
)
figure_kwargs = {"figsize": (10, 5), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, ncols=2)
# main_figure.please_add_dataset(dataset)
levels = {0: ["time"]}
main_figure.please_set_animation_levels(levels)

# Add plots to the figures
main_figure.please_add_plot(plot2d, col=0)
main_figure.please_add_plot(plot2d_conv, col=1)

# Set the animation values and then render
render_settings = {"frame_rate": 30}
# main_figure.please_preview([-1])
main_figure.please_render_video(
    Path("renders/exovista.mp4"), render_settings=render_settings
)

# system.planets = system.planets[:3]

# GETTING EXOVISTA POSITION DATA
# dataset = system.create_dataset(times)
# planet_inds = np.arange(len(system.planets))
# times64 = times.datetime64
# x_data = np.zeros((len(times), len(planet_inds)))
# y_data = np.zeros((len(times), len(planet_inds)))
# z_data = np.zeros((len(times), len(planet_inds)))
# x_pix_data = np.zeros((len(times), len(planet_inds)))
# y_pix_data = np.zeros((len(times), len(planet_inds)))
# # dataset = xr.Dataset(coords={"time": times64, "index": planet_inds, "body": "planet"})
# for i, planet_ind in enumerate(planet_inds):
#     planet = system.planets[planet_ind]
#     _x = np.interp(times.decimalyear, planet._t.decimalyear, planet._x.to(u.AU).value)
#     _y = np.interp(times.decimalyear, planet._t.decimalyear, planet._y.to(u.AU).value)
#     _z = np.interp(times.decimalyear, planet._t.decimalyear, planet._z.to(u.AU).value)
#     _x_pix = np.interp(times.decimalyear, planet._t.decimalyear, planet._x_pix.value)
#     _y_pix = np.interp(times.decimalyear, planet._t.decimalyear, planet._y_pix.value)
#     # _x, _y, _z = planet._x, planet._y, planet._z
#     x_data[:, i], y_data[:, i], z_data[:, i] = (_x, _y, _z)
#     x_pix_data[:, i], y_pix_data[:, i] = _x_pix, _y_pix
# coords = [
#     ("time", times64),
#     ("index", planet_inds),
#     ("object", ["planet"]),
#     ("frame", ["helio-sky"]),
# ]
# x_xr = xr.DataArray(x_data[..., np.newaxis, np.newaxis], coords=coords)
# y_xr = xr.DataArray(y_data[..., np.newaxis, np.newaxis], coords=coords)
# z_xr = xr.DataArray(z_data[..., np.newaxis, np.newaxis], coords=coords)
# x_pix_xr = xr.DataArray(x_pix_data[..., np.newaxis, np.newaxis], coords=coords)
# y_pix_xr = xr.DataArray(y_pix_data[..., np.newaxis, np.newaxis], coords=coords)
# # # Add the exovista data to the dataset
# dataset = dataset.assign(
#     ev_x=x_xr, ev_y=y_xr, ev_z=z_xr, ev_x_pix=x_pix_xr, ev_y_pix=y_pix_xr
# )
# dataset["ev_x_pix"].sel(object="planet", index=0, frame="helio-sky")


# wavelength = 500 * u.nm
# frac_bandwidth = 0.15
# bandpass = SpectralElement(
#     Gaussian1D,
#     mean=wavelength,
#     stddev=frac_bandwidth * wavelength / np.sqrt(2 * np.pi),
# )

# obs_scen = {
#     "diameter": 8 * u.m,
#     "central_wavelength": wavelength,
#     "time": times[0],
#     "exposure_time": 48 * u.hr,
#     "frame_time": 1 * u.hr,
#     "include_star": False,
#     "include_planets": True,
#     "include_disk": False,
#     "bandpass": bandpass,
#     "spectral_resolution": 100,
#     "return_spectrum": False,
#     "return_frames": False,
#     "return_sources": True,
#     "wavelength_resolved_flux": False,
#     "wavelength_resolved_transmission": False,
#     "detector_pixel_scale": 0.005 * u.arcsec / u.pix,
#     "detector_shape": (300, 300),
# }
# observing_scenario = observing_scenario.ObservingScenario(obs_scen)
# obs1 = observation.Observation(
#     coro1, system, observing_scenario, logging_level="WARNING"
# )
# obs2 = observation.Observation(
#     coro2, system, observing_scenario, logging_level="WARNING"
# )


# # Define plots
# # planet_params_3d = {"project": {"point": ["x", "y", "z"], "trail": ["x", "y", "z"]}}
# # plot3d = Orbit(system, planet_params=planet_params_3d)
# # plot_image1 = Image(system, coro1, observing_scenario, ax_kwargs={"title": "Vector"})
# # plot_image2 = Image(system, coro2, observing_scenario, ax_kwargs={"title": "APLC"})

# # planet_params = {"planets_to_plot": [0], "planet_plot_kwargs": {"color": "blue"}}
# planet_params = {"project": {"point": ["z"], "trail": ["z"]}}
# plot3d = Orbit(
#     system,
#     # plane_2d="z",
#     # ax_kwargs={"aspect": "equal", "lims": {"x": (-10, 10), "y": (-10, 10)}},
#     orbit_params={
#         "ref_frame": "helio-sky",
#         "propagation": "nbody",
#         "unit": u.AU,
#     },
#     planet_params=planet_params,
# )
# plot2d = Orbit(
#     system,
#     plane_2d="z",
#     orbit_params={
#         "ref_frame": "helio-sky",
#         "propagation": "nbody",
#         "unit": u.AU,
#     },
# )
# # plot2d_helio = Orbit(
# #     system,
# #     plane_2d="z",
# #     ax_kwargs={"aspect": "equal"},
# #     orbit_params={
# #         "frame": "helio-sky",
# #         "propagation": "nbody",
# #         "unit": u.AU,
# #     },
# # )
# # plot2d_sky = Orbit(
# #     system,
# #     plane_2d="z",
# #     ax_kwargs={"aspect": "equal"},
# #     orbit_params={"frame": "sky", "unit": u.arcsec, "distance": 10 * u.pc},
# # )
# # plot2d_exovista = Orbit(
# #     system,
# #     plane_2d="z",
# #     axis_keys={"x": "ev_x", "y": "ev_y"},
# #     ax_kwargs={"aspect": "equal", "lims": {"x": (-10, 10), "y": (-10, 10)}},
# #     orbit_params={"frame": "helio-sky", "propagation": "nbody", "unit": u.AU},
# # )
# # plot_image1 = Image(system, coro1, observing_scenario, ax_kwargs={"title": "Vector"})
# plot_image1 = Image(
#     system,
#     coro2,
#     observing_scenario,
#     ax_kwargs={"title": "APLC (coro)"},
#     imaging_params={"plane": "coro"},
# )
# plot_image2 = Image(
#     system,
#     coro2,
#     observing_scenario,
#     ax_kwargs={"title": "APLC (det)"},
#     imaging_params={"plane": "det", "object": "planet"},
# )
# # plot2d_exovista_pix = Orbit(
# #     system,
# #     plane_2d="z",
# #     axis_keys={"x": "ev_x_pix", "y": "ev_y_pix"},
# #     ax_kwargs={"aspect": "equal", "lims": {"x": (-1000, 1000), "y": (-1000, 1000)}},
# #     orbit_params={"frame": "helio-sky", "propagation": "nbody"},
# # )

# # Create a figure and add the plots
# figure_kwargs = {"figsize": (10, 10), "layout": None}
# main_figure = Figure(fig_kwargs=figure_kwargs)  # , ncols=2, nrows=2)
# # main_figure.please_add_dataset(dataset)
# main_figure.please_set_animation_values(times, "time")

# # Add plots to the figures
# main_figure.please_add_plot(plot3d)
# # main_figure.please_add_plot(plot2d, col=1)
# # main_figure.please_add_plot(plot2d_exovista, col=1)
# # main_figure.please_add_plot(plot2d_exovista_pix, col=1)
# # main_figure.please_add_plot(plot2d_helio, col=1)
# # main_figure.please_add_plot(plot2d_sky, col=1)
# # main_figure.please_add_plot(plot_image1, row=1)
# # main_figure.please_add_plot(plot_image2, col=1, row=1)

# # Set the animation values and then render
# render_settings = {"animation_duration": 10}
# # main_figure.please_preview([-1])
# main_figure.please_render_video(
#     Path("renders/exovista.mp4"), render_settings=render_settings
# )
# # main_figure.please_render_images(Path("renders/exovista"))
