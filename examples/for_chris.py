import copy
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from coronagraphoto import (
    Coronagraph,
    Observation,
    ObservingScenario,
    Settings,
)
from exoverses.exovista import ExovistaSystem
from lod_unit import lod_eq
from matplotlib.colors import LogNorm, Normalize
from synphot import SpectralElement
from synphot.models import Gaussian1D

from pleaserender.core import Figure
from pleaserender.exoplanet import ObservationFrames, Orbit

plt.style.use("dark_background")

# Create some simple data for the plots
start_time = 2004.35
end_time = 2004.65
ntimes = 200
times = Time(np.linspace(start_time, end_time, ntimes), format="decimalyear")
# Image plot data
obs_times = Time(
    np.linspace(times[ntimes // 2].decimalyear, end_time, 1), format="decimalyear"
)
generation_data = {"start_time": obs_times}

# Input files
coronagraph1_dir = Path("input/coronagraphs/LUVOIR-B-VC6_timeseries/")
scene = Path("input/scenes/999-HIP_-TYC_SUN-mv_4.83-L_1.00-d_10.00-Teff_5778.00.fits")
# scene = Path("input/scenes/999-HIP_-TYC_SUN-mv_3.32-L_1.00-d_5.00-Teff_5778.00.fits")

# Create a system and coronagraph
coro1 = Coronagraph(coronagraph1_dir)
system = ExovistaSystem(scene)
breakpoint()

wavelength = 500 * u.nm
frac_bandwidth = 0.15
bandpass = SpectralElement(
    Gaussian1D,
    mean=wavelength,
    stddev=frac_bandwidth * wavelength / np.sqrt(2 * np.pi),
)

obs_scen = {
    "diameter": 8 * u.m,
    "central_wavelength": wavelength,
    "start_time": obs_times[0],
    "exposure_time": 30 * u.day,
    "frame_time": 0.25 * u.day,
    "bandpass": bandpass,
    "spectral_resolution": 50,
    "detector_pixel_scale": 0.01 * u.arcsec / u.pix,
    "detector_shape": (300, 300),
}
_settings = {
    "include_star": False,
    "include_planets": True,
    "include_disk": True,
    "return_frames": True,
    "return_spectrum": False,
    "return_sources": True,
    "wavelength_resolved_transmission": False,
    "time_invariant_planets": False,
}
observing_scenario = ObservingScenario(custom_scenario=obs_scen)
settings1 = Settings(custom_settings=_settings)

obs1 = Observation(
    coro1, system, observing_scenario, settings1, logging_level="WARNING"
)
settings2 = copy.deepcopy(settings1)
settings2.include_star = True
settings2.include_disk = False
settings2.include_planets = False
obs2 = Observation(
    coro1, system, observing_scenario, settings2, logging_level="WARNING"
)

img_params = {"plane": "coro", "unit": u.arcsec}


def title_gen(plot, plot_initialized):
    if not plot_initialized:
        title = "30 day exposure not started"
    else:
        key = plot.animation_kwargs["title_key"]
        val = plot.state.context[key]
        current_time = Time(val)
        initial_time = Time(plot.data[key].values[0])
        days_passed = (current_time - initial_time).to(u.d)
        if days_passed.value >= 29.5:
            title = "Exposure complete"
        else:
            title = f"Exposure in progress, {days_passed.value:.0f}/30 (day)"
    return title


obsplot = ObservationFrames(
    obs1,
    gen_data=generation_data,
    cumulative=True,
    imaging_params=img_params,
    ax_kwargs={
        "xlabel": "x (arcsec)",
        "ylabel": "",
        "auto_title": True,
        "title_function": title_gen,
    },
)
img_params1 = {"plane": "coro", "unit": u.arcsec, "object": "planet"}
planets = ObservationFrames(
    obs1,
    gen_data=generation_data,
    # cumulative=True,
    imaging_params=img_params1,
    ax_kwargs={
        "title": "Planet signal",
        "ylabel": "y (arcsec)",
        "xlabel": "x (arcsec)",
    },
    plot_kwargs={"norm": LogNorm(vmin=1e-6, vmax=1e7)},
)
img_params2 = {"plane": "coro", "unit": u.arcsec, "object": "disk"}
disk = ObservationFrames(
    obs1,
    gen_data=generation_data,
    imaging_params=img_params2,
    ax_kwargs={
        "title": "Disk frames",
        "ylabel": "",
        "xlabel": "x (arcsec)",
    },
    plot_kwargs={"norm": Normalize(vmin=1e-6, vmax=7e4)},
)
img_params3 = {"plane": "coro", "unit": u.arcsec, "object": "star"}
star = ObservationFrames(
    obs2,
    gen_data=generation_data,
    imaging_params=img_params3,
    ax_kwargs={
        "title": "Speckle pattern",
        "ylabel": "",
        "xlabel": "",
    },
    plot_kwargs={"norm": LogNorm(vmin=1e-6, vmax=1e10)},
)
# frames = SpectralCube(
#     obs1, cumulative=True, gen_data=generation_data, imaging_params=img_params
# )
# spectra = SpectralCube(obs1, gen_data=generation_data, imaging_params=img_params)
# bandpass_plot = Bandpass(bandpass, obs=obs1)
plot3d = Orbit(
    system,
    gen_data={"time": times},
    orbit_params={
        "ref_frame": "helio-sky",
        "propagation": "nbody",
        "unit": u.arcsec,
        "pixel_scale": obs_scen["detector_pixel_scale"],
        "distance": 5 * u.pc,
    },
)
pix_scale = (coro1.pixel_scale * u.pix).to(
    u.arcsec, lod_eq(wavelength, obs_scen["diameter"])
) / u.pix
plot2d = Orbit(
    system,
    gen_data={"time": times},
    plane_2d="z",
    orbit_params={
        "ref_frame": "helio-sky",
        "propagation": "nbody",
        "unit": u.arcsec,
        # "pixel_scale": obs_scen["detector_pixel_scale"],
        "pixel_scale": pix_scale,
        "distance": 5 * u.pc,
    },
    ax_kwargs={
        "lims": {"x": [-0.4, 0.4], "y": [-0.4, 0.4]},
        "aspect": "equal",
        "xlabel": "",
        "title": "Solar system at 5 pc",
    },
)

levels = {1: ["time"]}

# Create a figure and add the plots
render_settings = {"framerate": 25}
# blowup_kwargs = {"figsize": (8, 8), "layout": None}
# blowup_figure = Figure(fig_kwargs=blowup_kwargs, ncols=1, nrows=1)
# blowup_figure.please_set_animation_levels(levels)
# blowup_figure.please_add_plot(obsplot)
# blowup_figure.please_render_video(
#     Path("renders/blowup.mp4"), render_settings=render_settings
# )

figure_kwargs = {"figsize": (16, 8), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, ncols=4, nrows=2)

main_figure.please_set_animation_levels(levels)

main_figure.please_add_plot(obsplot, col=1, colspan=2, rowspan=2)

main_figure.please_add_plot(plot2d)
main_figure.please_add_plot(planets, row=1, shared_plot_data=obsplot)

main_figure.please_add_plot(star, row=0, col=3)
main_figure.please_add_plot(disk, row=1, col=3, shared_plot_data=obsplot)


# Set the animation values and then render
main_figure.please_render_video(
    Path("renders/reduced.mp4"), render_settings=render_settings
)
