from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from coronagraphoto import Coronagraph, Observation, ObservingScenario
from exoverses.exovista import ExovistaSystem
from lod_unit import lod_eq
from synphot import SpectralElement
from synphot.models import Gaussian1D

from pleaserender.core import Figure
from pleaserender.exoplanet import ObservationFrames, Orbit

plt.style.use("dark_background")

# Create some simple data for the plots
start_time = 2004.25
end_time = 2004.75
ntimes = 50
times = Time(np.linspace(start_time, end_time, ntimes), format="decimalyear")
# Image plot data
obs_times = Time(
    np.linspace(times[ntimes // 2].decimalyear, end_time, 1), format="decimalyear"
)
generation_data = {"start_time": obs_times}

# Input files
coronagraph1_dir = Path("input/coronagraphs/LUVOIR-B-VC6_timeseries/")
# scene = Path("input/scenes/999-HIP_-TYC_SUN-mv_4.83-L_1.00-d_10.00-Teff_5778.00.fits")
scene = Path("input/scenes/999-HIP_-TYC_SUN-mv_3.32-L_1.00-d_5.00-Teff_5778.00.fits")

# Create a system and coronagraph
coro1 = Coronagraph(coronagraph1_dir)
system = ExovistaSystem(scene)

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
    "frame_time": 3 * u.day,
    "include_star": False,
    "include_planets": True,
    "include_disk": False,
    "bandpass": bandpass,
    "spectral_resolution": 50,
    "return_frames": True,
    "return_spectrum": False,
    "return_sources": True,
    "wavelength_resolved_transmission": False,
    "time_invariant_planets": False,
    "detector_pixel_scale": 0.01 * u.arcsec / u.pix,
    "detector_shape": (300, 300),
}
observing_scenario = ObservingScenario(obs_scen)
obs1 = Observation(coro1, system, observing_scenario, logging_level="WARNING")

img_params = {"plane": "coro", "unit": u.arcsec}


def title_gen(plot, plot_initialized):
    if not plot_initialized:
        title = "Exposure not started"
    else:
        key = plot.animation_kwargs["title_key"]
        val = plot.state.context[key]
        current_time = Time(val)
        initial_time = Time(plot.data[key].values[0])
        days_passed = (current_time - initial_time).to(u.d)
        if days_passed.value == 29.5:
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
        "ylabel": "y (arcsec)",
        "auto_title": True,
        "title_function": title_gen,
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
        "title": "Helio-sky",
    },
)
plot2d_2 = Orbit(
    system,
    gen_data={"time": times},
    plane_2d="z",
    orbit_params={
        "ref_frame": "bary",
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
        "title": "Bary",
    },
)

# Create a figure and add the plots
figure_kwargs = {"figsize": (15, 7.5), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, ncols=3, nrows=2)
# main_figure = Figure(fig_kwargs=figure_kwargs, ncols=2)
# levels = {1: ["time"], 2: ["spectral_wavelength(nm)"]}
levels = {1: ["time"]}

main_figure.please_set_animation_levels(levels)

main_figure.please_add_plot(plot2d)
main_figure.please_add_plot(plot2d_2, row=1)
main_figure.please_add_plot(obsplot, col=1, colspan=2, rowspan=2)


# Set the animation values and then render
render_settings = {"framerate": 20}
main_figure.please_render_video(
    Path("renders/orbit.mp4"), render_settings=render_settings
)
