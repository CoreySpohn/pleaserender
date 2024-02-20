import copy
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto import (Coronagraph, Observation, Observations,
                            ObservingScenario)
from exoverses.base import Planet, Star, System
from exoverses.exovista import ExovistaSystem
from synphot import SpectralElement
from synphot.models import Gaussian1D
from tqdm import main

from pleaserender.core import Figure, Plot, Scatter
from pleaserender.exoplanet import Bandpass, Image
from pleaserender.exoplanet import Observation as ObsPlot
from pleaserender.exoplanet import ObservationFrames, Orbit, SpectralCube

plt.style.use("dark_background")

# Create some simple data for the plots
times = Time(np.linspace(2000, 2005, 10), format="decimalyear")

# Input files
coronagraph1_dir = Path("input/coronagraphs/LUVOIR-B-VC6_timeseries/")
coronagraph2_dir = Path(
    "input/coronagraphs/LUVOIR-A_APLC_18bw_medFPM_2021-05-07_Dyn10pm-nostaticabb"
)
coronagraph3_dir = Path("input/coronagraphs/usort_offaxis_optimal_order_6/")
scene = Path("input/scenes/999-HIP_-TYC_SUN-mv_4.83-L_1.00-d_10.00-Teff_5778.00.fits")

# Create a system and coronagraph
coro1 = Coronagraph(coronagraph1_dir)
# coro2 = coronagraph.Coronagraph(coronagraph2_dir)
# coro3 = coronagraph.Coronagraph(coronagraph3_dir)
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
    "start_time": times[0],
    "exposure_time": 30 * u.day,
    "frame_time": 10 * u.day,
    "include_star": False,
    "include_planets": True,
    "include_disk": False,
    "bandpass": bandpass,
    "spectral_resolution": 20,
    # "return_frames": True,
    "return_frames": True,
    "return_spectrum": True,
    # "separate_sources": False,
    # "wavelength_resolved_flux": False,
    "wavelength_resolved_transmission": True,
    "time_invariant_planets": False,
    "detector_pixel_scale": 0.01 * u.arcsec / u.pix,
    "detector_shape": (300, 300),
}
observing_scenario = ObservingScenario(obs_scen)
obs1 = Observation(coro1, system, observing_scenario, logging_level="WARNING")


obs_times = Time(np.linspace(2000, 2005, 1), format="decimalyear")
generation_data = {"start_time": obs_times}
obs = ObsPlot(obs1, gen_data=generation_data, imaging_params={"plane": "coro"})
frames = ObservationFrames(
    obs1, gen_data=generation_data, imaging_params={"plane": "coro"}
)
spectra = SpectralCube(obs1, gen_data=generation_data, imaging_params={"plane": "coro"})
bandpass_plot = Bandpass(bandpass, obs=obs1)
plot3d = Orbit(
    system,
    gen_data={"time": times},
    orbit_params={
        "ref_frame": "helio-sky",
        "propagation": "nbody",
        "unit": u.AU,
    },
)
plot2d = Orbit(
    system,
    gen_data={"time": times},
    plane_2d="z",
    orbit_params={
        "ref_frame": "helio-sky",
        "propagation": "nbody",
        "unit": u.AU,
    },
)

# Create a figure and add the plots
figure_kwargs = {"figsize": (10, 10), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, ncols=3, nrows=2)
# animation_info = {"time": {"method": "value", "initial": times[0]}, "frame": None}
levels = {0: ["start_time"], 1: ["time"], 2: ["spectral_wavelength(nm)"]}
main_figure.please_set_animation_levels(levels)
# main_figure.please_set_secondary_animation_keys(["frame"])

main_figure.please_add_plot(obs)
main_figure.please_add_plot(frames, col=1, shared_plot_data=obs)
main_figure.please_add_plot(spectra, col=2, shared_plot_data=obs)
main_figure.please_add_plot(bandpass_plot, row=1, col=2, colspan=1, draw_with=spectra)
# main_figure.please_add_plot(plot3d, row=1, col=0)
# main_figure.please_add_plot(plot2d, row=1, col=1)


# Set the animation values and then render
render_settings = {"animation_duration": 5, "framerate": 25}
main_figure.please_render_video(
    Path("renders/new.mp4"), render_settings=render_settings
)
# main_figure.please_preview([-1])
# main_figure.please_render_images(Path("renders/non_static_exposure"))
# main_figure.please_render_images(Path("renders/exovista"))
