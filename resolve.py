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
from tqdm import main

from pleaserender.core import Figure, Plot, Scatter
from pleaserender.exoplanet import Image, ObservationFrames, Orbit

plt.style.use("dark_background")

# Create some simple data for the plots
times = Time(np.linspace(2000, 2005, 10), format="decimalyear")

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
    "time": times[0],
    "exposure_time": 30 * 24 * u.hr,
    "frame_time": 24 * u.hr,
    "include_star": False,
    "include_planets": True,
    "include_disk": False,
    # "bandpass": bandpass,
    "bandpass_model": Gaussian1D,
    "frac_bandwidth": frac_bandwidth,
    "spectral_resolution": 100,
    "return_spectrum": False,
    "return_frames": True,
    "return_sources": False,
    "separate_sources": False,
    "wavelength_resolved_flux": False,
    "wavelength_resolved_transmission": False,
    "prop_during_exposure": True,
    "detector_pixel_scale": 0.01 * u.arcsec / u.pix,
    "detector_shape": (300, 300),
}
observing_scenario = observing_scenario.ObservingScenario(obs_scen)
obs = observation.Observation(
    coro1, system, observing_scenario, logging_level="WARNING"
)
# obs.create_count_rates()
# obs.count_photons()


plot_image = ObservationFrames(
    obs, ax_kwargs={"title": "Current frame"}, imaging_params={"plane": "coro"}
)
plot_image_cumu = ObservationFrames(
    obs,
    cumulative=True,
    ax_kwargs={"title": "Cumulative"},
    imaging_params={"plane": "coro"},
)
nframes, *_ = obs.calc_frame_info()
frame_inds = np.arange(nframes)

# Create a figure and add the plots
wavelengths = np.linspace(400, 1000, 2) * u.nm
figure_kwargs = {"figsize": (10, 5), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, ncols=2)
main_figure.please_set_animation_values(frame_inds, "frame")

# Add plots to the figures
# main_figure.please_add_plot(plot3d)
# main_figure.please_add_plot(plot2d, col=1)
# main_figure.please_add_plot(plot2d_exovista, col=1)
# main_figure.please_add_plot(plot2d_exovista_pix, col=1)
# main_figure.please_add_plot(plot2d_helio, col=1)
# main_figure.please_add_plot(plot2d_sky, col=1)
main_figure.please_add_plot(plot_image)
main_figure.please_add_plot(plot_image_cumu, col=1)
# main_figure.please_add_plot(plot_image2, row=1, col=1)

# Set the animation values and then render
render_settings = {"animation_duration": 5}
# main_figure.please_preview([-1])
main_figure.please_render_video(
    Path("renders/non_static_exposure.mp4"), render_settings=render_settings
)
# main_figure.please_render_images(Path("renders/exovista"))
