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
from pleaserender.exoplanet_plots import Image, Orbit

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
    "wavelength": wavelength,
    "time": times[0],
    "exposure_time": 48 * u.hr,
    "frame_time": 1 * u.hr,
    "include_star": False,
    "include_planets": True,
    "include_disk": False,
    "bandpass": bandpass,
    "spectral_resolution": 100,
    "return_spectrum": True,
    "return_frames": False,
    "return_sources": True,
    "separate_sources": False,
    "wavelength_resolved_flux": True,
    "wavelength_resolved_transmission": True,
    "detector_pixel_scale": 0.001 * u.arcsec / u.pix,
    "detector_shape": (300, 300),
}
observing_scenario = observing_scenario.ObservingScenario(obs_scen)
obs = observation.Observation(coro1, system, observing_scenario)
obs.create_count_rates()
obs.count_photons()


plot_image1 = Image(system, coro1, observing_scenario, ax_kwargs={"title": "Vector"})

breakpoint()
# Create a figure and add the plots
figure_kwargs = {"figsize": (15, 15), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs)
main_figure.please_set_animation_values(times, "time")

# Add plots to the figures
# main_figure.please_add_plot(plot3d)
# main_figure.please_add_plot(plot2d, col=1)
# main_figure.please_add_plot(plot2d_exovista, col=1)
# main_figure.please_add_plot(plot2d_exovista_pix, col=1)
# main_figure.please_add_plot(plot2d_helio, col=1)
# main_figure.please_add_plot(plot2d_sky, col=1)
main_figure.please_add_plot(plot_image1, row=1)
# main_figure.please_add_plot(plot_image2, row=1, col=1)

# Set the animation values and then render
render_settings = {"animation_duration": 5}
# main_figure.please_preview([-1])
main_figure.please_render_video(
    Path("renders/exovista.mp4"), render_settings=render_settings
)
# main_figure.please_render_images(Path("renders/exovista"))