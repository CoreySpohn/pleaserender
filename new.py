import copy
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from astropy.time import Time
from coronagraphoto import (Coronagraph, Observation, Observations,
                            ObservingScenario, render_engine)
from exoverses.base import Planet, Star, System
from exoverses.exovista import ExovistaSystem
from synphot import SourceSpectrum, SpectralElement
from synphot.models import (BlackBodyNorm1D, Box1D, Empirical1D, Gaussian1D,
                            GaussianFlux1D)
from tqdm import main

from pleaserender.core import Figure, Plot, Scatter
from pleaserender.exoplanet import Image, ObservationFrames, Orbit

plt.style.use("dark_background")

# Create some simple data for the plots
times = Time(np.linspace(2000, 2005, 3), format="decimalyear")

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
    "time": times[0],
    "exposure_time": 30 * u.day,
    "frame_time": 5 * u.day,
    "include_star": False,
    "include_planets": True,
    "include_disk": False,
    "bandpass": bandpass,
    "spectral_resolution": 100,
    "return_frames": True,
    # "separate_sources": False,
    # "wavelength_resolved_flux": False,
    # "wavelength_resolved_transmission": False,
    # "time_invariant_planets": False,
    "detector_pixel_scale": 0.01 * u.arcsec / u.pix,
    "detector_shape": (300, 300),
}
observing_scenario = ObservingScenario(obs_scen)
obs1 = Observation(coro1, system, observing_scenario, logging_level="WARNING")
# obs_scen2 = copy.copy(obs_scen)
# obs_scen2["time_invariant_planets"] = True
# observing_scenario2 = ObservingScenario(obs_scen2)


plot_image1 = ObservationFrames(
    obs1, ax_kwargs={"title": "Current frame"}, imaging_params={"plane": "coro"}
)
# plot_image_cumu1 = ObservationFrames(
#     obs1,
#     cumulative=True,
#     ax_kwargs={"title": "Cumulative"},
#     imaging_params={"plane": "coro"},
# )
plot2d = Orbit(
    system,
    plane_2d="z",
    orbit_params={
        "ref_frame": "helio-sky",
        "propagation": "nbody",
        "unit": u.AU,
    },
)

nframes, *_ = obs1.calc_frame_info()
frame_inds = np.arange(nframes)

# Create a figure and add the plots
figure_kwargs = {"figsize": (10, 10), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs)
generation_data = {"time": times}
animation_info = {"time": {"method": "value", "initial": times[0]}, "frame": None}
main_figure.please_set_animation_info(generation_data, animation_info)
# main_figure.please_set_secondary_animation_keys(["frame"])

main_figure.please_add_plot(plot_image1)
# main_figure.please_add_plot(
#     plot_image_cumu1, col=1, row=1, shared_plot_data=plot_image1
# )


# Set the animation values and then render
render_settings = {"animation_duration": 5}
main_figure.please_render_video(
    Path("renders/new.mp4"), render_settings=render_settings
)
# main_figure.please_preview([-1])
# main_figure.please_render_images(Path("renders/non_static_exposure"))
# main_figure.please_render_images(Path("renders/exovista"))
