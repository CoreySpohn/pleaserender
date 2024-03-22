from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from coronagraphoto import Coronagraph, Observation, ObservingScenario
from exoverses.exovista import ExovistaSystem
from synphot import SpectralElement
from synphot.models import Gaussian1D

from pleaserender.core import Figure
from pleaserender.exoplanet import Observation as ObsPlot
from pleaserender.exoplanet import ObservationFrames, Orbit

plt.style.use("dark_background")

# Create some simple data for the plots
start_time = 2000
end_time = 2002
# Orbit plot data
times = Time(np.linspace(start_time, 2002.5, 100), format="decimalyear")
# Image plot data
# obs_times = Time(np.linspace(start_time, end_time, 1), format="decimalyear")
obs_times = Time(np.linspace(start_time, end_time, 5), format="decimalyear")
generation_data = {"start_time": obs_times}

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
    "frame_time": 2 * u.day,
    "include_star": False,
    "include_planets": True,
    "include_disk": False,
    "bandpass": bandpass,
    "spectral_resolution": 50,
    # "return_frames": True,
    "return_frames": True,
    "return_spectrum": False,
    # "separate_sources": False,
    # "wavelength_resolved_flux": False,
    "wavelength_resolved_transmission": False,
    "time_invariant_planets": False,
    "detector_pixel_scale": 0.01 * u.arcsec / u.pix,
    "detector_shape": (300, 300),
}
observing_scenario = ObservingScenario(obs_scen)
obs1 = Observation(coro1, system, observing_scenario, logging_level="WARNING")


img_params = {"plane": "coro", "unit": u.arcsec}
obs = ObsPlot(
    obs1,
    gen_data=generation_data,
    imaging_params=img_params,
    ax_kwargs={"title": "Full image"},
)
frames = ObservationFrames(
    obs1,
    gen_data=generation_data,
    cumulative=True,
    imaging_params=img_params,
    ax_kwargs={"title": "Current frame"},
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
        "distance": 10 * u.pc,
    },
)
plot2d = Orbit(
    system,
    gen_data={"time": times},
    plane_2d="z",
    orbit_params={
        "ref_frame": "helio-sky",
        "propagation": "nbody",
        "unit": u.arcsec,
        "pixel_scale": obs_scen["detector_pixel_scale"],
        "distance": 10 * u.pc,
    },
    ax_kwargs={"lims": {"x": [-0.4, 0.4], "y": [-0.4, 0.4]}, "aspect": "equal"},
)

# Create a figure and add the plots
figure_kwargs = {"figsize": (8, 8), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, ncols=2, nrows=2)
# main_figure = Figure(fig_kwargs=figure_kwargs, ncols=2)
# levels = {1: ["time"], 2: ["spectral_wavelength(nm)"]}
levels = {1: ["time"]}

main_figure.please_set_animation_levels(levels)

main_figure.please_add_plot(obs)
main_figure.please_add_plot(frames, col=1, shared_plot_data=obs)
# main_figure.please_add_plot(spectra, col=2, shared_plot_data=obs)
# main_figure.please_add_plot(bandpass_plot, row=1, col=2, draw_with=spectra)
main_figure.please_add_plot(plot3d, row=1, col=0)
main_figure.please_add_plot(plot2d, row=1, col=1, shared_plot_data=plot3d)


# Set the animation values and then render
render_settings = {"animation_duration": 5, "framerate": 25}
main_figure.please_render_video(
    Path("renders/new.mp4"), render_settings=render_settings
)
