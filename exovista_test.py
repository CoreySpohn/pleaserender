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

tp = Time(system_conv.getpattr("T_p")).unix
per = system_conv.getpattr("T").to(u.s).value
e = np.array(system_conv.getpattr("e"))
w = system_conv.getpattr("w").to(u.rad).value
K = system_conv.getpattr("K").to(u.m / u.s).value

rv = calc_RV_from_time(times.unix, tp, per, e, w, K)
rv_df = system_conv.propagate_rv(times)
prop = system_conv.propagate(times, prop="nbody", ref_frame="bary", clean=True)
true_rv = prop.sel(object="star", index=0, ref_frame="bary", prop="nbody")["vz"].values
breakpoint()

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
