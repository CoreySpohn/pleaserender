from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from coronagraphoto import coronagraph
from exoverses.exovista import ExovistaSystem

from pleaserender.core import Figure
from pleaserender.exoplanet_plots import Orbit

# Input files
coronagraph_dir = Path("input/coronagraphs/LUVOIR-B-VC6_timeseries/")
scene = Path("input/scenes/999-HIP_-TYC_SUN-mv_4.83-L_1.00-d_10.00-Teff_5778.00.fits")

# Create a system and coronagraph
coro = coronagraph.Coronagraph(coronagraph_dir)
system = ExovistaSystem(scene)

# Create some simple data for the plots
t_values = Time(np.linspace(2000, 2020, 200), format="decimalyear")

# Create Plot instances
animation_kwargs = {}
# animation_kwargs = {"rotate": {"azim": (1, 89)}, "elev": 0}
# animation_kwargs = {"rotate": {"azim": (44, 46)}}
# animation_kwargs = {"rotate": None, "azim": 270, "elev": 90}
# planet_params = {}
# planet_params_3d = {"planets_to_plot": [1], "project": ["x", "y", "z"]}
planet_params_3d = {"project": {"point": ["x", "y", "z"], "trail": ["x", "y", "z"]}}
# planet_params_2d = {"planets_to_plot": [1]}
planet_params_2d = {}

# system.planets[1].mass = 1 * u.Mjup
# system.star.midplane_I = 0 * u.deg
# system.star.midplanet_PA = 0 * u.deg
# system.planets[0].inc = 90 * u.deg
# system.planets[0].w = 0 * u.deg
# system.planets[0].W = 0 * u.deg
ax_kwargs = {}
plot1 = Orbit(
    system,
    planet_params=planet_params_3d,
    animation_kwargs=animation_kwargs,
    ax_kwargs=ax_kwargs,
)

plot2 = Orbit(
    system, planet_params=planet_params_2d, plane_2d="z", ax_kwargs={"aspect": "equal"}
)

# Create a figure and add the plots
figure_kwargs = {"figsize": (15, 10), "layout": None}
main_figure = Figure(fig_kwargs=figure_kwargs, ncols=2)

# Add plots to the figures
main_figure.please_add_plot(plot1)
main_figure.please_add_plot(plot2, col=1)


# Set the animation values and then render
plt.style.use("dark_background")

main_figure.please_set_animation_values(t_values, "time")
render_settings = {"animation_duration": 5}
# main_figure.please_render_images(Path("renders/test"))
main_figure.please_render_video(
    Path("renders/test.mp4"), render_settings=render_settings
)
