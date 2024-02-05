import astropy.units as u
import numpy as np
from astropy.time import Time

from pleaserender.exoplanet import Image


class ObservationFrames(Image):
    def __init__(self, observation, cumulative=False, *args, **kwargs):
        super().__init__(observation, *args, **kwargs)
        self.observation = observation
        self.cumulative = cumulative

    def add_extra_sel_call(self, sel_call, obs, animation_value, animation_key):
        if self.cumulative:
            sel_call["frame"] = np.arange(animation_value)
        else:
            sel_call["frame"] = animation_value

        return sel_call

    def process_photons(self, photons):
        if self.cumulative:
            photons = photons.sum("frame")
        photons = photons[self.imsel].data
        return photons
