import astropy.units as u
import numpy as np
from astropy.time import Time

from pleaserender.exoplanet import Image


class ObservationFrames(Image):
    def __init__(self, observation, cumulative=False, *args, **kwargs):
        super().__init__(observation, *args, **kwargs)
        self.observation = observation
        self.cumulative = cumulative

    def add_extra_sel_call(self, sel_call, obs, draw_data):
        for key in self.valid_animation_keys:
            in_data = key in self.data.dims
            if in_data:
                key_val = draw_data[key]
                in_draw_data = key in draw_data
                in_sel_call = key in sel_call
                if in_draw_data and not in_sel_call:
                    if self.cumulative:
                        if self.data[key][0].data == key_val:
                            sel_call[key] = key_val
                        else:
                            sel_call[key] = np.arange(key_val)
                    else:
                        sel_call[key] = key_val
        return sel_call

    def process_photons(self, photons):
        if self.cumulative:
            if "frame" in photons.dims:
                # Sum the photons over, if we have not already selected an
                # individual frame (such as the first frame)
                photons = photons.sum("frame")

        photons = photons[self.imsel].data
        return photons
