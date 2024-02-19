import copy

import astropy.units as u
import numpy as np
from astropy.time import Time

from pleaserender.exoplanet import Image


class ObservationFrames(Image):
    def __init__(self, observation, cumulative_keys=[], all_keys=[], *args, **kwargs):
        if "render_state_kwargs" not in kwargs:
            kwargs["render_state_kwargs"] = {"key_strategies": {}}
        if "access_kwargs" not in kwargs:
            kwargs["access_kwargs"] = {"sum_keys": []}
        if all_keys:
            for all_key in all_keys:
                kwargs["render_state_kwargs"]["key_strategies"][all_key] = "All"
                kwargs["access_kwargs"]["sum_keys"].append(all_key)
        if cumulative_keys:
            for cum_key in cumulative_keys:
                kwargs["render_state_kwargs"]["key_strategies"][cum_key] = "Cumulative"
                kwargs["access_kwargs"]["sum_keys"].append(cum_key)
        # else:
        #     required_keys = self.get_required_keys()
        #     for key in required_keys:
        #         kwargs["render_state_kwargs"]["key_strategies"][key] = "Value"

        super().__init__(observation, *args, **kwargs)
        self.observation = observation

    # def get_frame_data(self):
    #     # sel_call = copy.copy(self.state.context)
    #     # for key, val in sel_call.items():
    #     #     if key in self.cumulative_keys:
    #     #         sel_call[key] = slice(None, val)
    #     #     if key in self.sum_keys:
    #     #         sel_call[key] = self.data.coords[key].values
    #     # base_data = self.data.sel(**sel_call)
    #     base_data = self.state.context_data()
    #     photons = self.access_data(base_data)
    #     # photons = self.process_photons(base_data)
    #     return photons

    def process_photons(self, photons):
        for key in self.cumulative_keys:
            photons = photons.sum(key)
        for key in self.sum_keys:
            photons = photons.sum(key)

        photons = photons[self.imsel].data
        return photons
