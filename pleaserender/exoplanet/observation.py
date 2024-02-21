import numpy as np
from astropy.time import Time

from pleaserender.exoplanet import Image


class Observation(Image):
    def __init__(self, observation, *args, **kwargs):
        all_keys = []
        sum_keys = []
        if observation.return_spectrum:
            all_keys.append("spectral_wavelength(nm)")
            sum_keys.append("spectral_wavelength(nm)")

        if observation.return_frames:
            all_keys.append("time")
            sum_keys.append("time")

        if "render_state_kwargs" not in kwargs:
            kwargs["render_state_kwargs"] = {"key_strategies": {}}
        if "access_kwargs" not in kwargs:
            kwargs["access_kwargs"] = {"sum_keys": []}
        if all_keys:
            for all_key in all_keys:
                kwargs["render_state_kwargs"]["key_strategies"][all_key] = "All"
                kwargs["access_kwargs"]["sum_keys"].append(all_key)

        # Set title key
        multiple_start_times = len(kwargs["gen_data"]["start_time"]) > 1
        if multiple_start_times:
            if "animation_kwargs" not in kwargs:
                kwargs["animation_kwargs"] = {"title_key": "start_time"}
            else:
                if "title_key" not in kwargs["animation_kwargs"]:
                    kwargs["animation_kwargs"]["title_key"] = "start_time"

        super().__init__(observation, *args, **kwargs)
        self.observation = observation
