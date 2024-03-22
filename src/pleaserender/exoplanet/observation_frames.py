from pleaserender.exoplanet import Image


class ObservationFrames(Image):
    def __init__(self, observation, cumulative=False, *args, **kwargs):
        all_keys = []
        sum_keys = []
        if observation.settings.return_spectrum:
            all_keys.append("spectral_wavelength(nm)")
            sum_keys.append("spectral_wavelength(nm)")

        cumulative_keys = []
        if cumulative:
            cumulative_keys.append("time")
            sum_keys.append("time")

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

        # Set title key
        if "animation_kwargs" not in kwargs:
            kwargs["animation_kwargs"] = {"title_key": "time"}
        else:
            if "title_key" not in kwargs["animation_kwargs"]:
                kwargs["animation_kwargs"]["title_key"] = "time"

        super().__init__(observation, *args, **kwargs)
        self.observation = observation
