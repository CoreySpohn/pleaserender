from pleaserender.exoplanet import Image


class SpectralCube(Image):
    def __init__(self, observation, cumulative=False, *args, **kwargs):
        if "render_state_kwargs" not in kwargs:
            kwargs["render_state_kwargs"] = {"key_strategies": {}}
        if "access_kwargs" not in kwargs:
            kwargs["access_kwargs"] = {"sum_keys": []}

        # Set key strategies
        if cumulative:
            kwargs["render_state_kwargs"]["key_strategies"][
                "spectral_wavelength(nm)"
            ] = "Cumulative"
            kwargs["access_kwargs"]["sum_keys"].append("spectral_wavelength(nm)")
        else:
            kwargs["render_state_kwargs"]["key_strategies"][
                "spectral_wavelength(nm)"
            ] = "Value"

        # Set title key
        if "animation_kwargs" not in kwargs:
            kwargs["animation_kwargs"] = {"title_key": "spectral_wavelength(nm)"}
        else:
            if "title_key" not in kwargs["animation_kwargs"]:
                kwargs["animation_kwargs"]["title_key"] = "spectral_wavelength(nm)"

        # Initialize the Image class
        super().__init__(observation, *args, **kwargs)
        self.observation = observation
