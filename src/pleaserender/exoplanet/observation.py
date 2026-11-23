from pleaserender.exoplanet import Image


class Observation(Image):
    """Detector images over a set of start times, summed over the bands when asked."""

    def __init__(self, scene, instrument, bands, *args, sum_bands=False, **kwargs):
        if "render_state_kwargs" not in kwargs:
            kwargs["render_state_kwargs"] = {"key_strategies": {}}
        if "access_kwargs" not in kwargs:
            kwargs["access_kwargs"] = {"sum_keys": []}
        if sum_bands:
            kwargs["render_state_kwargs"]["key_strategies"][
                "spectral_wavelength(nm)"
            ] = "All"
            kwargs["access_kwargs"]["sum_keys"].append("spectral_wavelength(nm)")

        # Set title key
        gen_data = kwargs.get("gen_data") or {}
        if len(gen_data.get("start_time", [0.0])) > 1:
            if "animation_kwargs" not in kwargs:
                kwargs["animation_kwargs"] = {"title_key": "start_time"}
            else:
                if "title_key" not in kwargs["animation_kwargs"]:
                    kwargs["animation_kwargs"]["title_key"] = "start_time"

        super().__init__(scene, instrument, bands, *args, **kwargs)
