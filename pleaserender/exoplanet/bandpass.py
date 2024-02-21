import astropy.units as u
import xarray as xr

from pleaserender.core import Plot


class Bandpass(Plot):
    def __init__(self, bandpass, wavelengths=None, obs=None, **kwargs):
        if kwargs.get("axis_keys") is None:
            kwargs["axis_keys"] = {
                "x": "spectral_wavelength(nm)",
                "y": "transmission",
            }
        if "render_state_kwargs" not in kwargs:
            kwargs["render_state_kwargs"] = {
                "key_strategies": {"spectral_wavelength(nm)": "Cumulative"},
                "repeat": True,
            }
        super().__init__(**kwargs)
        self.bandpass = bandpass
        self.plot_method = "plot"
        self.name = "bandpass"
        assert (
            wavelengths is not None or obs is not None
        ), "Either wavelengths or an observation object must be provided."
        if wavelengths is not None:
            self.wavelengths = wavelengths
        elif obs is not None:
            self.wavelengths = obs.spectral_wavelength_grid

    def generate_data(self):
        self.transmission = self.bandpass(self.wavelengths)
        self.data = xr.DataArray(
            self.transmission,
            coords={"spectral_wavelength(nm)": self.wavelengths.to(u.nm).value},
            dims=["spectral_wavelength(nm)"],
        )
        self.data.name = "transmission"
        self.data = self.data.to_dataset()

    def get_required_keys(self):
        self.required_keys = ["spectral_wavelength(nm)"]

    def draw_plot(self):
        # val = self.state.context["spectral_wavelength(nm)"]
        # breakpoint()
        # self.generic_plot(self.data, val, "spectral_wavelength(nm)")
        # data = self.state.context_data()
        self.generic_plot()
