import astropy.units as u
import coronagraphoto as cp
import numpy as np
import xarray as xr
from matplotlib.colors import Normalize

from pleaserender import util
from pleaserender.core import Plot


class Image(Plot):
    """Detector images of a scene, rendered by ``coronagraphoto.expected``.

    Args:
        scene: The ``skyscapes.Scene`` to image.
        instrument: The ``optixstuff.Instrument`` to image it with.
        bands: The ``optixstuff.Band`` list; each band is one entry along the
            ``spectral_wavelength(nm)`` axis, labeled by its center.
        gen_data: Optional dict with ``start_time``, the observation times in
            days from the scene epoch. One time is used when omitted.
        imaging_params: Overrides of ``object`` (the data variable prefix),
            ``plane`` (only ``det`` exists) and ``unit`` (``u.pix`` or
            ``u.arcsec`` for the axes).
        roll_deg: Position angle of detector +y, north through east.
        order: Quadrature nodes per band piece.
    """

    def __init__(
        self,
        scene,
        instrument,
        bands,
        gen_data=None,
        imaging_params=None,
        roll_deg=0.0,
        order=2,
        **kwargs,
    ):
        default_imaging_params = {"object": "img", "plane": "det", "unit": u.pix}
        self.imaging_params = default_imaging_params
        if imaging_params is not None:
            self.imaging_params.update(imaging_params)
        self.imsel = f"{self.imaging_params['object']}({self.imaging_params['plane']})"

        if kwargs.get("axis_keys") is None:
            kwargs["axis_keys"] = {
                "x": f"xpix({self.imaging_params['plane']})",
                "y": f"ypix({self.imaging_params['plane']})",
            }
        if not kwargs.get("ax_kwargs"):
            kwargs["ax_kwargs"] = {}
        if "auto_title" not in kwargs["ax_kwargs"]:
            kwargs["ax_kwargs"]["auto_title"] = True
        if "aspect" not in kwargs["ax_kwargs"]:
            kwargs["ax_kwargs"]["aspect"] = "equal"

        # Set the default plot_kwargs
        if not kwargs.get("plot_kwargs"):
            kwargs["plot_kwargs"] = {}
        if "cmap" not in kwargs["plot_kwargs"]:
            kwargs["plot_kwargs"]["cmap"] = "viridis"
        if "norm" not in kwargs["plot_kwargs"]:
            kwargs["plot_kwargs"]["norm"] = Normalize(vmin=0, vmax=4e6)
        if "origin" not in kwargs["plot_kwargs"]:
            kwargs["plot_kwargs"]["origin"] = "lower"

        super().__init__(**kwargs)

        self.gen_data = gen_data
        self.scene = scene
        self.instrument = instrument
        self.bands = list(bands)
        self.roll_deg = roll_deg
        self.order = order
        self.plot_method = "imshow"
        self.name = "image"

    def generate_data(self):
        gen_data = self.gen_data or {}
        times = np.atleast_1d(
            np.asarray(gen_data.get("start_time", [0.0]), dtype=float)
        )
        centers = [band.center_nm for band in self.bands]
        frames = []
        for time_d in times:
            exposure = cp.build_exposure(
                self.scene,
                self.instrument,
                self.bands,
                time_d=float(time_d),
                roll_deg=self.roll_deg,
                order=self.order,
            )
            rates = cp.expected(self.scene, self.instrument, exposure)
            frames.append(np.asarray(rates.total()))
        ny, nx = frames[0].shape[-2:]
        plane = self.imaging_params["plane"]
        dims = [
            "start_time",
            "spectral_wavelength(nm)",
            f"ypix({plane})",
            f"xpix({plane})",
        ]
        cube = xr.DataArray(
            np.stack(frames),
            dims=dims,
            coords={
                "start_time": times,
                "spectral_wavelength(nm)": centers,
                f"ypix({plane})": np.arange(ny),
                f"xpix({plane})": np.arange(nx),
            },
            name=self.imsel,
        )
        # A single time or band is not an axis to animate over.
        for dim in ("start_time", "spectral_wavelength(nm)"):
            if cube.sizes[dim] == 1:
                cube = cube.squeeze(dim)
        self.data = cube.to_dataset()

    def get_required_keys(self):
        """
        Get the keys that are required to plot the data
        """
        all_possible_keys = ["start_time", "spectral_wavelength(nm)"]
        for key in all_possible_keys:
            if key in self.data.dims:
                self.required_keys.append(key)
        self.skipna = False

    def draw_plot(self, plot_kwargs=None):
        # photons = self.get_frame_data()
        if plot_kwargs is None:
            plot_kwargs = self.plot_kwargs
        data = self.get_plot_data()
        extent = self.get_extent()
        self.ax.imshow(data, extent=extent, **plot_kwargs)

    def get_plot_data(self, context=None):
        """
        Method to get the data required to plot. Not implemented in the core classes.
        """
        sel = self.state.context_sel(context)
        plot_data = self.data.sel(**sel)

        # Sum over keys if necessary
        for key in self.access_kwargs["sum_keys"]:
            if key in plot_data.dims and len(plot_data[key].shape) > 0:
                plot_data = plot_data.sum(dim=key, skipna=self.skipna)
        data = plot_data[self.imsel].data
        return data

    def check_valid_context(self, fig_context):
        context = self.state.extract_plot_context(fig_context)
        if context is None:
            return False
        for key, val in context.items():
            if val not in self.data[key]:
                return False
        plot_data = self.get_plot_data(context)
        return bool(np.any(~np.isnan(plot_data)))

    def access_data(self, data):
        return data[self.imsel].data

    def ax_lims_helper(self, necessary_axes, data=None, equal=False):
        extent = self.get_extent()
        xrange = [extent[0], extent[1]]
        yrange = [extent[2], extent[3]]
        self.ax_kwargs["lims"]["x"] = util.calculate_axis_limits_and_ticks(*xrange)
        self.ax_kwargs["lims"]["y"] = util.calculate_axis_limits_and_ticks(*yrange)

    def get_extent(self):
        ny, nx = self.instrument.detector.grid.shape
        if self.imaging_params["unit"] == u.pix:
            return [0, nx, 0, ny]
        pixel = self.instrument.detector.pixel_scale_arcsec * u.arcsec
        unit = self.imaging_params["unit"]
        half_x = (nx / 2 * pixel).to(unit).value
        half_y = (ny / 2 * pixel).to(unit).value
        return [-half_x, half_x, -half_y, half_y]
