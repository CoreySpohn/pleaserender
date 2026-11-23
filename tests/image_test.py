"""The Image plot builds its data cube from the coronagraphoto functional API."""

from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from optixstuff import Band

from pleaserender.exoplanet import Image, Observation, SpectralCube

NY, NX = 4, 6


class _Rates:
    def __init__(self, value):
        self.value = value

    def total(self):
        return self.value


@pytest.fixture
def instrument():
    grid = SimpleNamespace(shape=(NY, NX))
    detector = SimpleNamespace(grid=grid, pixel_scale_arcsec=0.5)
    return SimpleNamespace(detector=detector)


@pytest.fixture
def stub_rates(monkeypatch):
    """Rates equal to time_d plus the band index, so each cube entry is identifiable."""
    import coronagraphoto as cp

    def build_exposure(scene, instrument, bands, *, time_d, roll_deg, order):
        return SimpleNamespace(time_d=time_d, n_bands=len(bands))

    def expected(scene, instrument, exposure):
        index = np.arange(exposure.n_bands)[:, None, None]
        return _Rates(exposure.time_d + index + np.zeros((exposure.n_bands, NY, NX)))

    monkeypatch.setattr(cp, "build_exposure", build_exposure)
    monkeypatch.setattr(cp, "expected", expected)


BANDS = [Band(500.0, 520.0), Band(600.0, 640.0)]


def test_cube_axes_follow_times_and_bands(instrument, stub_rates):
    """Fails if a time or band lands on the wrong axis or carries the wrong label."""
    plot = Image(None, instrument, BANDS, gen_data={"start_time": [10.0, 20.0, 30.0]})
    plot.generate_data()
    cube = plot.data["img(det)"]
    assert cube.dims == (
        "start_time",
        "spectral_wavelength(nm)",
        "ypix(det)",
        "xpix(det)",
    )
    assert cube.shape == (3, 2, NY, NX)
    np.testing.assert_allclose(cube["spectral_wavelength(nm)"], [510.0, 620.0])
    np.testing.assert_allclose(
        cube.sel(start_time=20.0, **{"spectral_wavelength(nm)": 620.0}), 21.0
    )


def test_single_time_and_band_are_not_axes(instrument, stub_rates):
    """Fails if a length-one time or band stays as an axis imshow cannot draw."""
    plot = Image(None, instrument, BANDS[:1])
    plot.generate_data()
    assert plot.data["img(det)"].dims == ("ypix(det)", "xpix(det)")
    plot.get_required_keys()
    assert plot.required_keys == []


def test_extent_in_pixels_and_arcsec(instrument):
    """Fails if the arcsec extent is not centered with half the field on each side."""
    assert Image(None, instrument, BANDS).get_extent() == [0, NX, 0, NY]
    arcsec = Image(None, instrument, BANDS, imaging_params={"unit": u.arcsec})
    assert arcsec.get_extent() == [-1.5, 1.5, -1.0, 1.0]


def test_subclasses_set_summing_and_title(instrument, stub_rates):
    """Fails if sum_bands or the multi-time title key is not wired."""
    gen = {"start_time": [0.0, 1.0]}
    observation = Observation(None, instrument, BANDS, sum_bands=True, gen_data=gen)
    assert observation.access_kwargs["sum_keys"] == ["spectral_wavelength(nm)"]
    assert observation.animation_kwargs["title_key"] == "start_time"
    cube = SpectralCube(None, instrument, BANDS, cumulative=True)
    assert cube.access_kwargs["sum_keys"] == ["spectral_wavelength(nm)"]
