"""Image matching, and the self-calibration built on top of it.

The most important test in this file is the sign convention of the matcher. Getting a
displacement's sign backwards does not produce a weaker correction — it produces one
that **doubles** the error it was meant to remove, and every internal check still
passes, because the arithmetic is self-consistent about the wrong thing. It happened
once during development, and only the end-to-end residual caught it.
"""

from __future__ import annotations

import numpy as np
import pytest

from spp.refine.matcher import (
    Match,
    gradient_magnitude,
    has_texture,
    match_windows,
    phase_correlate,
    robust_median,
)


def _textured(size: int = 128, seed: int = 0) -> np.ndarray:
    """An image with structure at several scales, as a real scene has."""
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:size, 0:size]
    return (
        100.0
        + 30.0 * np.sin(x / 7.0) * np.cos(y / 11.0)
        + 20.0 * ((x // 16 + y // 16) % 2)
        + rng.normal(0.0, 2.0, (size, size))
    )


# -- the sign convention ----------------------------------------------------


def test_the_sign_convention_is_what_the_callers_assume():
    """`moving[l, c]` shows what `reference[l + dline, c + dcolumn]` shows.

    This is the one test that must never be quietly "fixed" by flipping an
    expectation. The self-calibration inverts this displacement through the sensor
    model; read it backwards and the estimated correction has the wrong sign, so
    applying it doubles the misregistration instead of removing it — while every
    intermediate number still looks reasonable.
    """
    reference = _textured()
    # Translate the content by +5 rows and +3 columns.
    moving = np.roll(np.roll(reference, 5, axis=0), 3, axis=1)

    dline, dcolumn, peak = phase_correlate(reference, moving)

    assert peak > 0.5
    assert dline == pytest.approx(-5.0, abs=0.2)
    assert dcolumn == pytest.approx(-3.0, abs=0.2)

    # And the property the callers actually rely on, stated directly: sampling the
    # reference at (l + dline, c + dcolumn) recovers what the moving image shows.
    # This is the line the self-calibration depends on.
    line, column = 40, 60
    assert moving[line, column] == pytest.approx(
        reference[round(line + dline), round(column + dcolumn)]
    )


def test_a_zero_shift_measures_zero():
    dline, dcolumn, peak = phase_correlate(_textured(), _textured())
    assert peak > 0.9
    assert dline == pytest.approx(0.0, abs=0.05)
    assert dcolumn == pytest.approx(0.0, abs=0.05)


def test_sub_pixel_shifts_are_resolved():
    """A whole-pixel-only matcher would be useless: the corrections are sub-pixel.

    The test scene is deliberately *aperiodic*. Phase correlation on a periodic pattern
    is genuinely ambiguous — it can lock onto any multiple of the period, and would
    report a large shift with a confident peak — so a periodic test image would be
    testing the wrong thing.
    """
    size = 128
    y, x = np.mgrid[0:size, 0:size].astype(float)
    # A single smooth blob: one unambiguous feature, shifted by a known fraction.
    def blob(dy: float, dx: float) -> np.ndarray:
        return np.exp(-(((x - 64 - dx) / 18.0) ** 2 + ((y - 64 - dy) / 22.0) ** 2))

    dline, dcolumn, _ = phase_correlate(blob(0.0, 0.0), blob(1.6, 0.4))
    assert dline == pytest.approx(-1.6, abs=0.25)
    assert dcolumn == pytest.approx(-0.4, abs=0.25)


# -- matching across spectral bands -----------------------------------------


def test_matching_survives_an_inverted_contrast():
    """Two bands can disagree about a scene's *sign* and still agree about its edges.

    Vegetation is dark in blue and bright in the near infrared; water is the reverse.
    Correlating raw radiance across bands would fight that inversion. Gradients do not
    care: an edge is an edge whichever way the contrast runs, which is why the matcher
    works on gradient magnitude.
    """
    reference = _textured()
    inverted = 200.0 - reference  # same structure, opposite polarity
    moving = np.roll(np.roll(inverted, 4, axis=0), -2, axis=1)

    dline, dcolumn, peak = phase_correlate(reference, moving)

    assert peak > 0.3
    assert dline == pytest.approx(-4.0, abs=0.3)
    assert dcolumn == pytest.approx(2.0, abs=0.3)


def test_gradient_magnitude_ignores_a_constant_offset():
    """A brightness difference between bands must not register as structure."""
    image = _textured()
    np.testing.assert_allclose(
        gradient_magnitude(image), gradient_magnitude(image + 500.0), atol=1e-9
    )


# -- refusing to measure what is not there ----------------------------------


def test_featureless_water_is_rejected():
    """Phase correlation always returns a peak. Over flat water, that peak is noise.

    A matcher that reports a displacement for a featureless window is not measuring;
    it is generating. The texture gate is what keeps "no evidence" from being averaged
    in as "zero shift".
    """
    rng = np.random.default_rng(1)
    flat_water = 50.0 + rng.normal(0.0, 0.05, (128, 128))  # ~0.1% contrast
    assert not has_texture(flat_water)
    assert has_texture(_textured())


def test_windows_without_texture_are_not_matched():
    """The gate must actually keep them out of the result, not merely flag them."""
    rng = np.random.default_rng(2)
    reference = np.full((768, 768), 50.0) + rng.normal(0.0, 0.05, (768, 768))
    assert match_windows(reference, reference.copy(), window=384) == []


def test_no_matches_means_no_estimate_not_a_zero_shift():
    """An empty match set is 'unknown', and must never collapse to 'no displacement'."""
    displacement, _, count = robust_median([])
    assert count == 0
    assert np.isnan(displacement).all()


# -- robustness -------------------------------------------------------------


def test_outliers_do_not_move_the_consensus():
    """A ship, a cloud edge or a glint patch will mismatch. The answer must survive it.

    A mean would carry the outliers into the estimate; the median plus a
    median-absolute-deviation gate does not.
    """
    good = [Match(0, 0, 2.0 + 0.05 * i, -1.0, 0.8) for i in range(20)]
    wild = [Match(0, 0, 40.0, 60.0, 0.3), Match(0, 0, -35.0, 25.0, 0.3)]

    displacement, _, used = robust_median(good + wild)

    assert displacement[0] == pytest.approx(2.5, abs=0.3)
    assert displacement[1] == pytest.approx(-1.0, abs=0.1)
    assert used == len(good)


def test_mismatched_shapes_are_an_error():
    with pytest.raises(ValueError, match="same shape"):
        match_windows(np.zeros((64, 64)), np.zeros((64, 32)))
