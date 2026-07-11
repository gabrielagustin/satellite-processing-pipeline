"""Terrain intersection: orthorectification, and the height datum.

Runs entirely on synthetic elevation, so it needs no network and no external data.
The physics being checked is exact and known in closed form — relief displaces a
pixel by `height x tan(view angle)` — which makes it testable without a real DEM.
"""

from __future__ import annotations

import numpy as np
import pytest

from spp.geometry.frames import ecef_to_geodetic, geodetic_to_ecef
from spp.geometry.terrain import (
    DEMTerrain,
    EllipsoidTerrain,
    GeoidModel,
    check_undulation,
)

ALTITUDE_M = 400_000.0


def _flat_geoid(undulation_m: float = -30.0) -> GeoidModel:
    """A constant-undulation geoid — enough to test the datum handling."""
    lons = np.linspace(-1.0, 1.0, 5)
    lats = np.linspace(-1.0, 1.0, 5)
    return GeoidModel(
        lons=lons,
        lats=lats,
        undulation=np.full((lats.size, lons.size), undulation_m),
    )


def _plateau_dem(height_m: float, geoid: GeoidModel | None = None) -> DEMTerrain:
    """A DEM of uniform orthometric height, over a small patch about (0, 0)."""
    lons = np.linspace(-1.0, 1.0, 201)
    lats = np.linspace(1.0, -1.0, 201)  # descending, as a north-up raster is
    return DEMTerrain(
        elevation=np.full((lats.size, lons.size), float(height_m)),
        lons=lons,
        lats=lats,
        geoid=geoid or _flat_geoid(0.0),
    )


def _ray_from(lon: float, lat: float, view_angle_rad: float) -> tuple[np.ndarray, np.ndarray]:
    """A ray from directly above (lon, lat), tilted `view_angle_rad` towards the east."""
    platform = geodetic_to_ecef(np.array([[lon, lat, ALTITUDE_M]]))

    down = -platform / np.linalg.norm(platform, axis=-1, keepdims=True)
    # A local east vector, orthogonal to the radius.
    east = np.cross(np.array([0.0, 0.0, 1.0]), platform)
    east /= np.linalg.norm(east, axis=-1, keepdims=True)

    direction = down * np.cos(view_angle_rad) + east * np.sin(view_angle_rad)
    return platform, direction / np.linalg.norm(direction, axis=-1, keepdims=True)


# -- the geoid guard --------------------------------------------------------


def test_a_real_undulation_is_accepted():
    check_undulation(np.array([[-31.0, -30.5], [-30.0, -29.5]]))


def test_an_all_zero_undulation_is_rejected():
    """The silent failure this module exists to prevent.

    Without its geoid grid, PROJ returns the height unchanged rather than raising.
    Every elevation would then be wrong by the local undulation — tens of metres,
    uniformly, with no other symptom. It must not be possible to build a geoid
    model out of that.
    """
    with pytest.raises(RuntimeError, match="ballpark"):
        check_undulation(np.zeros((4, 4)))


def test_a_non_finite_undulation_is_rejected():
    with pytest.raises(RuntimeError, match="non-finite"):
        check_undulation(np.array([[-30.0, np.nan]]))


# -- sampling ---------------------------------------------------------------


def test_dem_height_is_orthometric_plus_undulation():
    """The datum conversion, which is the whole reason the geoid is here.

    A DEM cell reading 100 m over a geoid 30 m *below* the ellipsoid is at an
    ellipsoidal height of 70 m — not 100. Getting this backwards, or skipping it,
    is a systematic error the size of the undulation.
    """
    terrain = _plateau_dem(100.0, geoid=_flat_geoid(-30.0))
    height = terrain.height(np.array([0.1]), np.array([0.2]))
    np.testing.assert_allclose(height, [70.0], atol=1e-9)


def test_sampling_outside_the_grid_clamps_rather_than_extrapolates():
    """A ray landing just off the loaded window gets the nearest real height.

    Extrapolating a DEM invents terrain, and the invented value feeds straight back
    into the intersection. Clamping is wrong by the edge gradient; extrapolation can
    be wrong without bound.
    """
    terrain = _plateau_dem(250.0)
    inside = terrain.height(np.array([0.0]), np.array([0.0]))
    outside = terrain.height(np.array([50.0]), np.array([-40.0]))
    np.testing.assert_allclose(outside, inside, atol=1e-9)


# -- intersection -----------------------------------------------------------


def test_ellipsoid_terrain_lands_at_its_own_height():
    terrain = EllipsoidTerrain(height_m=0.0)
    origin, direction = _ray_from(10.0, 20.0, 0.0)

    point, converged = terrain.intersect(origin, direction)
    assert converged.all()
    np.testing.assert_allclose(ecef_to_geodetic(point)[0, 2], 0.0, atol=0.1)


def test_intersection_lands_on_the_terrain_surface():
    """The ground point must sit *on* the DEM, not on the ellipsoid beneath it."""
    terrain = _plateau_dem(1000.0)
    origin, direction = _ray_from(0.0, 0.0, np.radians(1.0))

    point, converged = terrain.intersect(origin, direction)
    assert converged.all()

    geodetic = ecef_to_geodetic(point)[0]
    np.testing.assert_allclose(geodetic[2], 1000.0, atol=1.0)


def test_relief_displacement_matches_the_closed_form():
    """Terrain shifts a pixel by `height x tan(incidence angle)` — check it, don't assume.

    This is the quantity orthorectification exists to remove, so a model that gets it
    wrong is not orthorectifying, it is only moving pixels around plausibly.

    The angle in that formula is the incidence angle **at the ground**, not the view
    angle at the platform, and the two differ: Earth's curvature tilts the local
    vertical between the two, magnifying the angle by `(R + altitude) / R`. At 400 km
    that is a factor of 1.063 — a 6% effect, small but well above the accuracy this
    level claims, and the reason a flat-Earth `h x tan(view angle)` underestimates the
    displacement.
    """
    view_angle = np.radians(1.0)
    height = 2000.0
    earth_radius = 6_371_008.8

    origin, direction = _ray_from(0.0, 0.0, view_angle)
    flat = ecef_to_geodetic(EllipsoidTerrain(0.0).intersect(origin, direction)[0])[0]
    raised = ecef_to_geodetic(_plateau_dem(height).intersect(origin, direction)[0])[0]

    # The *horizontal* shift is the relief displacement. The 3D separation between the
    # two points is dominated by the 2 km of height between them, which is not what
    # orthorectification corrects.
    metres_per_degree = 111_320.0
    displacement = metres_per_degree * np.hypot(
        (raised[0] - flat[0]) * np.cos(np.radians(flat[1])), raised[1] - flat[1]
    )

    magnification = (earth_radius + ALTITUDE_M) / earth_radius
    incidence = np.arcsin(magnification * np.sin(view_angle))
    expected = height * np.tan(incidence)

    assert abs(displacement - expected) / expected < 0.02


def test_a_nadir_ray_is_untouched_by_terrain():
    """Straight down, relief changes the height but not the ground position.

    The counterpart to the test above: `tan(0) = 0`. If terrain moved a nadir pixel
    sideways, the intersection geometry would be wrong in a way the tilted case
    might hide.
    """
    origin, direction = _ray_from(0.0, 0.0, 0.0)
    flat = ecef_to_geodetic(EllipsoidTerrain(0.0).intersect(origin, direction)[0])[0]
    raised = ecef_to_geodetic(_plateau_dem(2000.0).intersect(origin, direction)[0])[0]

    np.testing.assert_allclose(raised[:2], flat[:2], atol=1e-6)
    assert abs(raised[2] - 2000.0) < 1.0


def test_intersection_converges_quickly():
    """Two iterations should be enough for a narrow-field instrument.

    Pinning this keeps the cost honest: the terrain loop runs once per geolocation
    node, and a model that needed many iterations would be a different performance
    proposition than the one the spec assumes.
    """
    terrain = _plateau_dem(1500.0)
    origin, direction = _ray_from(0.0, 0.0, np.radians(1.1))

    _, converged = terrain.intersect(origin, direction, max_iterations=2)
    assert converged.all()


def test_a_ray_that_misses_the_earth_does_not_converge():
    """Pointing at the sky must not yield a ground point, however many iterations."""
    origin = np.array([[10_000e3, 0.0, 0.0]])
    direction = np.array([[1.0, 0.0, 0.0]])

    point, converged = EllipsoidTerrain(0.0).intersect(origin, direction)
    assert not converged.any()
    assert np.isnan(point).all()
