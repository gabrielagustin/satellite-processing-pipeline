# L1C — Architecture

L1B is radiometric and streams in windows; L1C is **geometric** and works from a
sensor model. The design is in [`l1c_spec.md`](spec.md); this is the shape of the
code.

```text
L1B radiance (sensor coords)  +  ephemeris/attitude  +  DEM  +  calibration
        │
        ▼
   SensorModel.locate(band, lines, columns, terrain)  ──►  (lon, lat, height)
        │      composes: timing · ephemeris · camera · frames · terrain
        │      parameterised by: Convention  (how to read the telemetry)
        ▼
   geoloc_grid.build(model, band)  ──►  GeolocationGrid   (lattice, one per band)
        │
        ▼
   TargetGrid.covering(footprints)  ──►  one projected grid for all bands
        │
        ▼
   GeolocWarper.warp(band, geoloc, target)  ──►  L1C raster
```

### Layers

```text
src/spp/geometry/
  frames.py         WGS84, quaternion algebra, orbital frame, Earth rotation,
                    ray-ellipsoid intersection, geodetic conversion
  timing.py         imager ticks -> UTC; per-line exposure times; gap detection
  ephemeris.py      Ephemeris (Hermite) and Attitude (SLERP, rate-checked)
  camera.py         interior orientation: detector column -> body-frame ray
  terrain.py        TerrainModel seam; ellipsoid and DEM surfaces; geoid; the
                    iterative ray-terrain intersection
  sensor_model.py   composes the above; carries the Convention
  conventions.py    resolves the telemetry's undocumented conventions by measurement
  geoloc_grid.py    the sensor model, sampled onto a lattice
src/spp/resample/
  grid.py           the target grid: CRS, ground sample distance, extent
  warper.py         Warper seam; the geolocation-array warp
```

### The seams

```python
class TerrainModel(ABC):
    def height(self, lon, lat) -> np.ndarray: ...
    def intersect(self, origins, directions) -> tuple[np.ndarray, np.ndarray]: ...

class Warper(ABC):
    def warp(self, source_path, geoloc, target, output_path) -> Path: ...
```

`SensorModel` is concrete but *parameterised*: everything the delivered package leaves
undocumented — which frame the ephemeris is in, which frame and direction and component
order the attitude quaternion uses, how the detector's axes map into the camera frame,
which way the scan ran — lives in a `Convention` dataclass rather than in the code. It
is resolved once, by measurement, and the model states which reading of the telemetry
it depends on.

### Two properties worth calling out

**Each band is located independently.** Its own line times, its own view angle, its own
terrain intersection, its own geolocation grid. They are then resampled onto **one**
target grid. Co-registration is not a correction applied to the bands; it is what
happens when each is put where it belongs.

**The geolocation lattice tracks the DEM, not the orbit.** The obvious optimisation — a
coarse lattice, since the orbit and attitude are smooth — would low-pass the terrain
component of the geolocation field and silently undo part of the orthorectification. The
spacing is 8 detector samples (~30 m, the elevation model's own resolution), not 64.
