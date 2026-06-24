# References

Standards, conventions and tooling consulted while building this framework.

## Standards & product definitions

- **CEOS — Committee on Earth Observation Satellites.** Conventions for EO
  product levels, traceable processing chains, and calibration/validation
  practice. <https://ceos.org>
- **STAC — SpatioTemporal Asset Catalog.** Product and catalogue metadata model;
  the input package's item is a STAC 1.1 item. <https://stacspec.org>

## Data handling & raster I/O

- **Rasterio** — windowed raster I/O, overviews, GeoTIFF/COG writing.
  <https://rasterio.readthedocs.io>
- **GDAL / OGR** — underlying raster/vector engine; warp, VRT, COG.
  <https://gdal.org>
- **NumPy** — array computation for the calibration chain.
  <https://numpy.org>
- **Cloud-Optimized GeoTIFF (COG)** — tiled + overview raster layout.
  <https://www.cogeo.org>

## Geometric reference (L1C)

- **Copernicus DEM** — elevation reference for terrain correction.
  <https://spacedata.copernicus.eu>
- **GDAL warp / VRT** — resampling onto a target CRS/grid.
  <https://gdal.org>

## Solar spectral reference (radiometric plausibility, L2A ESUN)

The package solar reference (`spectral_solar.json`) is a **TSIS**-style spectrum;
these underpin the expected radiance magnitude and the visible→NIR ordering.

- **TSIS-1 Hybrid Solar Reference Spectrum** — Coddington et al. (2021),
  *Geophysical Research Letters*. <https://lasp.colorado.edu/lisird/data/tsis1_hsrs>
- **Thuillier et al. (2003)** — solar spectral irradiance (SOLSPEC), *Solar Physics* 214.
- **ASTM E490** — standard extraterrestrial solar spectral irradiance.

## Atmospheric reference (L2A)

- **6S / Py6S** — atmospheric radiative-transfer model and Python interface.
  <https://py6s.readthedocs.io>, <https://6s.ltdri.org>
- **ACOLITE** — atmospheric correction for aquatic scenes.
  <https://github.com/acolite>

## Visualisation & QA

- **Matplotlib** — RGB quicklook output. <https://matplotlib.org>
- **QGIS** — visual inspection of rasters. <https://qgis.org>

## Acronyms

Acronyms used across the documentation, in one place.

| Acronym | Expansion |
|---|---|
| ADC | Analogue-to-Digital Converter |
| AOD | Aerosol Optical Depth |
| CEOS | Committee on Earth Observation Satellites |
| CIS | CMOS Image Sensor |
| COG | Cloud-Optimized GeoTIFF |
| CPF | Calibration Parameter File |
| CRC | Cyclic Redundancy Check |
| CRS | Coordinate Reference System |
| CWL | Central Wavelength |
| DEM | Digital Elevation Model |
| DN | Digital Number (raw detector count) |
| ECEF | Earth-Centred, Earth-Fixed (coordinate frame) |
| EO | Earth Observation / Electro-Optical (STAC extension) |
| ESUN | Exo-atmospheric Solar Irradiance (band-mean solar irradiance) |
| GCP | Ground Control Point |
| GNSS | Global Navigation Satellite System |
| GSD | Ground Sample Distance |
| LoS | Line of Sight |
| NIR | Near Infrared |
| PAN | Panchromatic |
| RPC | Rational Polynomial Coefficients |
| RSR | Relative Spectral Response |
| RT | Radiative Transfer |
| SLERP | Spherical Linear Interpolation (of quaternions) |
| STAC | SpatioTemporal Asset Catalog |
| SWIR | Short-Wave Infrared |
| TDI | Time-Delay Integration |
| TOA | Top of Atmosphere |
| VNIR | Visible and Near Infrared |
