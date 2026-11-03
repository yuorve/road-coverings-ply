# Road Influence Zones PLY Generator

A complete Python and Flask-based web interface to parse road networks, calibrate kilometer points (PKs), handle physical road discontinuities, resolve geometric overlaps at intersections, and generate standard PLY files (GeoMedia-compatible) representing road influence buffers.

## Features
- **Robust Data Parsing**: Seamlessly processes GeoMedia `.mdb` databases using native `mdbtools` or GDAL/OGR (ODBC), as well as standard Shapefiles via `pyshp`.
- **Discontinuity Handling**: Detects and correctly splits routes that are physically disconnected or share bypass segments.
- **Advanced Geoprocessing**: Dynamically splits dual-carriageway buffers when lanes diverge and handles intersection/bifurcation overlaps with a customizable hierarchy.
- **Interactive GIS Map Preview**: Built-in interactive SVG-based web map utilizing panning, smooth zooming, boundary overlays (municipalities and provinces), and dynamic tooltips pointing directly to geometry deviations.
- **CSV & PLY Exporting**: Easily import source tables via CSV and export fully optimized GeoMedia-formatted `.ply` and summary `.csv` deliverables.

## Requirements & Dependencies

### System Packages
- `mdbtools` (specifically `mdb-export` and `mdb-tables` for fast, GDAL-free access to Access databases)
- `gdal-bin` & `libgdal-dev` (optional, needed only if reading directly from Access ODBC/OGR structures)

### Python Libraries
- `Flask` (Web service)
- `shapely>=2` (Computational geometry)
- `pyshp` (Shapefile parsing)
- `numpy` (Numerical operations)
- `fiona` (Optional, GIS formats)

## Installation

```bash
# Install system dependencies (Ubuntu/Debian example)
sudo apt-get update && sudo apt-get install -y mdbtools libgdal-dev gdal-bin

# Install python packages
pip install flask pyshp "shapely>=2" fiona

Running the Application
To start the Flask development server:

python app_web.py
It will launch the interface at http://127.0.0.1:5000 and automatically attempt to open your system's web browser.

File Structure
app_web.py: The Flask server defining REST API endpoints (/api/proyecto, /api/previsualizar, /api/conectar, /api/exportar/ply, etc.) and overall project state management.
fuentes.py: Low-level data connectors for extracting geometries from axes and kilometer marks. Features customized bin decoding for GeoMedia GDO formats.
ply_zonas.py: Core computational geometry and GIS logic. Performs calibration, linear reference projections, topology correction, overlap resolution, and file exporting.
index.html: Responsive single-page application built with modern vanilla JavaScript and reactive SVG mechanics for real-time visualization.
